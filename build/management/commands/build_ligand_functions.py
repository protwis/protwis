#from django.db import connection
#from django.db import IntegrityError
from django.db.models import Q, Count, Max
from django.utils.text import slugify
from django.db import transaction, IntegrityError

from ligand.models import Ligand, LigandID, LigandType, LigandRole
from common.models import WebResource
from common.tools import get_or_create_url_cache, fetch_from_web_api, save_to_cache, fetch_from_cache
from collections import defaultdict, Counter

import time
import os
import re
import requests
import xmltodict
import logging
import contextlib

import pandas as pd
import datamol as dm
import pubchempy as pcp
from pubchempy import BadRequestError, PubChemHTTPError
from rdkit import RDLogger
from rdkit import Chem
from chembl_structure_pipeline import standardizer
from typing import Dict, Any, Optional

from protein.models import Protein
from .ligand_normalization import normalize_for_parent, normalize_for_child, infer_stereo_status

# Disable the RDkit verbosity
RDLogger.DisableLog('rdApp.*')

log = logging.getLogger(__name__)

external_sources = ["pubchem", "gtoplig", "chembl_ligand", "chembl_peptide", "drugbank", "drug_central"]

_LIGAND_TYPE_CACHE: dict = {}
_WEB_RESOURCE_CACHE: dict = {}
_PDBE_MAX_LEN = Ligand._meta.get_field("pdbe").max_length
APO_NO_LIGAND_NAME = "Apo (no ligand)"

# Set once per worker process (see build_structures.main_func) to the shared
# multiprocessing.Lock used to serialize the DB-touching part of ligand
# resolution across parallel build workers. Left None (no locking) outside
# a parallel build context, e.g. single-process call sites/tests.
_ligand_lock = None


def set_ligand_lock(lock):
    global _ligand_lock
    _ligand_lock = lock


def _ligand_lock_cm():
    return _ligand_lock if _ligand_lock is not None else contextlib.nullcontext()

_GTP_CACHE_TTL = 7 * 24 * 3600  # 7 days

# Curated GtP ligands to force-include in the build even when GtP's own
# interactions/endogenous data doesn't link them to a GPCR target (e.g.
# covalently-bound chromophores like retinal, which GtP does not tabulate
# as a standard ligand-target interaction). Add new rules here as needed.
# Each rule supports "name_contains" (case-insensitive substrings matched
# against the ligand name) and/or "ligand_ids" (explicit GtP ligand_id list).
CUSTOM_GTP_LIGAND_INCLUSIONS = [
    {
        "reason": "retinal isomers (rhodopsin/opsin chromophore; not recorded by GtP as a ligand-target interaction)",
        "name_contains": ["retinal"],
    },
]


def get_custom_gtp_ligand_ids(gtp_complete_ligands):
    """Return GtP ligand_ids to force-include per CUSTOM_GTP_LIGAND_INCLUSIONS,
    regardless of GPCR-target linkage in interactions/endogenous data."""
    ids = set()
    for rule in CUSTOM_GTP_LIGAND_INCLUSIONS:
        mask = pd.Series(False, index=gtp_complete_ligands.index)
        for substr in rule.get("name_contains", []):
            mask |= gtp_complete_ligands["name"].str.contains(substr, case=False, na=False)
        if rule.get("ligand_ids"):
            mask |= gtp_complete_ligands["ligand_id"].isin(rule["ligand_ids"])
        ids.update(gtp_complete_ligands.loc[mask, "ligand_id"].dropna().unique().tolist())
    return ids


def load_gtp_source_data():
    """Fetch, cache, and normalise all Guide to Pharmacology source CSVs.

    Filters interactions and endogenous data to GPCR-linked ligands only,
    plus any ligand matched by CUSTOM_GTP_LIGAND_INCLUSIONS.
    Returns a dict with keys:
        iuphar_ids              - list[str] GtP target IDs mapping to GPCR proteins
        ligand_ids              - list[str] ligand IDs linked to GPCR targets
                                  (bioactivity + endogenous + custom inclusions, deduplicated)
        gtp_uniprot             - normalised DataFrame
        gtp_complete_ligands    - normalised DataFrame
        gtp_ligand_mapping      - normalised DataFrame
        gtp_interactions        - normalised DataFrame
        gtp_detailed_endogenous - normalised DataFrame
        peptide_ligand_ids      - list[str] ligand IDs from peptides linked to GPCR targets
        gtp_peptides            - normalised DataFrame (all peptides, unfiltered)
    """
    def _fetch(url):
        link = get_or_create_url_cache(url, _GTP_CACHE_TTL)
        df = pd.read_csv(link, dtype=str, header=1)
        df.columns = [c.strip().lower().replace(" ", "_") for c in df.columns]
        return df

    gtp_uniprot             = _fetch("https://www.guidetopharmacology.org/DATA/GtP_to_UniProt_mapping.csv")
    gtp_complete_ligands    = _fetch("https://www.guidetopharmacology.org/DATA/ligands.csv")
    gtp_ligand_mapping      = _fetch("https://www.guidetopharmacology.org/DATA/ligand_id_mapping.csv")
    gtp_interactions        = _fetch("https://www.guidetopharmacology.org/DATA/interactions.csv")
    gtp_detailed_endogenous = _fetch("https://www.guidetopharmacology.org/DATA/endogenous_ligand_detailed.csv")
    gtp_peptides            = _fetch("https://www.guidetopharmacology.org/DATA/peptides.csv")

    # GPCR target IDs: match GtP UniProt IDs against GPCRdb protein accessions
    gpcr_accessions_base = {
        acc.split("-")[0]
        for acc in Protein.objects.filter(
            family__slug__startswith="0", sequence_type__slug="wt"
        ).values_list("accession", flat=True)
    }
    iuphar_ids = list(
        gtp_uniprot.loc[
            gtp_uniprot["uniprotkb_id"].isin(gpcr_accessions_base),
            "gtopdb_iuphar_id",
        ]
        .dropna()
        .unique()
    )

    def _ligand_ids_for(df):
        matched = set(df["target_id"].unique()).intersection(set(iuphar_ids))
        ids = df.loc[df["target_id"].isin(matched), "ligand_id"].unique().tolist()
        return [x for x in ids if x == x]  # drop NaN

    bioactivity_ids = _ligand_ids_for(gtp_interactions)
    endogenous_ids  = _ligand_ids_for(gtp_detailed_endogenous)
    custom_ids      = get_custom_gtp_ligand_ids(gtp_complete_ligands)
    ligand_ids      = list(set(bioactivity_ids + endogenous_ids) | custom_ids)

    # Peptides linked to GPCR targets: filter by uniprot_id column (pipe-separated)
    def _any_gpcr_uniprot(val):
        if not val or (isinstance(val, float)):
            return False
        return any(u.strip() in gpcr_accessions_base for u in str(val).split("|"))

    gpcr_peptide_mask = gtp_peptides["uniprot_id"].apply(_any_gpcr_uniprot)
    peptide_ligand_ids = (
        gtp_peptides.loc[gpcr_peptide_mask, "ligand_id"].dropna().unique().tolist()
    )

    return {
        "iuphar_ids":               iuphar_ids,
        "ligand_ids":               ligand_ids,
        "peptide_ligand_ids":       peptide_ligand_ids,
        "gtp_uniprot":              gtp_uniprot,
        "gtp_complete_ligands":     gtp_complete_ligands,
        "gtp_ligand_mapping":       gtp_ligand_mapping,
        "gtp_interactions":         gtp_interactions,
        "gtp_detailed_endogenous":  gtp_detailed_endogenous,
        "gtp_peptides":             gtp_peptides,
    }


def _iter_sids(source_ids):
    """Yield (slug, str_value) pairs, expanding list values accumulated during dedup."""
    for slug, val in source_ids.items():
        if isinstance(val, list):
            for v in val:
                yield slug, str(v)
        else:
            yield slug, str(val)


def _inchikeys_compatible(entry_ck: Optional[str], child_inchikey: Optional[str]) -> bool:
    """Return True when it is safe to link a source ID from an entry to an existing child.

    Rules (both sides must be non-empty for rejection to fire):
    - Connectivity block (first dash-segment) must match.
    - If both sides have a stereo block (second segment), they must also match.
    A None on either side means no structural evidence to compare — allow by default.
    """
    if not entry_ck or not child_inchikey:
        return True
    e_parts = entry_ck.split("-")
    c_parts = child_inchikey.split("-")
    if e_parts[0] != c_parts[0]:
        return False
    if len(e_parts) >= 2 and len(c_parts) >= 2 and e_parts[1] != c_parts[1]:
        # Allow when either side is truly achiral (UHFFFAOYSA = no stereo specified).
        # This preserves the achiral-short-circuit design where a racemic entry can
        # reuse any existing child with the same connectivity.
        # Only reject when BOTH sides carry specific stereo information that differs.
        if e_parts[1] != "UHFFFAOYSA" and c_parts[1] != "UHFFFAOYSA":
            return False
    return True


_MW_TOL_FRAC = 0.05   # 5 % relative tolerance on molecular weight
_MW_TOL_ABS  = 10.0   # absolute floor (Da) — avoids over-rejection for small molecules


def _props_compatible(entry_props: dict, candidate: "Ligand") -> bool:
    """Return False only when both sides have a value AND it differs beyond tolerance.

    Checks mw (within 5 % or 10 Da), hacc, and hdon (integer exact match).
    logp is deliberately excluded — algorithm differences (clogp vs alogp) make
    it unreliable as a rejection criterion.
    """
    if not entry_props:
        return True
    for field in ("mw", "hacc", "hdon"):
        ev = entry_props.get(field)
        cv = getattr(candidate, field, None)
        if ev is None or cv is None:
            continue
        ev, cv = float(ev), float(cv)
        if field == "mw":
            if abs(ev - cv) > max(_MW_TOL_ABS, _MW_TOL_FRAC * max(ev, cv, 1.0)):
                return False
        else:
            if round(ev) != round(cv):
                return False
    return True


def _prepare_ligand_entries(entries):
    """Steps A/A.5 of ligand-batch processing: RDKit normalization/descriptors and
    UniChem source-ID enrichment. Pure computation + network I/O, no shared DB
    state — safe (and important for throughput) to run without holding the
    ligand DB lock, unlike the resolve/create steps that follow.
    """
    # ── Step A: pre-normalize all entries in memory ───────────────────────────
    for entry in entries:
        entry["_pk"] = None   # parent connectivity key (first InChIKey block)
        entry["_ck"] = None   # child identity key (full normalized InChIKey)
        entry["_parent_norm"] = None

        raw_smiles = entry.get("raw_smiles")
        raw_inchikey = entry.get("raw_inchikey")

        if raw_smiles and not entry.get("skip_normalization"):
            pnorm = normalize_for_parent(raw_smiles)
            cnorm = normalize_for_child(raw_smiles)
            if pnorm:
                entry["_pk"] = pnorm["clean_inchikey"]
                entry["_parent_norm"] = pnorm
            if cnorm:
                if cnorm.get("inchikey"):
                    entry["_ck"] = cnorm["inchikey"]
                entry["_child_mol"] = cnorm.get("mol")
        # Always fall back to raw_inchikey when normalization didn't produce keys
        if not entry["_ck"] and raw_inchikey:
            entry["_ck"] = raw_inchikey
        if not entry["_pk"] and raw_inchikey:
            entry["_pk"] = raw_inchikey.split("-")[0]

        # Pre-compute physicochemical properties for child-matching discrimination.
        # Raw SMILES is preferred so values are consistent with the structure shown to users;
        # parent-normalised SMILES is the fallback when raw SMILES is unavailable or unparseable.
        entry["_entry_props"] = {}
        _smiles_for_props = entry.get("raw_smiles") or (entry.get("_parent_norm") or {}).get("smiles")
        if _smiles_for_props:
            try:
                _pmol = dm.to_mol(_smiles_for_props, sanitize=True)
                if _pmol:
                    entry["_entry_props"] = {
                        "mw":   dm.descriptors.mw(_pmol),
                        "hacc": dm.descriptors.n_hba(_pmol),
                        "hdon": dm.descriptors.n_hbd(_pmol),
                    }
            except Exception:
                pass

        # Pre-compute isotope label so Steps C/D/E/G can distinguish different
        # radioactive variants of the same compound (e.g. [³H] vs [¹²⁵I]).
        entry["_radioactive_label"] = None
        if entry.get("radioactive"):
            _rl = get_radioactive_info(entry.get("name", ""), entry.get("synonyms", "") or "")
            entry["_radioactive_label"] = _rl[0] if _rl else "Yes"

    # ── Step A.5: UniChem source-ID enrichment for structure-less entries ─────
    # Entries with no InChIKey, SMILES, or HELM cannot be matched structurally.
    # Query UniChem to find crossover IDs that may already exist on children in
    # the DB, turning a structure-less dead-end into a successful source-ID match.
    # raw_sequence is intentionally NOT a skip condition: Step C builds no
    # existing_by_sequence child lookup, so sequence-only entries (e.g. peptides
    # with ambiguous_alias=False) can only be matched via source-ID.
    for entry in entries:
        if entry.get("_ck") or entry.get("raw_smiles") or entry.get("raw_helm"):
            continue  # InChIKey/SMILES/HELM matching paths cover this entry
        sids = entry.setdefault("source_ids", {})
        for src_slug, src_val in list(_iter_sids(sids)):
            if src_slug not in unichem_src_types_inv:
                continue
            for hit in match_id_via_unichem(src_slug, str(src_val)):
                if hit["type"] in external_sources and hit["type"] not in sids:
                    sids[hit["type"]] = hit["id"]


def get_or_create_ligands_batch(entries, source=None, lig_type="small-molecule", extended_matching=True):
    """Process a list of ligand entries in batch, minimizing DB round-trips.

    entries: list of dicts with keys:
        name, raw_smiles, raw_inchikey, raw_helm, raw_sequence,
        source_ids (dict), radioactive (bool), synonyms (str)

    Returns list of (parent, child) tuples, one per entry (either may be None on failure).
    """
    _prepare_ligand_entries(entries)

    with _ligand_lock_cm():
        return _resolve_ligand_entries(entries, source=source, lig_type=lig_type)


def _resolve_ligand_entries(entries, source=None, lig_type="small-molecule", extended_matching=True):
    """Steps B-I: DB read/resolve/create for a batch of already-prepared entries.

    Runs under the shared ligand lock during a parallel build (see set_ligand_lock)
    since the check-then-insert logic below relies on that lock, not a DB
    constraint, to avoid two workers creating duplicate Ligand rows.
    """
    if lig_type not in _LIGAND_TYPE_CACHE:
        _LIGAND_TYPE_CACHE[lig_type] = LigandType.objects.get_or_create(
            slug=slugify(lig_type), defaults={"name": lig_type}
        )[0]
    lig_type_obj = _LIGAND_TYPE_CACHE[lig_type]

    # ── Step B: batch parent lookup ───────────────────────────────────────────
    all_clean_keys = {e["_pk"] for e in entries if e["_pk"]}
    all_parent_smiles = {e["_parent_norm"]["smiles"] for e in entries
                         if e.get("_parent_norm") and e["_parent_norm"].get("smiles")}
    all_helms = {e.get("raw_helm") for e in entries if e.get("raw_helm")}
    all_sequences = {e.get("raw_sequence") for e in entries
                     if e.get("raw_sequence") and len(e.get("raw_sequence", "")) > 3}
    all_uniprot_ids = {
        e.get("raw_uniprot") for e in entries
        if e.get("raw_uniprot") and not e.get("raw_sequence") and not e.get("raw_helm")
    }
    # Truly structure-less entries (no InChIKey, helm, sequence, or uniprot) are matched
    # back to existing parents by name, but only when those parents are flagged ambiguous.
    all_ambiguous_names = {
        e["name"] for e in entries
        if not e.get("_pk") and not e.get("raw_helm")
        and not e.get("raw_sequence") and not e.get("raw_uniprot")
    }

    parent_qs = Ligand.objects.filter(parent__isnull=True)
    q = Q()
    if all_clean_keys:
        q |= Q(clean_inchikey__in=all_clean_keys)
    if all_parent_smiles:
        q |= Q(smiles__in=all_parent_smiles)
    if all_helms:
        q |= Q(helm__in=all_helms)
    if all_sequences:
        q |= Q(sequence__in=all_sequences)
    if all_uniprot_ids:
        q |= Q(uniprot__in=all_uniprot_ids)
    if all_ambiguous_names:
        q |= Q(name__in=all_ambiguous_names, ambiguous_alias=True)

    existing_parents = {}
    if q:
        for p in parent_qs.filter(q):
            if p.clean_inchikey:
                existing_parents[p.clean_inchikey] = p
            if p.smiles:
                existing_parents.setdefault(p.smiles, p)
            if p.helm:
                existing_parents.setdefault(p.helm, p)
            if p.sequence:
                existing_parents.setdefault(p.sequence, p)
            if p.uniprot and not p.sequence and not p.helm:
                existing_parents.setdefault(p.uniprot, p)
            if p.ambiguous_alias:
                existing_parents.setdefault(p.name, p)

    # ── Step C: batch child lookup ────────────────────────────────────────────
    all_child_norm_keys = {e["_ck"] for e in entries if e["_ck"]}
    all_raw_iks = {e.get("raw_inchikey") for e in entries if e.get("raw_inchikey")}

    existing_by_norm_ik = {}
    if all_child_norm_keys:
        for c in Ligand.objects.filter(clean_inchikey__in=all_child_norm_keys, parent__isnull=False):
            existing_by_norm_ik[(c.clean_inchikey, c.radioactive)] = c

    existing_by_raw_ik = {}
    if all_raw_iks:
        for c in Ligand.objects.filter(inchikey__in=all_raw_iks, parent__isnull=False):
            existing_by_raw_ik[(c.inchikey, c.radioactive)] = c

    # HELM-based child index — for peptide entries that have a HELM string but no InChIKey.
    # Keyed by helm string; lowest gpcrdb_id wins as the canonical record.
    existing_by_helm = {}  # helm → child Ligand
    all_helms_no_ck = {
        e.get("raw_helm") for e in entries
        if e.get("raw_helm") and not e.get("_ck")
    }
    if all_helms_no_ck:
        for c in (Ligand.objects
                  .filter(helm__in=all_helms_no_ck, parent__isnull=False)
                  .select_related("parent")
                  .order_by("gpcrdb_id")):
            existing_by_helm.setdefault(c.helm, c)

    # SMILES-based child index — fallback for children already in DB with clean_inchikey=None.
    # Only built for entries that still have no _ck after Step A (normalization failed AND
    # no raw_inchikey was available). Skipping entries that have an InChIKey avoids a
    # potentially slow query for every protein/peptide SMILES (can be thousands of chars).
    # Keyed by (smiles, parent_id); lowest gpcrdb_id wins as the canonical record.
    existing_by_smiles = {}  # (smiles, parent_id) → child Ligand
    all_raw_smiles_no_ck = {
        e.get("raw_smiles") for e in entries
        if e.get("raw_smiles") and not e.get("_ck")
    }
    if all_raw_smiles_no_ck:
        for c in (Ligand.objects
                  .filter(smiles__in=all_raw_smiles_no_ck, parent__isnull=False)
                  .select_related("parent")
                  .order_by("gpcrdb_id")):
            existing_by_smiles.setdefault((c.smiles, c.parent_id), c)

    # ── Step D: batch source ID lookup ───────────────────────────────────────
    existing_by_source_id = {}  # (slug, index, is_radioactive) → child Ligand
    for src_slug in external_sources:
        ids_for_src = {
            str(e["source_ids"][src_slug])
            for e in entries
            if src_slug in e.get("source_ids", {})
        }
        if ids_for_src:
            for lid in (LigandID.objects
                        .filter(index__in=ids_for_src, web_resource__slug=src_slug)
                        .select_related("ligand", "ligand__parent")):
                existing_by_source_id[(src_slug, lid.index, lid.ligand.radioactive)] = lid.ligand

    # ── Step D.4: direct name match for the "Apo (no ligand)" sentinel ───────
    # These entries carry no identifying data at all (no pdbe/sequence/source
    # ids), so every occurrence across every structure would otherwise create
    # its own duplicate child under the shared ambiguous parent. Exact (not
    # iexact) match, so it can never be confused with a real ligand.
    existing_apo_child = None
    if any(e.get("name") == APO_NO_LIGAND_NAME for e in entries):
        existing_apo_child = (Ligand.objects
                               .filter(name=APO_NO_LIGAND_NAME, parent__isnull=False)
                               .order_by("gpcrdb_id")
                               .first())

    # ── Step D.5: batch pdbe lookup ──────────────────────────────────────────
    # Checked ahead of every other identifier below (Step E) - a PDB
    # chemical-component code is a strong, unambiguous match on its own.
    existing_by_pdbe = {}  # pdbe value → child Ligand
    all_pdbe_ids = {e["raw_pdbe"] for e in entries if e.get("raw_pdbe")}
    if all_pdbe_ids:
        for c in Ligand.objects.filter(pdbe__in=all_pdbe_ids, parent__isnull=False):
            existing_by_pdbe.setdefault(c.pdbe, c)

    # ── Step D.6: batch sequence lookup (opt-in via seq_and_name_lookup) ─────
    # Only for entries explicitly opting in (peptide/protein ligands with no pdbe
    # id - see get_or_create_ligand's seq_and_name_lookup flag). "X" marks a
    # non-standard/unresolved residue, so a sequence containing it is not a safe
    # unique key - two different molecules could share the same X-containing string.
    existing_by_seq_lookup = {}  # sequence → child Ligand
    all_seq_lookup_seqs = {
        e["raw_sequence"] for e in entries
        if e.get("seq_and_name_lookup") and e.get("raw_sequence") and "X" not in e["raw_sequence"]
    }
    if all_seq_lookup_seqs:
        for c in Ligand.objects.filter(sequence__in=all_seq_lookup_seqs, parent__isnull=False):
            existing_by_seq_lookup.setdefault(c.sequence, c)

    # ── Step E: classify each entry ──────────────────────────────────────────
    entries_needing_new_child = []
    entries_needing_new_parent = []  # implies new child too
    new_ligand_ids_to_create = []

    resolved = []  # (parent, child) per entry — filled below

    def _resolve_parent_for_entry(entry):
        pk = entry.get("_pk")
        pnorm = entry.get("_parent_norm")
        if pk and pk in existing_parents:
            return existing_parents[pk]
        if pnorm and pnorm.get("smiles") and pnorm["smiles"] in existing_parents:
            return existing_parents[pnorm["smiles"]]
        raw_helm = entry.get("raw_helm")
        if raw_helm and raw_helm in existing_parents:
            return existing_parents[raw_helm]
        raw_seq = entry.get("raw_sequence")
        if raw_seq and raw_seq in existing_parents:
            return existing_parents[raw_seq]
        raw_uniprot = entry.get("raw_uniprot")
        if (raw_uniprot and not entry.get("raw_sequence") and not entry.get("raw_helm")
                and raw_uniprot in existing_parents):
            return existing_parents[raw_uniprot]
        # Last resort: name-based match for truly structure-less entries, restricted to
        # ambiguous_alias=True parents so we don't merge distinct ligands by name.
        if (not entry.get("_pk") and not entry.get("raw_helm")
                and not entry.get("raw_sequence") and not entry.get("raw_uniprot")
                and entry["name"] in existing_parents):
            return existing_parents[entry["name"]]
        return None

    for entry in entries:
        entry_radioactive = bool(entry.get("radioactive"))
        entry_radioactive_label = entry.get("_radioactive_label")
        # When the isotope type is unknown ("Yes"), structural matching cannot
        # distinguish different radioactive variants — rely on source-ID only.
        _skip_structural = (entry_radioactive_label == "Yes")

        found_child = None
        if entry.get("name") == APO_NO_LIGAND_NAME:
            # No identifying data exists for this sentinel (no pdbe/sequence/
            # source ids), so every other mechanism below would be a no-op
            # anyway; match/reuse the single existing child directly by exact
            # name instead of creating a duplicate for every apo structure.
            found_child = existing_apo_child
        else:
            # PDBe match? — checked ahead of every other identifier below
            raw_pdbe = entry.get("raw_pdbe")
            if raw_pdbe:
                found_child = existing_by_pdbe.get(raw_pdbe)

        # Sequence / iexact-name match? (opt-in, for entries with no pdbe id)
        if found_child is None and entry.get("seq_and_name_lookup"):
            raw_seq = entry.get("raw_sequence")
            if raw_seq and "X" not in raw_seq:
                found_child = existing_by_seq_lookup.get(raw_seq)

            if found_child is None and entry.get("name"):
                _name_matches = list(
                    Ligand.objects.filter(name__iexact=entry["name"].strip(), parent__isnull=False)
                )
                name_hit = None
                if len(_name_matches) == 1:
                    name_hit = _name_matches[0]
                elif len(_name_matches) > 1:
                    _exact = [c for c in _name_matches if c.name == entry["name"].strip()]
                    name_hit = _exact[0] if _exact else _name_matches[0]

                if name_hit is not None:
                    if not name_hit.sequence:
                        # Nothing on file contradicts the name match — accept it,
                        # and back-fill the sequence for future lookups.
                        if raw_seq:
                            name_hit.sequence = raw_seq
                            name_hit.save(update_fields=["sequence"])
                        found_child = name_hit
                    elif raw_seq and raw_seq == name_hit.sequence and "X" not in raw_seq:
                        found_child = name_hit
                    else:
                        # Same name, different molecule (sequences differ, or either
                        # side is ambiguous) — not a match, but keep it in the family:
                        # a new child will be created under this ligand's parent
                        # instead of via generic parent resolution below.
                        entry["_resolved_parent"] = name_hit.parent

        # Source ID match?
        _source_id_that_found: Optional[tuple] = None
        if found_child is None:
            for src_slug, src_val in _iter_sids(entry.get("source_ids", {})):
                key = (src_slug, src_val, entry_radioactive_label)
                if key in existing_by_source_id:
                    found_child = existing_by_source_id[key]
                    _source_id_that_found = (src_slug, src_val)
                    break

        # If source-ID match found a structurally incompatible child, fall through
        # to InChIKey matching / new child creation so the entry's IDs are not silently lost.
        # Use clean_inchikey (normalized) in preference to inchikey (raw stored) — they can
        # differ for the same compound when different sources use different stereo conventions.
        if (found_child is not None
                and _source_id_that_found is not None
                and (found_child.clean_inchikey or found_child.inchikey)
                and not _inchikeys_compatible(
                    entry.get("_ck") or entry.get("raw_inchikey"),
                    found_child.clean_inchikey or found_child.inchikey)):
            log.info(
                f"Source-ID {_source_id_that_found} found incompatible child "
                f"(entry {entry.get('_ck')!r} vs child {found_child.clean_inchikey or found_child.inchikey!r}); falling through"
            )
            found_child = None

        entry_props = entry.get("_entry_props", {})

        # Child-normalized InChIKey match?
        if not _skip_structural and found_child is None and entry.get("_ck"):
            candidate = existing_by_norm_ik.get((entry["_ck"], entry_radioactive_label))
            if candidate is not None and _props_compatible(entry_props, candidate):
                # Guard: if entry has a specific raw_inchikey (e.g. an E/Z stereo block that
                # tautomer canonicalization stripped from _ck), verify it is compatible with
                # the candidate's inchikey before accepting the match.
                _entry_raw_ik = entry.get("raw_inchikey") or ""
                _erp = _entry_raw_ik.split("-")
                if (len(_erp) >= 2
                        and _erp[1] not in ("", "UHFFFAOYSA")
                        and not _inchikeys_compatible(_entry_raw_ik, candidate.inchikey or "")):
                    log.info(
                        f"norm-ik match rejected: entry raw_inchikey {_entry_raw_ik!r} "
                        f"incompatible with candidate inchikey {candidate.inchikey!r}"
                    )
                else:
                    found_child = candidate
            elif candidate is not None and not _props_compatible(entry_props, candidate):
                log.warning(
                    f"InChIKey match rejected by property mismatch: entry {entry.get('name')!r} "
                    f"vs child gpcrdb_id={candidate.gpcrdb_id} — "
                    f"entry mw={entry_props.get('mw')}, child mw={candidate.mw}"
                )

        # Raw InChIKey match?
        if not _skip_structural and found_child is None and entry.get("raw_inchikey"):
            candidate = existing_by_raw_ik.get((entry["raw_inchikey"], entry_radioactive_label))
            if candidate is not None and _props_compatible(entry_props, candidate):
                found_child = candidate

        # SMILES direct match — last resort when no InChIKey lookup succeeded but an
        # existing child with the same SMILES and parent exists in the DB.
        if not _skip_structural and found_child is None and entry.get("raw_smiles") and entry.get("_pk"):
            for (smi, par_id), candidate in existing_by_smiles.items():
                if (smi == entry["raw_smiles"]
                        and candidate.parent
                        and candidate.parent.clean_inchikey == entry["_pk"]
                        and candidate.radioactive == entry_radioactive_label
                        and _props_compatible(entry_props, candidate)):
                    found_child = candidate
                    break

        # HELM match — for peptides with HELM but no InChIKey
        if not _skip_structural and found_child is None and entry.get("raw_helm") and not entry.get("_ck"):
            candidate = existing_by_helm.get(entry["raw_helm"])
            if candidate is not None and candidate.radioactive == entry_radioactive_label:
                found_child = candidate

        # Achiral short-circuit
        if not _skip_structural and found_child is None and entry.get("_ck") and entry.get("_pk"):
            mol = entry.get("_child_mol")
            if mol is not None:
                from rdkit import Chem as _Chem
                _Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
                _no_chiral_centers = not _Chem.FindMolChiralCenters(mol, includeUnassigned=True)
                _no_ez_bonds = not any(
                    b.GetStereo() not in (_Chem.BondStereo.STEREONONE, _Chem.BondStereo.STEREOANY)
                    for b in mol.GetBonds()
                )
                _raw_ik_ac = entry.get("raw_inchikey") or entry.get("_ck") or ""
                _raw_ik_ac_parts = _raw_ik_ac.split("-")
                _has_specific_stereo_ik = (
                    len(_raw_ik_ac_parts) >= 2
                    and _raw_ik_ac_parts[1] not in ("", "UHFFFAOYSA")
                )
                if _no_chiral_centers and _no_ez_bonds and not _has_specific_stereo_ik:
                    if entry_radioactive_label:
                        _rac_filter = {"radioactive": entry_radioactive_label}
                    else:
                        _rac_filter = {"radioactive__isnull": True}
                    candidates = list(Ligand.objects.filter(
                        clean_inchikey__startswith=entry["_pk"],
                        parent__isnull=False,
                        **_rac_filter,
                    )[:1])
                    if candidates and _props_compatible(entry_props, candidates[0]):
                        found_child = candidates[0]

        if found_child is not None:
            # Reassignment check
            correct_pk = entry.get("_pk")
            if (correct_pk and found_child.parent is not None and
                    found_child.parent.clean_inchikey and
                    found_child.parent.clean_inchikey != correct_pk):
                correct_parent = existing_parents.get(correct_pk)
                if correct_parent:
                    found_child.parent = correct_parent
                    found_child.save(update_fields=["parent"])
                    log.warning(
                        f"Reassigned child {found_child.name!r} (gpcrdb_id={found_child.gpcrdb_id}) "
                        f"to parent {correct_pk!r}"
                    )
                else:
                    log.error(
                        f"Child {found_child.name!r} expected parent {correct_pk!r} but not found"
                    )

            # Back-fill missing name from entry, or fall back to parent name
            if not found_child.name:
                new_name = (entry.get("name") or "").strip()
                if not new_name and found_child.parent:
                    new_name = found_child.parent.name or ""
                if new_name:
                    found_child.name = new_name
                    found_child.save(update_fields=["name"])

            # Back-fill uniprot on existing children that lack it
            if not found_child.uniprot and entry.get("raw_uniprot"):
                found_child.uniprot = entry["raw_uniprot"]
                found_child.save(update_fields=["uniprot"])

            # Back-fill pdbe on existing children that lack it (however they were matched)
            if not found_child.pdbe and entry.get("raw_pdbe") and len(entry["raw_pdbe"]) <= _PDBE_MAX_LEN:
                found_child.pdbe = entry["raw_pdbe"]
                found_child.save(update_fields=["pdbe"])

            # Queue new LigandIDs — skip the DB fetch entirely when the entry has no source IDs
            entry_source_ids = entry.get("source_ids", {})
            if entry_source_ids:
                existing_child_id_keys = set(
                    found_child.ids.values_list("web_resource__slug", "index")
                )
                for src_slug, src_val in _iter_sids(entry_source_ids):
                    if (src_slug, src_val) not in existing_child_id_keys:
                        new_ligand_ids_to_create.append((found_child, src_slug, src_val))

            resolved.append((found_child.parent, found_child))
        else:
            # A name match with an incompatible sequence (above) already pins this
            # entry's parent so the new child stays in that ligand's family.
            existing_parent = entry.get("_resolved_parent") or _resolve_parent_for_entry(entry)
            if existing_parent is None:
                entries_needing_new_parent.append(entry)
            else:
                entry["_resolved_parent"] = existing_parent
                entries_needing_new_child.append(entry)

    # ── Step F: bulk-create new parents ──────────────────────────────────────
    unique_new_parents = {}  # _pk → entry (first occurrence wins)
    for entry in entries_needing_new_parent:
        pk = entry.get("_pk") or entry.get("name")
        if pk not in unique_new_parents:
            unique_new_parents[pk] = entry

    new_parent_objects = []
    pk_to_entry = {}
    for pk, entry in unique_new_parents.items():
        pnorm = entry.get("_parent_norm") or {}
        p = Ligand()
        p.name = entry["name"]
        p.ambiguous_alias = False
        p.ligand_type = lig_type_obj
        p.smiles = pnorm.get("smiles")
        p.inchikey = pnorm.get("inchikey")
        p.clean_inchikey = pnorm.get("clean_inchikey")
        if not p.clean_inchikey and entry.get("_pk"):
            p.clean_inchikey = entry["_pk"]
        raw_seq = entry.get("raw_sequence")
        if raw_seq:
            p.sequence = raw_seq
        raw_helm = entry.get("raw_helm")
        if raw_helm:
            p.helm = raw_helm
        p.uniprot = entry.get("raw_uniprot")
        p.ambiguous_alias = not bool(
            p.clean_inchikey or p.smiles or p.helm or p.sequence or p.uniprot
        )
        p.parent = None
        new_parent_objects.append(p)
        pk_to_entry[pk] = entry

    if new_parent_objects:
        next_id = allocate_next_gpcrdb_id()
        for i, p in enumerate(new_parent_objects):
            p.gpcrdb_id = next_id + i

        # Build pk → gpcrdb_id before bulk_create (ids already assigned above).
        # unique_new_parents and new_parent_objects share the same insertion order
        # (Python 3.7+ dicts preserve it), so zip is safe.
        pk_to_gpcrdb_id = dict(
            zip(unique_new_parents.keys(), [p.gpcrdb_id for p in new_parent_objects])
        )

        Ligand.objects.bulk_create(new_parent_objects)

        # Re-fetch ALL new parents by gpcrdb_id — works even when clean_inchikey is NULL.
        # Index by every available structural key so _resolve_parent_for_entry finds them.
        _fresh_gids = list(pk_to_gpcrdb_id.values())
        _gpcrdb_to_parent = {}
        for p in Ligand.objects.filter(gpcrdb_id__in=_fresh_gids, parent__isnull=True):
            _gpcrdb_to_parent[p.gpcrdb_id] = p
            if p.clean_inchikey:
                existing_parents[p.clean_inchikey] = p
            if p.helm:
                existing_parents.setdefault(p.helm, p)
            if p.sequence:
                existing_parents.setdefault(p.sequence, p)
            if p.uniprot and not p.sequence and not p.helm:
                existing_parents.setdefault(p.uniprot, p)
            if p.ambiguous_alias:
                existing_parents.setdefault(p.name, p)

        # Direct pk → parent mapping for entries with no structural identifiers at all.
        pk_to_new_parent = {
            pk: _gpcrdb_to_parent[gid]
            for pk, gid in pk_to_gpcrdb_id.items()
            if gid in _gpcrdb_to_parent
        }
    else:
        pk_to_new_parent = {}

    # Resolve parents for entries that needed a new parent
    for entry in entries_needing_new_parent:
        # Prefer structural-key lookup (helm / sequence / uniprot populated above)
        resolved_parent = _resolve_parent_for_entry(entry)
        # Fall back to the gpcrdb_id-backed direct mapping for structure-less entries
        if resolved_parent is None:
            pk = entry.get("_pk") or entry.get("name")
            resolved_parent = pk_to_new_parent.get(pk)
        entry["_resolved_parent"] = resolved_parent
        entries_needing_new_child.append(entry)

    # ── Step G: bulk-create new children ─────────────────────────────────────
    # Deduplicate within this batch: key is (identity, is_radioactive) so that
    # same-molecule entries collapse to one, while radioactive/non-radioactive
    # variants of the same molecule remain distinct.
    seen_child_cks: dict = {}  # ck_key → surviving entry
    unique_entries_needing_new_child = []
    for entry in entries_needing_new_child:
        _rl = entry.get("_radioactive_label")
        ck = (
            entry.get("_ck") or entry.get("raw_smiles") or entry.get("name"),
            (_rl, entry.get("name")) if _rl == "Yes" else _rl,
        )
        if ck not in seen_child_cks:
            seen_child_cks[ck] = entry
            unique_entries_needing_new_child.append(entry)
        else:
            # Merge source IDs from dropped duplicate into the surviving entry.
            # If the slug already exists with a different value, accumulate as a list
            # so that both IDs (e.g. two pubchem CIDs for the same molecule) survive.
            surviving = seen_child_cks[ck]
            _sids = surviving.setdefault("source_ids", {})
            for slug, val in entry.get("source_ids", {}).items():
                if slug not in _sids:
                    _sids[slug] = val
                elif _sids[slug] != val:
                    existing = _sids[slug]
                    if isinstance(existing, list):
                        if val not in existing:
                            existing.append(val)
                    else:
                        _sids[slug] = [existing, val]
    entries_needing_new_child = unique_entries_needing_new_child

    # Pre-flight: catch any _ck that already exists in the DB regardless of
    # parent status (e.g. legacy records stored with a full InChIKey as
    # clean_inchikey and parent_id=NULL that Step C's parent__isnull=False
    # filter silently skips).
    if entries_needing_new_child:
        _preflight_cks = {e["_ck"] for e in entries_needing_new_child if e.get("_ck")}
        if _preflight_cks:
            _already_in_db = {
                c.clean_inchikey: c
                for c in Ligand.objects.filter(clean_inchikey__in=_preflight_cks)
            }
            if _already_in_db:
                _surviving = []
                for entry in entries_needing_new_child:
                    ck = entry.get("_ck")
                    if ck and ck in _already_in_db:
                        existing = _already_in_db[ck]
                        # Mirror the radioactive guard from Step E: a radioactive entry
                        # must not be merged into a non-radioactive child (or vice-versa)
                        # even when they share an InChIKey (isotope not in SMILES).
                        if bool(entry.get("radioactive")) != bool(existing.radioactive):
                            if existing.radioactive and not entry.get("radioactive"):
                                # Radioactive variant was created first but should NOT own
                                # clean_inchikey — transfer ownership to the non-radioactive
                                # entry by clearing the radioactive row's clean_inchikey now,
                                # before the non-radioactive child is bulk-created.
                                existing.clean_inchikey = None
                                existing.save(update_fields=["clean_inchikey"])
                                # Don't clear _ck — non-radioactive child will own it.
                            else:
                                # Non-radioactive exists first; radioactive must not reuse
                                # its clean_inchikey.
                                entry["_ck"] = None
                            _surviving.append(entry)
                            continue
                        elif (entry.get("_radioactive_label") is not None
                              and entry.get("_radioactive_label") != existing.radioactive):
                            # Different isotope label — create a separate child rather than
                            # merging. The new child does not own clean_inchikey.
                            entry["_ck"] = None
                            _surviving.append(entry)
                            continue
                        resolved.append((existing.parent or existing, existing))
                        for src_slug, src_val in _iter_sids(entry.get("source_ids", {})):
                            new_ligand_ids_to_create.append((existing, src_slug, src_val))
                    else:
                        _surviving.append(entry)
                entries_needing_new_child = _surviving

    new_child_objects = []
    new_child_entries = []
    for entry in entries_needing_new_child:
        resolved_parent = entry.get("_resolved_parent")
        pnorm = entry.get("_parent_norm")
        mol = entry.get("_child_mol")

        # Compute RDKit properties from raw SMILES for consistency with the structure
        # displayed to users; fall back to parent-normalised SMILES if raw SMILES fails.
        rdkit_props = {}
        _smiles_for_props = entry.get("raw_smiles") or (pnorm or {}).get("smiles")
        if _smiles_for_props:
            try:
                pmol = dm.to_mol(_smiles_for_props, sanitize=True)
                if pmol:
                    rdkit_props = {
                        "mw":              dm.descriptors.mw(pmol),
                        "rotatable_bonds": dm.descriptors.n_rotatable_bonds(pmol),
                        "hacc":            dm.descriptors.n_hba(pmol),
                        "hdon":            dm.descriptors.n_hbd(pmol),
                        "logp":            dm.descriptors.clogp(pmol),
                    }
            except Exception:
                pass

        c = Ligand()
        _name_raw = entry.get("name")
        if isinstance(_name_raw, float):
            _name_raw = str(int(_name_raw)) if _name_raw.is_integer() else str(_name_raw)
        _entry_name = (_name_raw or "").strip()
        if not _entry_name and resolved_parent:
            _entry_name = resolved_parent.name or ""
        c.name = _entry_name
        c.ambiguous_alias = False
        c.ligand_type = lig_type_obj
        c.smiles = entry.get("raw_smiles")
        c.inchikey = entry.get("raw_inchikey") or entry.get("_ck")
        c.clean_inchikey = entry.get("_ck")
        c.sequence = entry.get("raw_sequence")
        c.helm = entry.get("raw_helm")
        c.uniprot = entry.get("raw_uniprot")
        _raw_pdbe = entry.get("raw_pdbe")
        if _raw_pdbe and len(_raw_pdbe) <= _PDBE_MAX_LEN:
            c.pdbe = _raw_pdbe
        c.ambiguous_alias = not bool(
            c.clean_inchikey or c.smiles or c.helm or c.sequence or c.uniprot
        )
        c.parent = resolved_parent
        c.source = source
        c.stereo_status = infer_stereo_status(mol)
        c.radioactive = entry.get("_radioactive_label")
        for field, val in rdkit_props.items():
            setattr(c, field, val)

        new_child_objects.append(c)
        new_child_entries.append(entry)

    if new_child_objects:
        next_id = allocate_next_gpcrdb_id()
        for i, c in enumerate(new_child_objects):
            c.gpcrdb_id = next_id + i
        Ligand.objects.bulk_create(new_child_objects)

        # Refresh to get DB-assigned PKs — clean_inchikey-indexed for most children
        fresh_cks = [c.clean_inchikey for c in new_child_objects if c.clean_inchikey]
        new_children_by_ck = {
            c.clean_inchikey: c
            for c in Ligand.objects.filter(clean_inchikey__in=fresh_cks, parent__isnull=False)
        } if fresh_cks else {}

        # Also re-fetch children with clean_inchikey=None (radioactive variants whose _ck
        # was cleared by the pre-flight check) by gpcrdb_id so their pk is guaranteed set.
        none_ck_gids = [c.gpcrdb_id for c in new_child_objects if not c.clean_inchikey and c.gpcrdb_id]
        new_children_by_gpcrdb = {}
        if none_ck_gids:
            for c in Ligand.objects.filter(gpcrdb_id__in=none_ck_gids, parent__isnull=False):
                new_children_by_gpcrdb[c.gpcrdb_id] = c

        for entry, child_obj in zip(new_child_entries, new_child_objects):
            found_child = (
                new_children_by_ck.get(child_obj.clean_inchikey)
                or new_children_by_gpcrdb.get(child_obj.gpcrdb_id)
                or child_obj
            )
            resolved.append((found_child.parent, found_child))

            for src_slug, src_val in _iter_sids(entry.get("source_ids", {})):
                new_ligand_ids_to_create.append((found_child, src_slug, src_val))

    # ── Step H: bulk-create new LigandIDs ────────────────────────────────────
    if not _WEB_RESOURCE_CACHE:
        for wr in WebResource.objects.filter(slug__in=list(external_sources) + ["uniprot"]):
            _WEB_RESOURCE_CACHE[wr.slug] = wr
    lid_objects = []
    for child, src_slug, src_val in new_ligand_ids_to_create:
        wr = _WEB_RESOURCE_CACHE.get(src_slug)
        if wr and child.pk:
            lid_objects.append(LigandID(ligand=child, index=src_val, web_resource=wr))
    if lid_objects:
        LigandID.objects.bulk_create(lid_objects, ignore_conflicts=True)

    # ── Step I: update parents that still lack structural data ────────────────
    parents_to_update = {p.pk: (p, c) for p, c in resolved if p and c}
    for parent, child in parents_to_update.values():
        if not parent.smiles or not parent.clean_inchikey:
            update_parent(parent, child)

    return resolved

def get_or_create_ligand(name, ids={}, lig_type="small-molecule", unichem=False,
                          extended_matching=True, source=None, helm=None,
                          radioactive=False, synonyms=False, seq_and_name_lookup=False):
    """Thin wrapper around get_or_create_ligands_batch() for single-ligand call sites."""
    # Normalise / clean the ids dict (CAS → PubChem, type coercions, etc.)
    cas_to_cid_url = "http://www.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pccompound&retmax=100&term=$index"
    cas_cache_dir = ["CAS", "cas_codes"]
    ids = dict(ids)  # don't mutate caller's dict

    for type_id in list(ids.keys()):
        if ids[type_id] is None or ids[type_id] == "None" or ids[type_id] == "":
            del ids[type_id]
        elif type_id == "pdbe":
            # Keep as a string (never numeric-coerced) - some PDB chemical-component
            # codes are all-digit (e.g. "140", "836") and must not lose type/leading
            # zeros before being compared against the Ligand.pdbe CharField.
            ids[type_id] = str(ids[type_id]).strip()
        elif isinstance(ids[type_id], str) and ids[type_id].isnumeric():
            ids[type_id] = int(ids[type_id])
        elif isinstance(ids[type_id], str) and is_float(ids[type_id]):
            ids[type_id] = int(float(ids[type_id]))
        elif isinstance(ids[type_id], str):
            ids[type_id] = ids[type_id].strip()
            if type_id == "chembl_ligand":
                ids[type_id] = ids[type_id].upper()
            if type_id == "CAS":
                data = fetch_from_web_api(cas_to_cid_url, ids[type_id], cas_cache_dir, xml=True)
                if data:
                    try:
                        ids["pubchem"] = int(data[3][0].text)
                    except Exception:
                        pass
                del ids["CAS"]
        elif isinstance(ids[type_id], list) and len(ids[type_id]) == 1:
            ids[type_id] = ids[type_id][0]

    _source_ids = {k: v for k, v in ids.items() if k in external_sources}
    if ids.get("uniprot"):
        _source_ids["uniprot"] = ids["uniprot"]

    raw_pdbe = ids.get("pdbe")
    if raw_pdbe is not None and str(raw_pdbe).lower() in ("pep", "apo"):
        raw_pdbe = None  # sentinel names in the annotation data, not real PDBe codes

    entry = {
        "name": name,
        "raw_smiles": ids.get("smiles"),
        "raw_inchikey": ids.get("inchikey"),
        "raw_helm": helm if helm not in (None, "None", "nan") else None,
        "raw_sequence": ids.get("sequence"),
        "raw_uniprot": ids.get("uniprot"),
        "raw_pdbe": raw_pdbe,
        "source_ids": _source_ids,
        "radioactive": radioactive,
        "synonyms": synonyms,
        "seq_and_name_lookup": seq_and_name_lookup,
    }
    results = get_or_create_ligands_batch(
        [entry], source=source, lig_type=lig_type, extended_matching=extended_matching
    )
    if results:
        return results[0][1]  # return child
    return None

unichem_src_types = {"1": "chembl_ligand", "2": "drugbank", "3": "pdb", "4": "gtoplig", "22": "pubchem", "34": "drug_central"}
unichem_src_types_inv = {v: k for k, v in unichem_src_types.items()}

def get_radioactive_info(name, synonyms):
    # Regex for radioactive isotopes
    regexes = {"isotope":"\[(<sup>\d+</sup>[A-Za-z]{1,2})\]"}#, "other":"\[(?=[^\]]*<sup>)([^\]]+)\]"}
    for _, regex in regexes.items():
        isotope_search = re.findall(regex, name)
        if len(isotope_search)>0:
            break

    isotopes = []
    for i in isotope_search:
        isotopes.append(i)
    # Check ligand name synonyms
    found = False
    if len(isotopes)==0 and synonyms:
        for syn in synonyms.split("|"):
            for _, regex in regexes.items():
                isotope_search = re.findall(regex, syn)
                for i in isotope_search:
                    print(i)
                    isotopes.append(i)
                    found = True
                if found:
                    break
            if found:
                break
    return isotopes

def filter_ligand_more_on_name(queryset, name):
    '''Tries to match on name or iexact name from queryset, otherwise returns first element in queryset'''
    result_name_iexact_filtered = queryset.filter(name__iexact=name, parent__isnull = False)
    if result_name_iexact_filtered.count()>1:
        result_name_filtered = queryset.filter(name=name, parent__isnull = False)
        if result_name_filtered.count()>0:
            return result_name_filtered.first()
    elif result_name_iexact_filtered.count()==1:
        return result_name_iexact_filtered.first()
    else:
        return queryset.first()

def match_id_via_unichem(type, id):
    """Query UniChem API v1 for crossover identifiers.

    type: source slug ('chembl_ligand', 'pubchem', 'drugbank', 'gtoplig',
                       'drug_central') or 'inchikey'
    id:   compound identifier in that source database, or an InChIKey string.
    Returns list of {"type": <slug>, "id": <str>} dicts.
    """
    results = []
    cache_path = ['unichem', f'id_match_{type}']
    cached = fetch_from_cache(cache_path, id)
    if cached is not None:
        data = cached
    else:
        if type == "inchikey":
            payload = {"type": "inchikey", "compound": id}
        elif type in unichem_src_types_inv:
            payload = {
                "type": "sourceID",
                "compound": str(id),
                "sourceID": int(unichem_src_types_inv[type]),
            }
        else:
            return results
        try:
            resp = requests.post(
                "https://www.ebi.ac.uk/unichem/api/v1/compounds",
                json=payload,
                timeout=10,
            )
            if resp.status_code != 200:
                return results
            data = resp.json()
            save_to_cache(cache_path, id, data)
        except Exception:
            return results

    compounds = data if isinstance(data, list) else (data or {}).get("compounds", [])
    for compound in compounds:
        for source in compound.get("sources", []):
            src_id_key = str(source.get("srcId", ""))
            src_compound_id = source.get("srcCompoundId")
            if src_id_key in unichem_src_types and src_compound_id:
                results.append({
                    "type": unichem_src_types[src_id_key],
                    "id": str(src_compound_id),
                })
    return results


def get_ligand_by_id(type, id, uniprot = None, forced=True, name=None):
    if forced:
        if uniprot is None:
            result = Ligand.objects.filter(ids__index=str(id), ids__web_resource__slug=type, parent__isnull=False)
        else:
            result = Ligand.objects.filter(ids__index=str(id), ids__web_resource__slug=type, uniprot__contains=uniprot.upper(), parent__isnull=False)
    else:
        if uniprot is None:
            result = Ligand.objects.filter(ids__index=str(id), ids__web_resource__slug=type)
        else:
            result = Ligand.objects.filter(ids__index=str(id), ids__web_resource__slug=type, uniprot__contains=uniprot.upper())

    if result.count() > 1 and name and name != "None":
        # Same source id can be shared by a compound and its radiolabeled variant
        # (e.g. "DAMGO" vs "[3H]DAMGO") - narrow down by name when possible.
        name_matched = result.filter(name__iexact=name)
        if name_matched.exists():
            result = name_matched

    if result.count() > 0:
        # For drugs we allow multiple entries because of stereochemistry if drug is racemic
        if result.count() > 1 and type not in ["drugbank", "drug_central", "pubchem"]:
            print("Multiple entries for the same ID - This should never happen - error", type, id)
        return result.first()
    else:
        return None

def create_ligand_from_id(name, type, id, lig_type, radioactive=False):
    # create new ligand
    ligand = Ligand()
    ligand.name = name
    ligand.ambiguous_alias = False
    ligand.ligand_type = LigandType.objects.get_or_create(slug=slugify(lig_type), defaults={'name': lig_type})[0]

    # Obtain external data (if needed)
    try:

        # TODO - continue here and add peptide/protein support for web services

        if type == "sequence":
            ligand.sequence = sequence
        elif type == "uniprot":
            ligand.uniprot = id
        elif type == "smiles":
            ligand.smiles = id
        elif type == "pubchem":
            # check cache
            cache_dir = ['pubchem', 'cid', 'property']
            url = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/$index/property/IsomericSMILES,InChIKey/json'
            if ligand.name == "":
                url = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/$index/property/Title,IsomericSMILES,InChIKey/json'
            pubchem = fetch_from_web_api(url, id, cache_dir)

            if pubchem != False and pubchem['PropertyTable']['Properties'][0]:
                properties = pubchem['PropertyTable']['Properties'][0]
                if len(properties['IsomericSMILES'])>0:
                    ligand.smiles = properties['IsomericSMILES']
                if len(properties['InChIKey'])>0:
                    ligand.inchikey = properties['InChIKey']
                if ligand.name == "":
                    ligand.name = properties['Title']
        elif type == "gtoplig":
            # check cache
            cache_dir = ['guidetopharmacology', 'ligand_structure']
            url = 'http://www.guidetopharmacology.org/services/ligands/$index/structure'
            gtop = fetch_from_web_api(url, id, cache_dir)

            if gtop:
                if len(gtop["smiles"])>0:
                    ligand.smiles = gtop["smiles"]
                if len(gtop["inchiKey"])>0:
                    ligand.inchikey = gtop["inchiKey"]
                if len(gtop["oneLetterSeq"])>0 and len(gtop["oneLetterSeq"])<1000:
                    ligand.sequence = gtop["oneLetterSeq"]
                    url = 'https://www.guidetopharmacology.org/services/ligands/$index/databaseLinks?database=UniProtKB'
                    gtop_uniprot = fetch_from_web_api(url, id, cache_dir)
                    if gtop_uniprot and len(gtop["accession"])>0:
                        ligand.uniprot = gtop["accession"]
        elif type == "chembl_ligand":
            # check cache
            cache_dir = ['chembl', 'ligand_entry']
            url = 'https://www.ebi.ac.uk/chembl/api/data/molecule/$index/?format=json'
            chembl = fetch_from_web_api(url, id, cache_dir)

            # NOTE: canonical SMILES from ChEMBL are isomeric smiles
            if chembl and "molecule_structures" in chembl:
                if len(chembl["molecule_structures"]["canonical_smiles"]) > 0:
                    ligand.smiles = chembl["molecule_structures"]["canonical_smiles"]
                if len(chembl["molecule_structures"]["standard_inchi_key"]) > 0:
                    ligand.inchikey = chembl["molecule_structures"]["standard_inchi_key"]
            if ligand.name == "" and chembl and "pref_name" in chembl and chembl["pref_name"] is not None:
                ligand.name = chembl["pref_name"]
            elif ligand.name == "":
                ligand.name = id
    except:
        skip=True

    # perform RDkit calculations
    if lig_type == "small-molecule" and ligand.smiles is not None and ligand.smiles != "":
        input_mol = dm.to_mol(ligand.smiles, sanitize=True)
        if input_mol:
            # Check if InChIKey has been set
            if ligand.inchikey is None:
                ligand.inchikey = dm.to_inchikey(input_mol)

            # Check if cleaned InChIKey has been set
            # if ligand.clean_inchikey is None:
                # ligand.clean_inchikey = get_cleaned_inchikey(ligand.smiles)

            # Calculate RDkit properties
            ligand.mw = dm.descriptors.mw(input_mol)
            ligand.rotatable_bonds = dm.descriptors.n_rotatable_bonds(input_mol)
            ligand.hacc = dm.descriptors.n_hba(input_mol)
            ligand.hdon = dm.descriptors.n_hbd(input_mol)
            ligand.logp = dm.descriptors.clogp(input_mol)
        else:
            print("Issues with molecule", name)
            # TODO - try again if we have other IDs

        # Before storing - check one more time if we already have this ligand
        if ligand.inchikey is not None and not radioactive:
            checkligand = get_ligand_by_inchikey(ligand.inchikey, ligand.name)
            if checkligand is not None:
                return checkligand

        # if ligand.clean_inchikey is not None and ligand.clean_inchikey != ligand.inchikey and not radioactive:
        #     checkligand = get_ligand_by_inchikey(ligand.clean_inchikey, ligand.name)
        #     if checkligand is not None:
        #         return checkligand

        return ligand
    elif lig_type != "small-molecule" or type == "":
        if type == "smiles" and ligand.inchikey is None:
            input_mol = dm.to_mol(ligand.smiles, sanitize=True)
            if input_mol:
                ligand.inchikey = dm.to_inchikey(input_mol)
                if ligand.inchikey is not None and not radioactive:
                    checkligand = get_ligand_by_inchikey(ligand.inchikey, ligand.name)
                    if checkligand is not None:
                        return checkligand

        # Before storing - check one more time if we already have this ligand based on sequence
        # TODO implement this double check
        # https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/71728433/XML
        # USE LXML for parsing => https://stackoverflow.com/questions/21628290/parsing-xml-with-xpath-in-python-3


        # https://www.ebi.ac.uk/chembl/api/data/molecule/CHEMBL59?format=xml&_=1643906691331
        # https://www.ebi.ac.uk/chembl/api/data/molecule/CHEMBL266481?format=xml
        #=> PEPTIDE parse HELM notation => no uniprot

        return ligand
    else:
        # Ligand could not be created
        return None

def get_ligand_by_inchikey(inchikey, name=None):
    if not inchikey:
        return None
    result = Ligand.objects.filter((Q(inchikey=inchikey) & Q(parent__isnull=False)))
    if result.count() > 0:
        if result.count() > 1:
            print("Multiple entries for the same InChIkey - This should never happen - error", inchikey)
            if name:
                return filter_ligand_more_on_name(result, name)
        return result.first()
    else:
        return None

def get_cleaned_inchikey(smiles):
    inchikey = None
    try:
        input_mol = dm.to_mol(smiles, sanitize=True)
        if input_mol:
            # Sligthly odd set of steps to get optimal cleaned and dominant tautomer
            input_mol = dm.to_neutral(input_mol)
            input_mol = dm.fix_mol(input_mol)
            cleaned_smiles = dm.standardize_smiles(dm.to_smiles(input_mol), tautomer=True)
            cleaned_mol = dm.to_mol(cleaned_smiles, sanitize=False)
            inchikey = dm.to_inchikey(cleaned_mol)
    except:
        skip = True
    return inchikey

def resolve_pubchem_SID(sid):
    cache_dir = ['pubchem', 'ligand_sid']
    url = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/substance/sid/$index/cids/json'
    pubchem = fetch_from_web_api(url, sid, cache_dir)
    if pubchem and "InformationList" in pubchem and "Information" in pubchem["InformationList"]:
        for entry in pubchem["InformationList"]["Information"]:
            if "CID" in entry:
                return entry["CID"][0]
    return None

def is_float(element):
    try:
        float(element)
        return True
    except ValueError:
        return False

def standardize_smiles(smiles):
    """
    Converts a SMILES string to an RDKit molecule, standardizes it using
    the chembl_structure_pipeline, and returns the standardized molecule.
    """
    # Create an RDKit molecule from the SMILES:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        print(f"Invalid SMILES: {smiles}")
        return None

    # Standardize the molecule.
    # The function returns a tuple: (standardized_mol, standardization_status)
    std_mol = standardizer.standardize_mol(mol)

    # Optionally, you can remove stereochemistry (if that’s your goal):
    Chem.RemoveStereochemistry(std_mol)

    smiles = Chem.MolToSmiles(std_mol, isomericSmiles=False)

    return smiles

# Example usage:
# smiles_example = "CC(=O)OC1=CC=CC=C1C(=O)O"  # Aspirin (for example)
# std_mol = standardize_smiles(smiles_example)
# if status == 'OK':
#     print("Standardized SMILES:", Chem.MolToSmiles(std_mol))
# else:
#     print("Standardization status:", status)
#
# ### HOW TO USE IT
#
# # Two SMILES strings that depict the same molecule but with different stereochemistry:
# smiles1 = "CN1[C@H]2C[C@H](OC(=O)[C@H](CO)c3ccccc3)C[C@@H]1[C@@H]1O[C@H]12"
# smiles2 = "CN1C2CC(CC1C3C2O3)OC(=O)C(CO)C4=CC=CC=C4"
#
# # Standardize the molecules:
# std_mol1 = standardize_smiles(smiles1)
# std_mol2 = standardize_smiles(smiles2)
#
# # Convert standardized molecules to canonical SMILES:
# std_smiles1 = Chem.MolToSmiles(std_mol1, isomericSmiles=False)
# std_smiles2 = Chem.MolToSmiles(std_mol2, isomericSmiles=False)

def allocate_next_gpcrdb_id():
    """
    Returns the next incremental gpcrdb_id (max + 1) under a row-level lock.
    Requires a DB that supports SELECT ... FOR UPDATE (most do).
    """
    with transaction.atomic():
        # Lock all rows touched by the aggregate to serialize concurrent writers
        # (Most DBs will lock the index range under the hood.)
        current_max = (
            Ligand.objects.select_for_update()
            .aggregate(m=Max("gpcrdb_id"))["m"]
        ) or 0
        return current_max + 1

def ensure_gpcrdb_id(ligand):
    if ligand is not None and ligand.gpcrdb_id is None:
        ligand.gpcrdb_id = allocate_next_gpcrdb_id()
        ligand.save(update_fields=["gpcrdb_id"])

def generate_parent(name, ids, lig_type):
    ligand = Ligand()
    ligand.name = name
    ligand.ambiguous_alias = False
    ligand.ligand_type = LigandType.objects.get_or_create(slug=slugify(lig_type), defaults={'name': lig_type})[0]

    if "sequence" in ids:
        ligand.sequence = ids["sequence"]
    if "uniprot" in ids:
        ligand.uniprot = ids["uniprot"]
    if "smiles" in ids:
        ligand.smiles = standardize_smiles(ids["smiles"])
    if "inchikey" in ids:
        head_inchi = ids["inchikey"].split('-')[0]
        ligand.clean_inchikey = head_inchi
        ligand.inchikey = ids["inchikey"]

    # perform RDkit calculations
    if lig_type == "small-molecule" and ligand.smiles is not None and ligand.smiles != "":
        input_mol = dm.to_mol(ligand.smiles, sanitize=True)
        if input_mol:
            # Check if InChIKey has been set
            if ligand.inchikey is None:
                inchi = dm.to_inchikey(input_mol)
                ligand.clean_inchikey = inchi.split('-')[0]
                ligand.inchikey = inchi

            # Calculate RDkit properties
            ligand.mw = dm.descriptors.mw(input_mol)
            ligand.rotatable_bonds = dm.descriptors.n_rotatable_bonds(input_mol)
            ligand.hacc = dm.descriptors.n_hba(input_mol)
            ligand.hdon = dm.descriptors.n_hbd(input_mol)
            ligand.logp = dm.descriptors.clogp(input_mol)

    ligand.parent = None

    # >>> NEW: assign gpcrdb_id here
    ligand.gpcrdb_id = allocate_next_gpcrdb_id()

    ligand.save()
    return ligand

def try_get_parent(query_params):
    try:
        parent_obj = Ligand.objects.get(**query_params)
        return parent_obj
    except Ligand.MultipleObjectsReturned:
        if 'name__iexact' in query_params:
            query_params['name'] = query_params.pop('name__iexact')
            try:
                parent_obj = Ligand.objects.get(**query_params)
                return parent_obj
            except Ligand.DoesNotExist:
                return None
        else:
            return None
    except Ligand.DoesNotExist:
        return None


def check_name(smiles: str,
               inchi: str,
               name: str,
               retries: int = 5,
               backoff_factor: float = 1.0):
    """
    Look up the IUPAC name for `smiles` on PubChem if any of
    (smiles, inchi, mol_weight) differs from the stored parent.
    Retries on PUGREST.ServerBusy with exponential backoff.

    Returns the new name if found, otherwise the original `name`.
    """
    # 1) Fetch parent object or bail
    try:
        parent = Ligand.objects.get(name=name, parent__isnull=True)
    except Ligand.DoesNotExist:
        return name

    # 2) Compute RDKit molecular weight (or bail)
    try:
        mol = dm.to_mol(smiles, sanitize=True)
        if mol is None:
            return name
        mol_weight = dm.descriptors.mw(mol)
    except Exception:
        return name

    # 3) If nothing has changed, return early
    if (parent.smiles == smiles and
        parent.inchikey == inchi and
        abs(parent.mw - mol_weight) < 1e-6):
        return name

    # 4) Otherwise, fetch IUPAC name with retries
    for attempt in range(retries):
        try:
            compounds = pcp.get_compounds(smiles, namespace='smiles')
            if not compounds:
                return name
            new_name = compounds[0].iupac_name
            return new_name or name
        except PubChemHTTPError as e:
            if 'PUGREST.ServerBusy' in str(e):
                wait = backoff_factor * (2 ** attempt)
                time.sleep(wait)
                continue
            # any other HTTP error: give up
            return name
        except BadRequestError:
            return name

    # 5) All retries failed
    return name

def update_parent(parent, child):
    """Fill missing structural fields on parent from child. Returns the parent.

    Collision detection / reassignment is handled upstream in the batch function.
    This function only propagates data from child to a structure-less parent.
    """
    updated = False

    # Derive correct parent clean_inchikey from child-normalized InChIKey
    correct_clean_key = child.clean_inchikey.split("-")[0] if child.clean_inchikey else None

    if correct_clean_key and not parent.clean_inchikey:
        parent.clean_inchikey = correct_clean_key
        updated = True

    if not parent.smiles and child.smiles:
        norm = normalize_for_parent(child.smiles)
        if norm:
            parent.smiles = norm["smiles"]
            parent.inchikey = norm["inchikey"]
            updated = True

    for field in ["sequence", "logp", "mw", "helm", "uniprot", "pdbe"]:
        if not getattr(parent, field) and getattr(child, field):
            setattr(parent, field, getattr(child, field))
            updated = True

    if updated:
        try:
            parent.save()
        except IntegrityError as e:
            log.error(f"update_parent failed for {parent.name!r}: {e!r}")

    return parent


#### Block for fixing mismatched LigandType assignment in Ligand model

# map LigandType IDs:
SMALL_MOLECULE = LigandType.objects.get(slug='small-molecule')
PEPTIDE        = LigandType.objects.get(slug='peptide')
PROTEIN        = LigandType.objects.get(slug='protein')
UNKNOWN        = LigandType.objects.get(slug='na')
LIPID          = LigandType.objects.get(slug='lipid')

from rdkit import Chem as _RDChem
from rdkit.Chem import Descriptors as _RDDesc

# 8+ consecutive non-ring sp3 CH2 — characteristic of fatty acid / acylglycerol / sphingolipid chains
_LIPID_CHAIN_SMARTS = _RDChem.MolFromSmarts('[CH2;!R][CH2;!R][CH2;!R][CH2;!R][CH2;!R][CH2;!R][CH2;!R][CH2;!R]')
# Steroid tetracyclic nucleus (cholesterol, bile acids, steroid hormones)
_STEROID_SMARTS = _RDChem.MolFromSmarts('[#6]1~[#6]~[#6]~[#6]2~[#6]~[#6]~[#6]3~[#6]~[#6]~[#6]~[#6]~[#6]3~[#6]~[#6]2~[#6]~1')

def predict_type(ligand):
    """
    Predict small-molecule vs peptide vs protein vs lipid from
    the presence/length of sequence *or* the pattern of amide
    bonds and lipid structural features in SMILES if sequence is missing.
    """
    seq = (ligand.sequence or '').strip()
    smiles = (ligand.smiles or '').strip()

    # 1) If we actually have a sequence, use length
    if seq:
        return PEPTIDE if len(seq) <= 15 else PROTEIN

    # 2) No sequence: look at SMILES
    if smiles:
        mol = _RDChem.MolFromSmiles(smiles)
        if mol is not None:
            frac_sp3 = _RDDesc.FractionCSP3(mol)
            n_aromatic_rings = sum(
                1 for r in mol.GetRingInfo().AtomRings()
                if all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in r)
            )
            has_long_chain = mol.HasSubstructMatch(_LIPID_CHAIN_SMARTS)
            has_steroid = mol.HasSubstructMatch(_STEROID_SMARTS)
            if (has_long_chain and frac_sp3 >= 0.65 and n_aromatic_rings <= 1) or has_steroid:
                return LIPID

        # count amide bonds (one per residue) in either orientation
        smiles_up = smiles.upper()
        count = (
            smiles_up.count('NC(=O)')
          + smiles_up.count('C(=O)N')
        )

        if count == 0:
            return SMALL_MOLECULE
        if count <= 15:
            return PEPTIDE
        return PROTEIN

    # 3) nothing to go on
    return UNKNOWN

def score(item):
    """
    Same as before: (ligand, predicted_type) → numeric score
    """
    lig, pred = item
    if pred in (PROTEIN, PEPTIDE):
        return len(lig.sequence or '')      # longer is better
    if pred == SMALL_MOLECULE:
        return -(len(lig.smiles or ''))     # shorter is better
    return float('-inf')

def pick_canonical_type(children):
    """
    Given a list of Ligand children, returns the string
    of the best ligand type (one of SMALL_MOLECULE, PEPTIDE, PROTEIN).
    """
    # 1) compute predictions
    pairs = []
    for lig in children:
        pred = predict_type(lig)
        pairs.append((lig, pred))

    # 2) split into "matched" vs "unmatched" against actual
    matched = []
    for lig, pred in pairs:
        actual = (lig.ligand_type.name.lower() if lig.ligand_type else UNKNOWN)
        if pred == actual:
            matched.append((lig, pred))

    if matched:
        candidates = matched
    else:
        # no actual matches → majority vote on predictions
        counts = Counter(pred for _, pred in pairs)
        most_common_pred, _ = counts.most_common(1)[0]
        candidates = [(lig, pred) for lig, pred in pairs if pred == most_common_pred]

    # 3) tie‐break by length heuristics
    best_lig, best_pred = max(candidates, key=score)
    return best_pred

def resolve_all_parents():
    """
    Find all parents with children of mixed types,
    and for each, choose the canonical child.
    Returns a dict: parent_id → (chosen_ligand, [all_children])
    """
    from django.db.models import Count

    # 1) find parent IDs with >1 distinct ligand_type among children
    mixed_parents = (
        Ligand.objects
              .filter(parent__isnull=False)
              .values('parent')
              .annotate(cnt=Count('ligand_type', distinct=True))
              .filter(cnt__gt=1)
              .values_list('parent', flat=True)
    )

    # 2) fetch all children for those parents
    qs = Ligand.objects.filter(parent_id__in=mixed_parents).select_related('ligand_type')

    # 3) group in Python
    by_parent = defaultdict(list)
    for child in qs:
        by_parent[child.parent_id].append(child)

    # 4) pick canonical for each
    result = {}
    for pid, kids in by_parent.items():
        result[pid] = (pick_canonical_type(kids), kids)

    return result

def apply_canonical_ligand_types():
    """
    For each parent with mixed‐type children,
    determines the canonical LigandType (via pick_canonical_child),
    then updates *all* those children—and the parent record itself—
    to use that LigandType.
    Returns the total number of rows updated.
    """
    mapping = resolve_all_parents()    # { parent_id: (best_type, [kids]) }
    total_updated = 0

    with transaction.atomic():
        for parent_id, (best_type, kids) in mapping.items():
            # collect children+p arent IDs
            ids = [c.pk for c in kids] + [parent_id]
            updated = Ligand.objects.filter(pk__in=ids).update(ligand_type=best_type)
            total_updated += updated

    return total_updated


def classify_lipids():
    """
    Reclassify small-molecule children whose SMILES match lipid structural
    patterns (long aliphatic chain or steroid nucleus), and propagate the
    lipid type to their parent records.
    Returns count of updated records.
    """
    candidates = Ligand.objects.filter(
        ligand_type=SMALL_MOLECULE,
        smiles__isnull=False,
        parent__isnull=False,
    ).select_related('parent')

    lipid_ids, parent_ids = [], set()
    for lig in candidates.iterator(chunk_size=500):
        mol = _RDChem.MolFromSmiles(lig.smiles or '')
        if mol is None:
            continue
        frac_sp3 = _RDDesc.FractionCSP3(mol)
        n_ar = sum(
            1 for r in mol.GetRingInfo().AtomRings()
            if all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in r)
        )
        has_chain = mol.HasSubstructMatch(_LIPID_CHAIN_SMARTS)
        has_steroid = mol.HasSubstructMatch(_STEROID_SMARTS)
        if (has_chain and frac_sp3 >= 0.65 and n_ar <= 1) or has_steroid:
            lipid_ids.append(lig.pk)
            parent_ids.add(lig.parent_id)

    if not lipid_ids:
        return 0

    with transaction.atomic():
        updated = Ligand.objects.filter(
            pk__in=lipid_ids + list(parent_ids)
        ).update(ligand_type=LIPID)

    return updated
