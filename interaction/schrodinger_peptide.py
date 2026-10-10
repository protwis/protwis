"""
Import Schrodinger Engine 2 receptor x peptide interactions into GPCRdb.

Engine 2 computes every pair of polymer segments a structure declares, each
pair both ways round, and its delivered tree holds one directory per structure:

    {data_dir}/{PDB}/plan.json              segments and work items
    {data_dir}/{PDB}/{key}/{key}.json       one record per work item
    {data_dir}/{PDB}/{key}/{key}.yaml       the interaction rows, when done
    {data_dir}/{PDB}/peptide_map.tsv        build_schrodinger_peptide_maps

The producer reads no GPCRdb annotation, so which work item answers which
anchor is decided on this side. build_schrodinger_peptide_maps decides it from
files only and writes peptide_map.tsv; the import reads that file and never
recomputes it:

* the receptor's author chain is the one that carries GPCRdb's receptor chain
  (matched by CA coordinates), and the receptor side is every segment on it,
  as in Engine 1, since a receptor can be declared in several pieces; rows on a
  fusion partner find no GPCRdb residue and are dropped and counted. The
  segment covering most of GPCRdb's receptor residues is recorded as the
  primary one; it must exist, or the receptor is not resolved. Sequence
  references are not used: they often name the entry itself or another
  species;
* a peptide's GPCRdb chain is matched to an author chain by CA coordinates,
  and by all atom coordinates when the chain has no standard CA (peptides of
  D-amino acids are all HETATM);
* every segment on that author chain is the peptide (a C-terminal amide cap is
  a segment of its own), and the work items are those with a peptide segment
  as the ligand side and a receptor segment as the receptor side.

The import writes two places, both replaced per anchor:

* ResidueFragmentInteraction rows on the anchor, which the structure page
  reads, routed exactly as Engine 1 routes them;
* InteractingPeptideResiduePair and InteractionPeptide rows under the anchor's
  LigandPeptideStructure, which the REST API reads, in the vocabulary the
  legacy contact network wrote (receptor first, peptide second). Families that
  table has no type for are not written there.

The anchors served are those whose pdb_reference is "pep", whatever the
ligand type: peptides, peptide drugs typed as small molecules, and protein
partners. Engine 1 skips every "pep" anchor, so the two imports never meet.
"""

import collections
import csv
import hashlib
import json
import os
import re

from django.db import transaction

from interaction import schrodinger_chain_map as chain_map
from interaction import schrodinger_complex as complex_file
from interaction import schrodinger_import as si

from contactnetwork.models import InteractingPeptideResiduePair, InteractionPeptide
from interaction.models import (
    ResidueFragmentInteraction,
    ResidueFragmentInteractionType,
    StructureLigandInteraction,
)
from ligand.models import LigandPeptideStructure
from residue.models import Residue
from structure.models import Fragment, PdbData, Rotamer


# ---------------------------------------------------------------------------
# Scope
# ---------------------------------------------------------------------------

PEPTIDE_REFERENCE = "PEP"


def is_in_scope(sli):
    """True iff this anchor is a chain this import serves: every "pep" anchor."""
    reference = (sli.pdb_reference or "").strip().upper()
    return reference == PEPTIDE_REFERENCE


# ---------------------------------------------------------------------------
# The delivered tree (pure)
# ---------------------------------------------------------------------------

PLAN_NAME = "plan.json"
MAP_NAME = "peptide_map.tsv"

DONE = "done"
# Outcomes that answer the question with "no interface": no rows, not a gap.
NO_INTERFACE = frozenset({"selections_apart", "selection_empty", "selections_overlap"})
# Terminal outcomes of a run that looked and could not answer (outcome words
# Engine 2 writes in each item's record): the anchor gets no rows and is
# reported, as Engine 1 treats a structure it ran without a product. Anything
# else that is not done (started, an unknown word) is an unfinished or
# unreadable delivery and fails the structure.
FAILED = frozenset({"preparation_failed", "compute_failed", "timed_out", "crashed"})


class MalformedProduct(si.MalformedProduct):
    """A delivered file that cannot be taken at its word."""


def read_json(path):
    try:
        with open(path) as fh:
            doc = json.load(fh)
    except (OSError, ValueError) as exc:
        raise MalformedProduct("{}: {}".format(path, exc))
    if not isinstance(doc, dict):
        raise MalformedProduct("{}: not a JSON object".format(path))
    return doc


# The contract version of the Engine 2 deliveries this importer reads, as
# plan.json and every item record carry it (the map builder checks the plan,
# anchor_rows each record it reads). The producer changes it when the meaning
# or the required fields of a product record change.
PRODUCT_CONTRACT = "engine2/3.0"


def load_plan(data_dir, pdb):
    """(segments by name, items) of one structure's plan.json."""
    path = os.path.join(data_dir, pdb, PLAN_NAME)
    plan = read_json(path)
    if plan.get("contract_version") != PRODUCT_CONTRACT:
        raise MalformedProduct(
            "{}: contract_version {!r}, this importer reads {!r}".format(
                path, plan.get("contract_version"), PRODUCT_CONTRACT
            )
        )
    segments = plan.get("segments")
    items = plan.get("items")
    if not isinstance(segments, list) or not isinstance(items, list):
        raise MalformedProduct("{}: no segments or items list".format(path))
    by_name = {}
    for seg in segments:
        if not isinstance(seg, dict) or not seg.get("name") or seg["name"] in by_name:
            raise MalformedProduct("{}: a segment without a unique name".format(path))
        by_name[seg["name"]] = seg
    for it in items:
        if not isinstance(it, dict) or not all(
            it.get(k) for k in ("key", "ligand_segment", "receptor_segment")
        ):
            raise MalformedProduct("{}: an item without key or segments".format(path))
    return by_name, items


def item_paths(data_dir, pdb, key):
    d = os.path.join(data_dir, pdb, key)
    return os.path.join(d, key + ".json"), os.path.join(d, key + ".yaml")


def segment_residues(segment):
    """Residue numbers a segment declares, from its ranges."""
    out = set()
    for pair in segment.get("ranges") or []:
        lo, hi = int(pair[0]), int(pair[1])
        out.update(range(lo, hi + 1))
    return out


# ---------------------------------------------------------------------------
# Deciding segments and items (pure; used by the map builder)
# ---------------------------------------------------------------------------


def receptor_ca_numbers(auth_chain, preferred_chain, cif_atoms, gpcrdb_atoms):
    """Author residue numbers of the receptor CA atoms GPCRdb and the mmCIF share."""
    keys = {
        a["key"]
        for a in gpcrdb_atoms
        if a["group"] == "ATOM" and a["atom"] == "CA" and a["chain"] == preferred_chain
    }
    out = set()
    for a in cif_atoms:
        if a["auth_asym"] == auth_chain and a["atom"] == "CA" and a["key"] in keys:
            try:
                out.add(int(a["auth_seq"]))
            except ValueError:
                continue
    return out


def receptor_segment(segments, auth_chain, ca_numbers):
    """
    (segment name, residues covered, note) of the receptor's main segment.

    The main segment is the one on the receptor's author chain that covers most
    of GPCRdb's receptor residues.

    None when no segment covers any; a tie is refused rather than broken.
    """
    scored = sorted(
        (
            (len(segment_residues(s) & ca_numbers), name)
            for name, s in segments.items()
            if s.get("chain_id") == auth_chain
        ),
        reverse=True,
    )
    if not scored or scored[0][0] == 0:
        return (
            None,
            0,
            "no segment on chain {} covers a receptor residue".format(auth_chain),
        )
    if len(scored) > 1 and scored[1][0] == scored[0][0]:
        return (
            None,
            scored[0][0],
            "segments {} and {} cover the receptor equally".format(
                scored[0][1], scored[1][1]
            ),
        )
    return scored[0][1], scored[0][0], ""


# A chain matched on all atoms needs this share of GPCRdb's atoms of that chain.
ANY_ATOM_MIN_SHARE = 0.9


def peptide_author_chain(pdb, gpcrdb_chain, cif_atoms, gpcrdb_atoms):
    """
    (author chain, method, note) for one GPCRdb peptide chain.

    method "ca": the receptor rule of the Engine 1 chain map, applied to the
    peptide chain; a renumbered match is still a match (the peptide numbers
    written come from the product) and is noted. method "any_atom": the chain
    has no standard CA in GPCRdb's text (a peptide of D-amino acids is all
    HETATM), so every atom of the chain there is looked up among the
    author-side atoms of the coordinate index, which keeps every HETATM record,
    so an all-HETATM chain is complete in it. A chain stored as HETATM in GPCRdb
    but as ATOM with CA atoms in the mmCIF falls under the share and is
    refused.
    """
    res = chain_map.resolve_receptor(pdb, gpcrdb_chain, cif_atoms, gpcrdb_atoms)
    if res["status"] in ("ok", "renumbered") and res["auth_chain"]:
        note = "renumbered" if res["status"] == "renumbered" else res.get("note", "")
        return (
            res["auth_chain"],
            "ca" if res["method"] == "exact" else res["method"],
            note,
        )
    keys = {a["key"] for a in gpcrdb_atoms if a["chain"] == gpcrdb_chain}
    if not keys:
        return "", "", "GPCRdb text has no atom on chain {}".format(gpcrdb_chain)
    per = collections.Counter(a["auth_asym"] for a in cif_atoms if a["key"] in keys)
    if not per:
        return (
            "",
            "",
            "no atom of GPCRdb chain {} is in the coordinate index".format(
                gpcrdb_chain
            ),
        )
    ((best, n),) = per.most_common(1)
    if n < ANY_ATOM_MIN_SHARE * len(keys) or list(per.values()).count(n) > 1:
        return (
            "",
            "",
            "GPCRdb chain {}: {} of {} atoms on author chain {}".format(
                gpcrdb_chain, n, len(keys), best
            ),
        )
    return best, "any_atom", "{}/{} atoms".format(n, len(keys))


def receptor_chain_segments(segments, auth_chain):
    """Every segment on the receptor's author chain, sorted: the receptor side."""
    return sorted(
        name for name, s in segments.items() if s.get("chain_id") == auth_chain
    )


def peptide_items(segments, items, author_chain, receptor_segs):
    """
    [(peptide segment, item key)] with the peptide as the ligand side.

    The receptor side is a receptor segment; a segment is never both.
    """
    receptor_segs = set(receptor_segs)
    peptide_segs = {
        name
        for name, s in segments.items()
        if s.get("chain_id") == author_chain and name not in receptor_segs
    }
    return sorted(
        (it["ligand_segment"], it["key"])
        for it in items
        if it["ligand_segment"] in peptide_segs
        and it["receptor_segment"] in receptor_segs
    )


# ---------------------------------------------------------------------------
# peptide_map.tsv
# ---------------------------------------------------------------------------

MAP_SCHEMA = "engine2-peptide-map/1"
RECEPTOR_PREFIX = "receptor."
RECEPTOR_KEYS = (
    "preferred_chain",
    "auth_chain",
    "status",
    "method",
    "segment",
    "segments",
    "n_ca_gpcrdb",
    "n_ca_matched",
    "n_covered",
    "gpcrdb_text_sha256",
    "note",
)
MAP_COLUMNS = (
    "pdb",
    "gpcrdb_chain",
    "titles",
    "types",
    "auth_chain",
    "chain_method",
    "status",
    "segments",
    "items",
    "outcomes",
    "note",
)
PROVENANCE_KEYS = (
    "annotation_commit",
    "ligands_sha256",
    "structures_sha256",
    "builder_sha256",
)
LIST_SEP = ","

# Row statuses. Only "ok" carries items to import.
ROW_OK = "ok"
ROW_STATUSES = frozenset(
    {ROW_OK, "chain_unresolved", "no_receptor_segment", "no_items"}
)


class MapMismatch(si.MapMismatch):
    """peptide_map.tsv does not describe what it should."""


class UnresolvedAnchor(si.UnresolvedAnchor):
    """The map could not decide which items answer an anchor."""


def load_peptide_map(path):
    """(pdb, receptor dict, {gpcrdb_chain: row}, provenance, header) from one peptide_map.tsv."""
    header, fieldnames, rows = si._read_map(path)
    if header.get("schema") != MAP_SCHEMA:
        raise MapMismatch(
            "{}: schema {!r}, this importer reads {!r}".format(
                path, header.get("schema"), MAP_SCHEMA
            )
        )
    if tuple(fieldnames or ()) != MAP_COLUMNS:
        raise MapMismatch(
            "{}: columns {}, expected {}".format(
                path, list(fieldnames or []), list(MAP_COLUMNS)
            )
        )
    pdb = (header.get("pdb") or "").strip().upper()
    if not pdb:
        raise MapMismatch("{}: the header names no pdb".format(path))
    receptor = {}
    for key in RECEPTOR_KEYS:
        if RECEPTOR_PREFIX + key not in header:
            raise MapMismatch(
                "{}: the header has no {}{}".format(path, RECEPTOR_PREFIX, key)
            )
        receptor[key] = header[RECEPTOR_PREFIX + key]
    receptor["segments_list"] = [x for x in receptor["segments"].split(LIST_SEP) if x]
    table = {}
    for r in rows:
        if (r.get("pdb") or "").strip().upper() != pdb:
            raise MapMismatch(
                "{}: a row for {!r} in the map of {}".format(path, r.get("pdb"), pdb)
            )
        chain = r.get("gpcrdb_chain") or ""
        if not chain or chain in table:
            raise MapMismatch(
                "{}: missing or duplicate gpcrdb_chain {!r}".format(path, chain)
            )
        if r.get("status") not in ROW_STATUSES:
            raise MapMismatch(
                "{}: status {!r} for chain {}".format(path, r.get("status"), chain)
            )
        r["segments_list"] = [x for x in (r.get("segments") or "").split(LIST_SEP) if x]
        r["items_list"] = [x for x in (r.get("items") or "").split(LIST_SEP) if x]
        r["outcomes_list"] = [x for x in (r.get("outcomes") or "").split(LIST_SEP) if x]
        if len(r["items_list"]) != len(r["outcomes_list"]) or (
            (r["status"] == ROW_OK) != bool(r["items_list"])
        ):
            raise MapMismatch(
                "{}: chain {}: items and outcomes do not agree with status {}".format(
                    path, chain, r["status"]
                )
            )
        table[chain] = r
    provenance = {k: header.get(k, "") for k in PROVENANCE_KEYS}
    return pdb, receptor, table, provenance, header


def write_peptide_map(path, header, receptor, rows):
    """Write one peptide_map.tsv; header and receptor are ordered key/value pairs."""
    tmp = path + ".tmp"
    with open(tmp, "w", newline="") as fh:
        fh.write("# schema\t{}\n".format(MAP_SCHEMA))
        for key, value in header:
            fh.write("# {}\t{}\n".format(key, value))
        for key in RECEPTOR_KEYS:
            fh.write("# {}{}\t{}\n".format(RECEPTOR_PREFIX, key, receptor.get(key, "")))
        writer = csv.DictWriter(
            fh,
            fieldnames=MAP_COLUMNS,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="ignore",
        )
        writer.writeheader()
        for r in rows:
            writer.writerow(r)
    os.replace(tmp, path)


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


# ---------------------------------------------------------------------------
# Peptide atom lines (pure)
# ---------------------------------------------------------------------------


class MalformedPeptideLine(MalformedProduct):
    """A peptide atom line that does not have the producer's layout."""


# Residue names the structure preparation gives protonation and disulfide
# states; GPCRdb's own text, and the legacy peptide table, use the standard one.
PREPARED_NAMES = {
    "HID": "HIS",
    "HIE": "HIS",
    "HIP": "HIS",
    "CYX": "CYS",
    "ASH": "ASP",
    "GLH": "GLU",
    "LYN": "LYS",
    "ARN": "ARG",
}


def parse_peptide_line(line, product_chain):
    """
    One producer atom line of a peptide (or any multi-residue ligand side).

    The producer's layout is one column short of standard PDB (no altloc) and
    widens the residue name for five-character codes; the chain it writes is
    the product (author) chain, which is known, so the fields are read around
    it. Returns a dict, or raises MalformedPeptideLine.
    """
    record = line[:6].strip()
    if record not in ("ATOM", "HETATM"):
        raise MalformedPeptideLine("not an atom record: {!r}".format(line[:30]))
    serial = line[6:11].strip()
    name = line[12:16].strip()
    m = re.match(
        r"^(?P<resname>[A-Za-z0-9 ]{3,5}?) "
        + re.escape(product_chain)
        + r"(?P<num>[ -]*-?\d+)(?P<icode>[A-Za-z ]?)(?P<tail>.*)$",
        line[16:],
    )
    tail = si._LIGAND_TAIL_RE.match(m.group("tail")) if m else None
    if not serial.isdigit() or not name or m is None or tail is None:
        raise MalformedPeptideLine(
            "chain {}: cannot read {!r}".format(product_chain, line)
        )
    return {
        "record": record,
        "serial": int(serial),
        "name": name,
        "resname": PREPARED_NAMES.get(
            m.group("resname").strip().upper(), m.group("resname").strip().upper()
        ),
        "resnum": int(m.group("num").replace(" ", "")),
        "icode": m.group("icode").strip(),
        "x": float(tail.group("x")),
        "y": float(tail.group("y")),
        "z": float(tail.group("z")),
        "occ": float(tail.group("occ")),
        "b": float(tail.group("b")),
        "element": tail.group("element").upper(),
    }


def standard_peptide_line(atom, gpcrdb_chain):
    """
    (standard PDB v3.3 line, b_factor_was_capped) for one parsed atom.

    The residue name is cut to three characters and the chain is GPCRdb's, as
    in GPCRdb's own stored structure text.
    """
    capped = atom["b"] > si._MAX_PDB_B
    out = "{:<6}{:>5} {} {:>3} {:1}{:>4}{:1}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2}".format(
        atom["record"],
        atom["serial"] % 100000,
        si._pdb_atom_name(atom["name"], atom["element"]),
        atom["resname"][:3],
        gpcrdb_chain,
        atom["resnum"],
        atom["icode"],
        atom["x"],
        atom["y"],
        atom["z"],
        atom["occ"],
        min(atom["b"], si._MAX_PDB_B),
        atom["element"],
    )
    if len(out) != si._PDB_ATOM_LINE_WIDTH or len(gpcrdb_chain) != 1:
        raise MalformedPeptideLine(
            "a field does not fit standard PDB columns: {!r}".format(atom)
        )
    return out, capped


def peptide_atoms(block, product_chain):
    """Parsed atoms of one row's ligand_pdb_block."""
    return [
        parse_peptide_line(line, product_chain)
        for line in (block or "").splitlines()
        if line.strip()
    ]


# ---------------------------------------------------------------------------
# Rows for the peptide tables (pure)
# ---------------------------------------------------------------------------

# (family, direction) -> (interaction_type, specific_type) in the legacy
# contact-network vocabulary, which names the receptor side first: a peptide
# donor is "acceptor-donor", a negative peptide "positive-negative".
PEPTIDE_TYPES = {
    ("Donor", "ligand-donor"): ("polar", "acceptor-donor"),
    ("Acceptor", "ligand-acceptor"): ("polar", "donor-acceptor"),
    ("NegCharge", "neg-pos"): ("ionic", "positive-negative"),
    ("PosCharge", "pos-neg"): ("ionic", "negative-positive"),
    ("PiCat", "ligand-cation"): ("aromatic", "pi-cation"),
    ("PiCat", "receptor-cation"): ("aromatic", "cation-pi"),
    ("Aromatic", "face-to-face"): ("aromatic", "face-to-face"),
    ("Aromatic", "edge-to-face"): ("aromatic", "edge-to-face"),
    ("HPhob", ""): ("hydrophobic", ""),
    ("VdW", ""): ("van-der-waals", ""),
}

# Families the peptide tables have no type for. They still go to RFI, except
# Wat-HBond, which goes nowhere (si.EXCLUDED_FAMILIES).
NOT_IN_PEPTIDE_TABLES = frozenset(
    {"Accessible", "Covalent", "XBond", "Metal", "Wat-HBond"}
)

# The legacy contact network names an aromatic ring by this pseudo-atom.
RING_ATOM = "RN1"

# The legacy level: 0 is the normal definition. Engine 2 has no loosened one.
LEVEL = 0

ONE_LETTER = {
    "ALA": "A",
    "ARG": "R",
    "ASN": "N",
    "ASP": "D",
    "CYS": "C",
    "GLN": "Q",
    "GLU": "E",
    "GLY": "G",
    "HIS": "H",
    "ILE": "I",
    "LEU": "L",
    "LYS": "K",
    "MET": "M",
    "PHE": "F",
    "PRO": "P",
    "SER": "S",
    "THR": "T",
    "TRP": "W",
    "TYR": "Y",
    "VAL": "V",
}


class UnroutableRow(si.UnroutableRow):
    """A (family, direction) the peptide vocabulary does not know."""


def peptide_type(row):
    """
    (interaction_type, specific_type, receptor_is_ring, peptide_is_ring).

    None for a family the peptide tables do not hold.
    """
    family = row["feature_family"]
    if family in NOT_IN_PEPTIDE_TABLES:
        return None
    key = (family, row.get("direction") or "")
    if key not in PEPTIDE_TYPES:
        raise UnroutableRow("no peptide table type for {}".format(key))
    itype, detail = PEPTIDE_TYPES[key]
    if family == "Aromatic":
        return itype, detail, True, True
    if family == "PiCat":
        receptor_ring = key[1] == "ligand-cation"
        return itype, detail, receptor_ring, not receptor_ring
    return itype, detail, False, False


def plan_peptide_pairs(rows, receptor_chain_name, product_chain):
    """
    Plan the peptide-table rows for the product rows of one anchor.

    Returns (pairs, counts). ``pairs`` maps
    (peptide resnum, icode, resname, receptor seq, receptor one-letter) to a
    sorted list of distinct (peptide_atom, receptor_atom, interaction_type,
    specific_type). counts accounts for every input row::

        rows_in == not_in_peptide_tables + nonstandard_residue + other_chain + used

    A peptide atom whose residue has an insertion code (antibody numbering) is
    left out and counted as insertion_code_atoms: the table has no column for
    the code, and writing the bare number would name another residue. The RFI
    fragment text keeps the code.
    """
    counts = collections.Counter()
    counts["rows_in"] = len(rows)
    pairs = collections.defaultdict(set)
    for row in rows:
        typed = peptide_type(row)
        if typed is None:
            counts["not_in_peptide_tables"] += 1
            continue
        res = row["receptor_residue"]
        amino_acid = str(res["name_1_letter"]).upper()
        if amino_acid not in si.STANDARD_AMINO_ACIDS:
            counts["nonstandard_residue"] += 1
            continue
        if receptor_chain_name and str(res["chain_id"]) != receptor_chain_name:
            counts["other_chain"] += 1
            continue
        counts["used"] += 1
        itype, detail, receptor_ring, peptide_ring = typed
        receptor_atom = (
            RING_ATOM if receptor_ring else str(row.get("receptor_atom_name") or "")
        )
        seq = int(res["pdb_residue_number"])
        for atom in peptide_atoms(row.get("ligand_pdb_block"), product_chain):
            if atom["icode"]:
                counts["insertion_code_atoms"] += 1
                continue
            key = (atom["resnum"], atom["icode"], atom["resname"], seq, amino_acid)
            pairs[key].add(
                (
                    RING_ATOM if peptide_ring else atom["name"],
                    receptor_atom,
                    itype,
                    detail,
                )
            )
    out = {k: sorted(v) for k, v in pairs.items()}
    counts["pairs"] = len(out)
    counts["interactions"] = sum(len(v) for v in out.values())
    return out, counts


def standardise_blocks(rows, product_chain, gpcrdb_chain):
    """
    Rewrite every row's ligand_pdb_block in standard PDB columns, in place.

    Returns the set of output lines whose B-factor was capped.
    """
    capped = set()
    for row in rows:
        lines = []
        for atom in peptide_atoms(row.get("ligand_pdb_block"), product_chain):
            line, was_capped = standard_peptide_line(atom, gpcrdb_chain)
            lines.append(line)
            if was_capped:
                capped.add(line)
        row["ligand_pdb_block"] = "\n".join(lines)
    return capped


def anchor_rows(data_dir, pdb, map_row):
    """
    (rows, failed) of one anchor, from the items its map row names.

    Each item's record is read again and must still say what the map says.
    ``failed`` lists "key:outcome" for the items the run could not answer
    (FAILED); when it is not empty ``rows`` is empty too, so an anchor is
    never half imported. Otherwise ``rows`` are the rows of the done items,
    each of which must have its YAML. An item that is neither done, a
    no-interface answer nor failed leaves the question open, and raises.
    """
    rows, failed = [], []
    for key, expected in zip(map_row["items_list"], map_row["outcomes_list"]):
        rec_path, yaml_path = item_paths(data_dir, pdb, key)
        record = read_json(rec_path)
        provenance = record.get("provenance")
        version = (
            provenance.get("contract_version") if isinstance(provenance, dict) else None
        )
        if version != PRODUCT_CONTRACT:
            raise MalformedProduct(
                "{}: contract_version {!r}, this importer reads {!r}".format(
                    rec_path, version, PRODUCT_CONTRACT
                )
            )
        if record.get("work_item_key") != key or record.get("outcome") != expected:
            raise MalformedProduct(
                "{}: record says {}/{}, the map says {}/{}".format(
                    rec_path,
                    record.get("work_item_key"),
                    record.get("outcome"),
                    key,
                    expected,
                )
            )
        if expected == DONE:
            rows.extend(si.read_instance_rows(yaml_path))
        elif expected in FAILED:
            failed.append("{}:{}".format(key, expected))
        elif expected not in NO_INTERFACE:
            raise MalformedProduct(
                "{}: outcome {} leaves the question open".format(rec_path, expected)
            )
    return ([] if failed else rows), failed


# ---------------------------------------------------------------------------
# Import of one structure (database)
# ---------------------------------------------------------------------------


class AnchorOutcome(object):
    """What happened to one in-scope anchor."""

    __slots__ = (
        "sli_id",
        "chain",
        "mode",
        "items",
        "notes",
        "rfi_deleted",
        "rfi_written",
        "rfi_counts",
        "rfi_dropped",
        "fragments_created",
        "pairs_deleted",
        "interactions_deleted",
        "pairs_written",
        "interactions_written",
        "pair_counts",
        "pair_dropped",
        "complex_file",
    )

    def __init__(self, sli_id, chain):
        """Start with nothing done for this anchor."""
        self.sli_id = sli_id
        self.chain = chain
        self.mode = ""
        self.items = []
        self.notes = []
        self.rfi_deleted = 0
        self.rfi_written = 0
        self.rfi_counts = collections.Counter()
        self.rfi_dropped = collections.Counter()
        self.fragments_created = 0
        self.pairs_deleted = 0
        self.interactions_deleted = 0
        self.pairs_written = 0
        self.interactions_written = 0
        self.pair_counts = collections.Counter()
        self.pair_dropped = collections.Counter()
        self.complex_file = ""


class MissingPeptideStructure(ValueError):
    """An anchor has no LigandPeptideStructure of its own."""


def _receptor_residue(structure, seq, amino_acid, dropped):
    residues = list(
        Residue.objects.filter(
            protein_conformation=structure.protein_conformation, sequence_number=seq
        )[:2]
    )
    if not residues:
        dropped["residue_not_found"] += 1
        return None
    if len(residues) > 1:
        dropped["residue_ambiguous"] += 1
        return None
    if residues[0].amino_acid != amino_acid:
        dropped["amino_acid_mismatch"] += 1
        return None
    return residues[0]


def _write_rfi(structure, sli, records, types, outcome):
    """Write the RFI rows; returns the receptor residue numbers written."""
    written = set()
    for rec in records:
        residue = _receptor_residue(
            structure, rec["sequence_number"], rec["amino_acid"], outcome.rfi_dropped
        )
        if residue is None:
            continue
        rotamers = list(
            Rotamer.objects.filter(structure=structure, residue=residue)[:2]
        )
        if len(rotamers) != 1:
            outcome.rfi_dropped[
                "rotamer_not_found" if not rotamers else "rotamer_ambiguous"
            ] += 1
            continue
        text = si.fragment_text(rec["ligand_lines"])
        fragment = (
            Fragment.objects.filter(
                ligand=sli.ligand,
                structure=structure,
                residue=residue,
                pdbdata__pdb=text,
            )
            .order_by("id")
            .first()
        )
        if fragment is None:
            fragment = Fragment.objects.create(
                ligand=sli.ligand,
                structure=structure,
                residue=residue,
                pdbdata=PdbData.objects.create(pdb=text),
            )
            outcome.fragments_created += 1
        ResidueFragmentInteraction.objects.create(
            structure_ligand_pair=sli,
            rotamer=rotamers[0],
            fragment=fragment,
            interaction_type=types[rec["slug"]],
        )
        outcome.rfi_written += 1
        written.add(rec["sequence_number"])
    return written


def _write_pairs(structure, peptide, pairs, outcome):
    bulk = []
    for (resnum, _icode, resname, seq, amino_acid), interactions in sorted(
        pairs.items()
    ):
        residue = _receptor_residue(structure, seq, amino_acid, outcome.pair_dropped)
        if residue is None:
            continue
        pair = InteractingPeptideResiduePair.objects.create(
            peptide_amino_acid_three_letter=resname[:3],
            peptide_amino_acid=ONE_LETTER.get(resname, "X"),
            peptide_sequence_number=resnum,
            peptide=peptide,
            receptor_residue=residue,
        )
        outcome.pairs_written += 1
        for peptide_atom, receptor_atom, itype, detail in interactions:
            bulk.append(
                InteractionPeptide(
                    interacting_peptide_pair=pair,
                    peptide_atom=peptide_atom[:10],
                    receptor_atom=receptor_atom[:10],
                    interaction_type=itype,
                    specific_type=detail,
                    interaction_level=LEVEL,
                )
            )
    InteractionPeptide.objects.bulk_create(bulk)
    outcome.interactions_written += len(bulk)


def _clear_anchor(sli, peptide, outcome):
    _, deleted = ResidueFragmentInteraction.objects.filter(
        structure_ligand_pair=sli
    ).delete()
    si._only_deleted(deleted, {"interaction.ResidueFragmentInteraction"})
    outcome.rfi_deleted = deleted.get("interaction.ResidueFragmentInteraction", 0)
    _, deleted = InteractingPeptideResiduePair.objects.filter(peptide=peptide).delete()
    si._only_deleted(
        deleted,
        {
            "contactnetwork.InteractingPeptideResiduePair",
            "contactnetwork.InteractionPeptide",
        },
    )
    outcome.pairs_deleted = deleted.get(
        "contactnetwork.InteractingPeptideResiduePair", 0
    )
    outcome.interactions_deleted = deleted.get("contactnetwork.InteractionPeptide", 0)


def choose_peptide_structure(candidates, chain, label):
    """
    The anchor's LigandPeptideStructure among those of its ligand.

    The only one, or the one on the anchor's chain when the ligand has several.
    """
    found = list(candidates)
    if len(found) > 1:
        found = [p for p in found if p.chain == chain]
    if len(found) != 1:
        raise MissingPeptideStructure(
            "{}: {} LigandPeptideStructure rows".format(label, len(found))
        )
    return found[0]


def peptide_structure(structure, sli):
    return choose_peptide_structure(
        LigandPeptideStructure.objects.filter(structure=structure, ligand=sli.ligand),
        (sli.chain_res or "").strip(),
        "{} anchor {}".format(structure.pdb_code.index, sli.id),
    )


def check_receptor(pdb, receptor):
    """Refuse a map whose receptor is not resolved to a primary segment it lists."""
    if (
        receptor.get("status") != "ok"
        or not receptor.get("segment")
        or receptor["segment"] not in receptor.get("segments_list", [])
    ):
        raise MapMismatch(
            "{}: receptor not resolved ({} {})".format(
                pdb, receptor.get("status"), receptor.get("note")
            )
        )


def anchor_action(pdb, sli_id, chain, chain_rows):
    """
    (action, map row) for one anchor, before anything is written.

    "clear" for an anchor with no chain_res (it keeps no rows);
    "import" for an anchor whose chain has an ok row. A chain the map does
    not list, or a row that is not ok, raises: the map could not answer it.
    """
    if not chain:
        return "clear", None
    row = chain_rows.get(chain)
    if row is None:
        raise MapMismatch(
            "{}: the map has no row for anchor {} chain {!r}".format(pdb, sli_id, chain)
        )
    if row["status"] != ROW_OK:
        raise UnresolvedAnchor(
            "{}: anchor {} chain {}: {}: {}".format(
                pdb, sli_id, chain, row["status"], row.get("note", "")
            )
        )
    return "import", row


def claim_peptide_structure(used, peptide_id, sli_id, pdb):
    """
    Record that an anchor owns a LigandPeptideStructure.

    A second owner is refused: its clear would wipe the first anchor's pairs.
    """
    if peptide_id in used:
        raise MissingPeptideStructure(
            "{}: anchors {} and {} share LigandPeptideStructure {}".format(
                pdb, used[peptide_id], sli_id, peptide_id
            )
        )
    used[peptide_id] = sli_id


def check_receptor_and_fingerprints(pdb, data_dir, receptor, header, gpcrdb_text):
    """
    check_receptor and check_fingerprints, a stale map reported first.

    The rule of schrodinger_import.checked_receptor_chain: a map built from
    another stored text or plan is reported as stale whatever its receptor row
    says, unless the row has no text fingerprint (its build failed before the
    text was read), when its note is the reason. An ok row without one is
    refused as stale.
    """
    if not receptor.get("gpcrdb_text_sha256"):
        check_receptor(pdb, receptor)
    check_fingerprints(pdb, data_dir, receptor, header, gpcrdb_text)
    check_receptor(pdb, receptor)


def check_fingerprints(pdb, data_dir, receptor, header, gpcrdb_text):
    """Refuse a structure whose stored text or plan changed since the map was built."""
    if receptor.get("gpcrdb_text_sha256") != chain_map.text_sha256(gpcrdb_text):
        raise MapMismatch(
            "{}: GPCRdb structure text differs from the one the map was built "
            "from (new dump?); rebuild the maps".format(pdb)
        )
    plan = os.path.join(data_dir, pdb, PLAN_NAME)
    if not os.path.isfile(plan) or header.get("plan_sha256") != sha256_file(plan):
        raise MapMismatch(
            "{}: plan.json differs from the one the map was built from; check "
            "--data-dir or rebuild the maps".format(pdb)
        )


def import_structure(structure, data_dir, receptor, chain_rows, header):
    """
    Replace the peptide-import rows of one structure, in one transaction.

    ``receptor``, ``chain_rows`` and ``header`` are this structure's
    peptide_map.tsv as load_peptide_map returns them. Returns
    (outcomes, cleanup_counter). Any exception leaves the structure exactly as
    it was.

    An anchor with no chain_res names no peptide chain at all; it is cleared
    and reported. An anchor with an item the run failed gets no
    rows and is reported (mode no_product), as Engine 1 treats a structure it
    ran without a product. An anchor whose map row is anything but ok is a
    question the map could not answer, and fails the structure (as Engine 1's
    unresolved anchors do).
    """
    pdb = structure.pdb_code.index.upper()
    types = {t.slug: t for t in ResidueFragmentInteractionType.objects.all()}
    outcomes = []
    replaced_files = set()
    gpcrdb_text = structure.pdb_data.pdb if structure.pdb_data_id else ""
    with transaction.atomic():
        slis = [
            sli
            for sli in (
                StructureLigandInteraction.objects.filter(structure=structure)
                .select_related("ligand__ligand_type")
                .order_by("id")
            )
            if is_in_scope(sli)
        ]
        if not slis:
            return outcomes, collections.Counter()
        check_receptor_and_fingerprints(
            pdb,
            data_dir,
            receptor,
            header,
            structure.pdb_data.pdb if structure.pdb_data_id else "",
        )
        used = {}
        for sli in slis:
            chain = (sli.chain_res or "").strip()
            outcome = AnchorOutcome(sli.id, chain)
            peptide = peptide_structure(structure, sli)
            claim_peptide_structure(used, peptide.id, sli.id, pdb)
            action, row = anchor_action(pdb, sli.id, chain, chain_rows)
            _clear_anchor(sli, peptide, outcome)
            if action == "clear":
                # It names no peptide chain, so nothing can answer it, and an
                # in-scope anchor keeps no rows it was not given.
                outcome.mode = "cleared"
                outcome.notes.append("anchor names no chain")
                outcome.complex_file, replaced = complex_file.write_complex_file(
                    sli, ""
                )
                replaced_files.add(replaced)
                outcomes.append(outcome)
                continue
            outcome.items = list(row["items_list"])
            rows, failed = anchor_rows(data_dir, pdb, row)
            if failed:
                outcome.mode = "no_product"
                outcome.notes.append("items failed: " + LIST_SEP.join(failed))
                outcome.complex_file, replaced = complex_file.write_complex_file(
                    sli, ""
                )
                replaced_files.add(replaced)
                outcomes.append(outcome)
                continue
            outcome.mode = "imported"
            # Both tables are planned before either is written: a row the
            # vocabulary cannot route fails the structure, not half of it. The
            # peptide pairs read the producer's lines, so they come first.
            pairs, pair_counts = plan_peptide_pairs(
                rows, receptor["auth_chain"], row["auth_chain"]
            )
            outcome.pair_counts.update(pair_counts)
            capped = standardise_blocks(rows, row["auth_chain"], chain)
            outcome.rfi_counts["ligand_lines_bfactor_capped"] = len(capped)
            records, rfi_counts, _ = si.plan_rows(rows, receptor["auth_chain"])
            outcome.rfi_counts.update(rfi_counts)
            written_seqs = _write_rfi(structure, sli, records, types, outcome)
            _write_pairs(structure, peptide, pairs, outcome)
            # The anchor's 3D file: the peptide chain and the residues just written.
            text = ""
            if outcome.rfi_written:
                text = complex_file.complex_text(
                    gpcrdb_text,
                    receptor["preferred_chain"],
                    written_seqs,
                    ligand_chain=chain,
                )
            outcome.complex_file, replaced = complex_file.write_complex_file(sli, text)
            if outcome.rfi_written and not text:
                outcome.complex_file = "ligand_not_found"
            replaced_files.add(replaced)
            outcomes.append(outcome)
        cleanup = si.delete_orphan_fragments(structure)
        cleanup.update(si.delete_unreferenced_pdbdata(replaced_files))
    return outcomes, cleanup
