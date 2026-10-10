"""
Stored interactions of a structure GPCRdb has, for the online calculation page.

interaction.views.calculate computes interactions for a PDB file a user
uploads or a PDB code a user types. For a typed code that GPCRdb already has,
the answer is the interactions stored in the database, so the page agrees with
the structure page. An uploaded file is still calculated: it may not be the
deposited structure.

The results come in the shape interaction.views.calculate_interactions returns,
so the page needs no other change.
"""

import collections

from interaction.models import ResidueFragmentInteraction
from structure.models import Structure

PEPTIDE_REFERENCE = "PEP"

THREE_LETTER = {
    "A": "ALA",
    "R": "ARG",
    "N": "ASN",
    "D": "ASP",
    "C": "CYS",
    "Q": "GLN",
    "E": "GLU",
    "G": "GLY",
    "H": "HIS",
    "I": "ILE",
    "L": "LEU",
    "K": "LYS",
    "M": "MET",
    "F": "PHE",
    "P": "PRO",
    "S": "SER",
    "T": "THR",
    "W": "TRP",
    "Y": "TYR",
    "V": "VAL",
}


HIDDEN_TYPE = "hidden"


def ligand_key(pdb_reference, ligand_name):
    """(results key, is a "pep" chain) of an anchor: its HET code, or the ligand's name."""
    reference = (pdb_reference or "").strip().upper()
    if reference == PEPTIDE_REFERENCE:
        return ligand_name, True
    return reference, False


def anchor_keys(anchors):
    """
    {anchor id: (results key, is pep)} for (id, pdb_reference, ligand name, chain_res).

    Anchors that would share a key (two copies of one HET code at different
    sites) each get their chain_res appended, so one site is never built from
    two pockets.
    """
    base = {sid: ligand_key(ref, name) for sid, ref, name, _chain in anchors}
    count = collections.Counter(base.values())
    out = {}
    for sid, _ref, _name, chain_res in anchors:
        key, is_pep = base[sid]
        if count[(key, is_pep)] > 1:
            key = "{} {}".format(key, (chain_res or "").strip() or sid)
        out[sid] = (key, is_pep)
    return out


def build_results(rows, chain):
    """
    Results from database rows, in the shape calculate_interactions returns.

    The rows are ((ligand key, is pep), one-letter amino acid, residue number,
    slug, name, type, direction).

    Each interaction is [residue, fragment file, slug, name, type, direction]
    with residue as three-letter name, number and chain (ASP113A, what
    interaction.views.regexaa reads); the fragment file is left empty, the page
    does not read it. A residue that is not a standard amino acid is left out.
    calculate takes the first ligand as the main one, and the page opens on a
    HET ligand when there is one, so HET ligands come first, then "pep"
    chains, each by the number of visible rows, then by key. The score is that
    number.
    """
    per = collections.OrderedDict()
    for (key, is_pep), amino_acid, number, slug, name, type_, direction in rows:
        three = THREE_LETTER.get((amino_acid or "").upper())
        if three is None:
            continue
        per.setdefault((key, is_pep), []).append(
            [
                "{}{}{}".format(three, number, chain),
                "",
                slug,
                name,
                type_ or "",
                direction or "",
            ]
        )

    def visible(interactions):
        return sum(1 for i in interactions if i[4] != HIDDEN_TYPE)

    ordered = sorted(per.items(), key=lambda kv: (kv[0][1], -visible(kv[1]), kv[0][0]))
    return collections.OrderedDict(
        (key, {"score": visible(interactions), "interactions": interactions})
        for (key, _is_pep), interactions in ordered
    )


def stored_results(pdbname):
    """(results, stored structure text) for an experimental structure GPCRdb has, else None."""
    structure = (
        Structure.objects.filter(
            pdb_code__index__iexact=(pdbname or "").strip(),
            structure_type__origin="experiment",
        )
        .select_related("pdb_data")
        .first()
    )
    if structure is None or structure.pdb_data is None:
        return None
    chain = (structure.preferred_chain or "").split(",")[0].strip()
    rows = list(
        ResidueFragmentInteraction.objects.filter(
            structure_ligand_pair__structure=structure
        )
        .order_by(
            "structure_ligand_pair_id",
            "rotamer__residue__sequence_number",
            "interaction_type__slug",
            "id",
        )
        .values_list(
            "structure_ligand_pair_id",
            "structure_ligand_pair__pdb_reference",
            "structure_ligand_pair__ligand__name",
            "structure_ligand_pair__chain_res",
            "rotamer__residue__amino_acid",
            "rotamer__residue__sequence_number",
            "interaction_type__slug",
            "interaction_type__name",
            "interaction_type__type",
            "interaction_type__direction",
        )
    )
    keys = anchor_keys({(r[0], r[1], r[2], r[3]) for r in rows})
    return (
        build_results(((keys[r[0]],) + tuple(r[4:]) for r in rows), chain),
        structure.pdb_data.pdb,
    )
