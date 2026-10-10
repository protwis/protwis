"""
Reconcile chain names between Schrodinger products and GPCRdb.

The products name chains as the RCSB mmCIF does (author chain, up to four
characters). GPCRdb names them as its stored PDB-format structure text does
(one character, sometimes renamed, split or hand-edited by curators). This
module decides two kinds of rows, both written into one chainmap.tsv per
structure, from GPCRdb's text on one side and, on the other, the coordinate
index the producer delivers with the products (the atoms of the mmCIF they
were computed from, see below):

* anchor map: one row per GPCRdb ligand-anchor copy (pdb, HET, chain_res
  token) naming the product instance that is the same ligand copy;
* receptor map: one row per structure naming the product (author) chain whose
  atoms are GPCRdb's preferred chain.

Exact coordinate identity decides. The annotation's label_asym_id is a second,
independent witness: it confirms coordinate answers, flags annotation errors
when it disagrees, and is the only thing that can confirm a name-based
fallback where GPCRdb stores an older model whose coordinates no longer match.

No database access and no Django imports: the commands that build the maps
(build_schrodinger_chainmap_files and build_schrodinger_peptide_maps) feed this
module with text.
"""

import collections
import hashlib
import os
import re

# ---------------------------------------------------------------------------
# Coordinate index (author side)
# ---------------------------------------------------------------------------
#
# The producer delivers, beside each structure's products, the atoms of the
# mmCIF it read that this matching uses -- every HETATM record, every CA and
# every atom of a residue without a CA, from the first model, without hydrogen
# and water -- with the sha256 of that mmCIF. The producer writes the file and
# this side only reads it, so the rules for which atoms an mmCIF yields live in
# one place, next to the mmCIF itself.

WATER = frozenset({"HOH", "DOD", "WAT"})

INDEX_SCHEMA = "structure-index/1"
INDEX_SUFFIX = "_structure_index.tsv"
INDEX_COLUMNS = (
    "label_asym",
    "auth_asym",
    "comp",
    "auth_seq",
    "icode",
    "atom",
    "group",
    "x",
    "y",
    "z",
)
_SHA256 = re.compile(r"^[0-9a-f]{64}$")
# Exactly the producer's "%.3f": anything else is a producer regression, and a
# coordinate that does not match GPCRdb's text to the digit sends the matching
# down its name-based fallbacks without a sound.
_COORD = re.compile(r"^-?[0-9]+\.[0-9]{3}$")
_GROUPS = frozenset({"ATOM", "HETATM"})


class ParseError(ValueError):
    """An input file does not have the shape this module relies on."""


def coord_key(x, y, z):
    """Coordinates as a hashable key at the 3-decimal precision both formats use."""
    return "%.3f %.3f %.3f" % (float(x), float(y), float(z))


def index_path(index_dir, pdb):
    """Where the producer puts the coordinate index of ``pdb``."""
    return os.path.join(index_dir, pdb, pdb + INDEX_SUFFIX)


def parse_structure_index(text):
    """
    (cif_sha256, atoms) of a coordinate index.

    Each atom is a dict with label_asym, auth_asym, comp, auth_seq, icode, atom,
    group, key -- the shape the resolvers below take. The header must name this
    schema, a sha256 and the number of atoms that follow; the column line must be
    INDEX_COLUMNS; every row must have one value per column, a record type of
    ATOM or HETATM and coordinates with three decimals. Anything else is a
    ParseError, because a half-read or reformatted index would resolve to the
    wrong chain, or through a fallback, without a sound.
    """
    header = {}
    lines = text.splitlines()
    i = 0
    while i < len(lines) and lines[i].startswith("# "):
        key, _, value = lines[i][2:].partition("\t")
        header[key] = value
        i += 1
    if header.get("schema") != INDEX_SCHEMA:
        raise ParseError(
            "index schema {!r}, this reader reads {!r}".format(
                header.get("schema"), INDEX_SCHEMA
            )
        )
    sha = header.get("cif_sha256", "")
    if not _SHA256.match(sha):
        raise ParseError("index header has no cif_sha256")
    if not re.match(r"^[0-9]+$", header.get("atoms", "")):
        raise ParseError("index header has no atom count")
    if i >= len(lines) or tuple(lines[i].split("\t")) != INDEX_COLUMNS:
        raise ParseError("index columns are not {}".format(INDEX_COLUMNS))
    atoms = []
    for line in lines[i + 1 :]:
        fields = line.split("\t")
        if len(fields) != len(INDEX_COLUMNS):
            raise ParseError(
                "index row has {} fields, expected {}".format(
                    len(fields), len(INDEX_COLUMNS)
                )
            )
        row = dict(zip(INDEX_COLUMNS, fields))
        if not all(_COORD.match(row[c]) for c in ("x", "y", "z")):
            raise ParseError(
                "index row has a coordinate not written as -?d.ddd: {!r}".format(line)
            )
        if row["group"] not in _GROUPS:
            raise ParseError("index row has record type {!r}".format(row["group"]))
        key = coord_key(row["x"], row["y"], row["z"])
        atoms.append(
            {
                "label_asym": row["label_asym"],
                "auth_asym": row["auth_asym"],
                "comp": row["comp"],
                "auth_seq": row["auth_seq"],
                "icode": row["icode"],
                "atom": row["atom"],
                "group": row["group"],
                "key": key,
            }
        )
    if not atoms:
        raise ParseError("index lists no atom")
    if len(atoms) != int(header["atoms"]):
        raise ParseError(
            "index lists {} atoms, its header says {}: a truncated file".format(
                len(atoms), header["atoms"]
            )
        )
    return sha, atoms


# ---------------------------------------------------------------------------
# GPCRdb stored PDB-format text (GPCRdb side)
# ---------------------------------------------------------------------------


def parse_gpcrdb_pdb(text):
    """
    Return first-model, non-hydrogen, non-water atoms of GPCRdb's stored text.

    Fields are read exactly as build_structures reads them: chain = column 22
    (line[21]), residue number = columns 23-26, residue name = columns 18-20
    (so five-character CCD codes appear truncated to three characters).
    """
    atoms = []
    for line in text.splitlines():
        if line.startswith("ENDMDL"):
            break
        if not line.startswith(("ATOM", "HETATM")) or len(line) < 54:
            continue
        resname = line[17:20].strip()
        element = line[76:78].strip() if len(line) >= 78 else ""
        atom = line[12:16].strip()
        if (
            resname in WATER
            or element in ("H", "D")
            or (not element and atom.startswith("H"))
        ):
            continue
        atoms.append(
            {
                "chain": line[21],
                "resnum": line[22:26].strip(),
                "icode": line[26].strip() if len(line) > 26 else "",
                "resname": resname,
                "atom": atom,
                "group": line[:6].strip(),
                "key": coord_key(line[30:38], line[38:46], line[46:54]),
            }
        )
    return atoms


# ---------------------------------------------------------------------------
# Annotation (upstream ligands.tsv) -> label_asym_id per copy
# ---------------------------------------------------------------------------


def annotation_labels(rows):
    """
    Map (PDB, HET, token) -> label_asym_id from ligands.tsv rows.

    Residue_seq_id and label_asym_id are comma lists aligned copy for copy.
    Rows whose two lists differ in length contribute nothing.
    """
    out = {}
    for r in rows:
        tokens = [
            t.strip() for t in (r.get("Residue_seq_id") or "").split(",") if t.strip()
        ]
        labels = [
            t.strip() for t in (r.get("label_asym_id") or "").split(",") if t.strip()
        ]
        if not tokens or len(tokens) != len(labels):
            continue
        key = (
            (r.get("PDB") or "").strip().upper(),
            (r.get("Name") or "").strip().upper(),
        )
        for tok, lab in zip(tokens, labels):
            out[key + (tok.replace(" ", ""),)] = lab
    return out


# ---------------------------------------------------------------------------
# Anchor map
# ---------------------------------------------------------------------------

TOKEN_RE = re.compile(r"^(?P<chain>[A-Za-z0-9]):(?P<resnum>-?\d+)(?P<icode>[A-Za-z]?)$")

ANCHOR_COLUMNS = (
    "pdb",
    "het",
    "token",
    "instance",
    "status",
    "source",
    "n_exact",
    "n_gpcrdb",
    "label",
    "label_instance",
    "note",
)


def instance_name(het, auth_asym, auth_seq, icode):
    return "{}_{}_{}{}".format(het, auth_asym, auth_seq, icode)


def split_tokens(chain_res):
    """chain_res -> list of 'C:NUM' tokens, or [] when empty / chain only."""
    items = [t.strip() for t in (chain_res or "").split(",") if t.strip()]
    if items and all(TOKEN_RE.match(t) for t in items):
        return items
    return []


def resolve_anchor(pdb, het, token, cif_atoms, gpcrdb_atoms, product_instances, label):
    """
    Decide which product instance is the ligand copy GPCRdb calls `token`.

    Returns a dict with the ANCHOR_COLUMNS fields. Rules:

    1. exact coordinates: the GPCRdb residue (chain, number, name[:3]) is
       compared atom by atom with every product instance of the HET; exactly
       one instance with the most shared atoms decides (source coord_exact).
    2. the label route names an instance via label_asym_id; it confirms rule 1
       (source coord+label), or disagrees (status errata, rule 1 still wins),
       or is absent (source coord_only).
    3. no shared atoms at all (GPCRdb stores an older model): the name-based
       candidate HET_<chain>_<number> is accepted only when the label route
       names the same instance (source fallback+label); otherwise unresolved.
    4. the product has no instance of this HET: no_product.
    """
    het = het.upper()
    copies = sorted(i for i in product_instances if i.split("_", 1)[0].upper() == het)
    row = {
        "pdb": pdb,
        "het": het,
        "token": token,
        "instance": "",
        "status": "",
        "source": "",
        "n_exact": 0,
        "n_gpcrdb": 0,
        "label": label or "",
        "label_instance": "",
        "note": "",
    }
    m = TOKEN_RE.match(token)
    chain, resnum, icode = m.group("chain"), m.group("resnum"), m.group("icode")

    by_instance = collections.defaultdict(set)
    label_candidates = set()
    for a in cif_atoms:
        if a["comp"].upper() != het:
            continue
        name = instance_name(het, a["auth_asym"], a["auth_seq"], a["icode"])
        by_instance[name].add(a["key"])
        if label and a["label_asym"] == label:
            label_candidates.add(name)
    if len(label_candidates) == 1:
        row["label_instance"] = next(iter(label_candidates))
    elif len(label_candidates) > 1:
        row["note"] = "label names several residues"

    gkeys = {
        a["key"]
        for a in gpcrdb_atoms
        if a["chain"] == chain
        and a["resnum"] == resnum
        and a["icode"] == icode
        and a["resname"].upper() == het[:3]
    }
    row["n_gpcrdb"] = len(gkeys)

    if not copies:
        row["status"], row["note"] = (
            "no_product",
            (row["note"] or "product has no instance of this HET"),
        )
        return row

    scores = sorted(
        ((len(gkeys & by_instance.get(c, set())), c) for c in copies), reverse=True
    )
    best, best_name = scores[0]
    if best > 0:
        if len(scores) > 1 and scores[1][0] == best:
            row["status"], row["note"] = (
                "unresolved",
                "two instances share the same coordinates",
            )
            return row
        row["instance"], row["n_exact"] = best_name, best
        if not row["label_instance"]:
            row["status"], row["source"] = "ok", "coord_only"
        elif row["label_instance"] == best_name:
            row["status"], row["source"] = "ok", "coord+label"
        else:
            row["status"], row["source"] = "errata", "coord_exact"
            row["note"] = "annotation label_asym_id points at {}".format(
                row["label_instance"]
            )
        return row

    fallback = instance_name(het, chain, resnum, icode)
    if fallback in copies and row["label_instance"] == fallback:
        row["instance"], row["status"], row["source"] = fallback, "ok", "fallback+label"
        row["note"] = "no exact coordinates (older model in GPCRdb?)"
        return row
    row["status"] = "unresolved"
    row["note"] = (
        "no exact coordinates; name candidate {} {}; label candidate {}".format(
            fallback,
            "exists" if fallback in copies else "absent",
            row["label_instance"] or "-",
        )
    )
    return row


# ---------------------------------------------------------------------------
# Receptor map
# ---------------------------------------------------------------------------

RECEPTOR_COLUMNS = (
    "pdb",
    "preferred_chain",
    "auth_chain",
    "status",
    "method",
    "n_ca_gpcrdb",
    "n_ca_matched",
    "renumbered",
    "note",
    "gpcrdb_text_sha256",
    "product_instances_sha256",
)


def text_sha256(text):
    """sha256 of GPCRdb's stored structure text, as stored."""
    return hashlib.sha256((text or "").encode("utf-8")).hexdigest()


def instances_sha256(names):
    """sha256 of a structure's sorted product instance names."""
    return hashlib.sha256("\n".join(sorted(names)).encode("utf-8")).hexdigest()


def resolve_receptor(pdb, preferred_chain, cif_atoms, gpcrdb_atoms):
    """
    Name the product (author) chain whose CA atoms are GPCRdb's preferred chain.

    exact: the author chain sharing the most CA coordinates with the GPCRdb
    preferred chain; every matched CA must keep its residue number, otherwise
    status renumbered (the importer must refuse the structure).
    identity_drift: no CA matches at all (older model in GPCRdb); accepted only
    if an author chain of the same name carries every GPCRdb (number, name)
    pair of the preferred chain.
    """
    pref = (preferred_chain or "").split(",")[0].strip()
    row = {
        "pdb": pdb,
        "preferred_chain": pref,
        "auth_chain": "",
        "status": "",
        "method": "",
        "n_ca_gpcrdb": 0,
        "n_ca_matched": 0,
        "renumbered": 0,
        "note": "",
    }
    g_ca = {
        a["key"]: (a["resnum"], a["icode"], a["resname"])
        for a in gpcrdb_atoms
        if a["group"] == "ATOM" and a["atom"] == "CA" and a["chain"] == pref
    }
    row["n_ca_gpcrdb"] = len(g_ca)
    if not g_ca:
        row["status"], row["note"] = (
            "unresolved",
            "GPCRdb text has no CA on the preferred chain",
        )
        return row
    per_chain = collections.Counter()
    renumbered = collections.Counter()
    for a in cif_atoms:
        if a["group"] == "ATOM" and a["atom"] == "CA" and a["key"] in g_ca:
            per_chain[a["auth_asym"]] += 1
            if (a["auth_seq"], a["icode"]) != g_ca[a["key"]][:2]:
                renumbered[a["auth_asym"]] += 1
    if per_chain:
        ((chain, n),) = per_chain.most_common(1)
        row.update(
            auth_chain=chain,
            n_ca_matched=n,
            method="exact",
            renumbered=renumbered[chain],
        )
        row["status"] = "renumbered" if renumbered[chain] else "ok"
        if len(per_chain) > 1:
            row["note"] = "CA also matched on " + ",".join(
                "%s:%d" % kv for kv in sorted(per_chain.items()) if kv[0] != chain
            )
        return row
    wanted = {(v[0], v[1], v[2].upper()) for v in g_ca.values()}
    have = {
        (a["auth_seq"], a["icode"], a["comp"].upper())
        for a in cif_atoms
        if a["group"] == "ATOM" and a["atom"] == "CA" and a["auth_asym"] == pref
    }
    if have and wanted <= have:
        row.update(
            auth_chain=pref,
            method="identity_drift",
            status="ok",
            note="no CA coordinate matches; every GPCRdb (number, name) found on the same-named chain",
        )
    else:
        row["status"] = "unresolved"
        row["note"] = (
            "no CA coordinate matches and the same-named chain does not carry the GPCRdb residues"
        )
    return row
