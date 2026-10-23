r"""
Build one chainmap.tsv per PDB, from files only.

The per-anchor and per-receptor resolution is interaction.schrodinger_chain_map.
Every input comes from files, not the database, and the output is one file per
structure, written into that structure's product directory, where the importer
looks for it::

    python manage.py build_schrodinger_chainmap_files \
        --data-dir <DATA_DIR>/structure_data/schrodinger/engine1

build_all runs it before the imports, so the maps always match the annotation
and structure text of that build. The author side of the matching is the
coordinate index the producer delivers in each structure's product directory
(<PDB>/<PDB>_structure_index.tsv, the atoms of the mmCIF the products were
computed from); no mmCIF is read here.

Where each input comes from:

    database                              file
    ------------------------------------  ---------------------------------------
    StructureLigandInteraction            structure_data/annotation/ligands.tsv
      pdb_reference / chain_res             Name / Residue_seq_id
    Structure.preferred_chain             structure_data/annotation/structures.tsv
                                            ChainID (column 6)
    Structure.pdb_data.pdb                structure_data/pdbs/<PDB>.pdb
    structure_type.origin == experiment   structures.tsv holds experimental only
    author chains, numbering, coordinates <index-dir>/<PDB>/<PDB>_structure_index.tsv

The command issues no database query.

Writes <out-dir>/<PDB>/chainmap.tsv, one per PDB, each carrying its own
provenance header, so maps built at different times can sit in one tree.

An empty body means the annotation lists no ligand Engine 1 serves for that
structure. It does not mean the engine looked and found none: the body is built
from ligands.tsv, never from the product tree.
"""

import csv
import hashlib
import inspect
import io
import os
import subprocess

import yaml

from django.conf import settings
from django.core.management.base import BaseCommand, CommandError

from interaction import schrodinger_chain_map as cm
from interaction import schrodinger_import as si

# The format constants live in schrodinger_import, the module that reads these
# files back; one definition, so a bump cannot land on only one side.
SCHEMA = si.CHAINMAP_SCHEMA
RECEPTOR_PREFIX = si.CHAINMAP_RECEPTOR_PREFIX


def _sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def annotation_commit(given, gdata):
    """
    The gpcrdb_data commit to record: the one given, else the HEAD of ``gdata``.

    Only recorded in the map headers, never compared, so when git cannot answer
    for that very directory (not an enclosing repository) the value is
    "unknown" rather than a stopped build.
    """
    if given:
        return given

    def git(*args):
        return subprocess.run(
            ["git", "-C", gdata] + list(args),
            stdout=subprocess.PIPE,
            stderr=subprocess.DEVNULL,
            universal_newlines=True,
            check=True,
            timeout=30,
        ).stdout.strip()

    try:
        if os.path.realpath(git("rev-parse", "--show-toplevel")) != os.path.realpath(
            gdata
        ):
            return "unknown"
        return git("rev-parse", "--short=7", "HEAD") or "unknown"
    except (OSError, subprocess.SubprocessError):
        return "unknown"


def product_input_sha256(tree, pdb):
    """
    The input sha256 Engine 1 recorded in <tree>/<PDB>/summary.yaml.

    None when there is no summary (no recorded input sha256 to compare
    against); "" when a summary is there but unreadable or without the field,
    which the contract does not allow and index_mismatch refuses.
    """
    path = os.path.join(tree, pdb, si.PRODUCT_SUMMARY_NAME)
    if not os.path.exists(path):
        return None
    try:
        with open(path, encoding="utf-8") as fh:
            doc = yaml.safe_load(fh)
    except (OSError, UnicodeDecodeError, yaml.YAMLError):
        return ""
    value = doc.get("input_sha256") if isinstance(doc, dict) else None
    return value or ""


def index_mismatch(summary_sha, index_sha):
    """
    The reason to refuse the coordinate index, or None.

    Refused: a summary without input_sha256, and an index read from another
    mmCIF than the products.
    """
    if summary_sha == "":
        return (
            "summary.yaml is unreadable or has no input_sha256 to check the "
            "coordinate index against"
        )
    if summary_sha and summary_sha != index_sha:
        return (
            "the coordinate index was read from another mmCIF than the products "
            "(index {}, summary.yaml {})".format(index_sha[:12], summary_sha[:12])
        )
    return None


def _sha256_parts(parts):
    """
    sha256 over a list, length-prefixed so the parts cannot be re-cut.

    One of the parts is a function body and carries newlines, so joining on a
    separator would not be injective.
    """
    h = hashlib.sha256()
    for part in parts:
        raw = part.encode("utf-8")
        h.update(str(len(raw)).encode("ascii") + b":" + raw)
    return h.hexdigest()


# The header keys that are the same for every structure of a run, plus the two
# that are not, in the order they are written.
HEADER_KEYS = (
    "annotation_commit",
    "ligands_sha256",
    "structures_sha256",
    "gpcrdb_pdb_sha256",
    "cif_sha256",
    "builder_sha256",
)


def has_product_summary(data_dir, pdb):
    """
    Did the producer leave its per-structure summary in this directory?

    Read from the product tree, never from the out dir, which may be elsewhere.
    """
    path = os.path.join(data_dir, pdb, si.PRODUCT_SUMMARY_NAME)
    try:
        # An empty file is a truncated copy, not a run.
        return os.path.isfile(path) and os.path.getsize(path) > 0
    except OSError:
        return False


def chainmap_header(pdb, values, receptor, has_summary):
    """
    The '# key<TAB>value' lines of one chainmap, in order.

    Every field is about this structure or about an input common to all of
    them, so each file stays true when trees are merged.
    """
    header = [("schema", SCHEMA), ("pdb", pdb)]
    header += [(key, values[key]) for key in HEADER_KEYS]
    header.append((si.PRODUCT_SUMMARY_KEY, "yes" if has_summary else "no"))
    header += [
        (RECEPTOR_PREFIX + c, receptor.get(c, ""))
        for c in cm.RECEPTOR_COLUMNS
        if c != "pdb"
    ]
    return header


def _some(names, limit=20):
    """A list for a log line, honest about what it left out."""
    head = ", ".join(names[:limit])
    return (
        head
        if len(names) <= limit
        else "{} (+{} more)".format(head, len(names) - limit)
    )


def _flat(value):
    """
    A header value may hold free text (note); keep it on one line, one field.

    Only the characters that would break the format are replaced, so the value
    still reads as what it was.
    """
    text = str(value if value is not None else "")
    for bad in ("\t", "\r", "\n"):
        text = text.replace(bad, " ")
    return text


def read_tsv(path):
    """
    Rows of a gpcrdb_data annotation TSV, keys and values stripped.

    A row wider than the header arrives under the key None as a list: refuse
    the file rather than guess which field went where. A row that is too short
    is padded with empty strings, which drops the anchor from the map and makes
    the importer refuse that structure -- loud, and the safe direction.
    """
    rows = []
    with open(path, newline="") as fh:
        for n, r in enumerate(csv.DictReader(fh, delimiter="\t"), start=2):
            if None in r:
                raise CommandError(
                    "{} line {}: more fields than the header has columns".format(
                        path, n
                    )
                )
            rows.append({(k or "").strip(): (v or "").strip() for k, v in r.items()})
    return rows


def preferred_chains(rows):
    """
    PDB -> preferred chain, as structure.functions.ParseStructureCSV stores it.

    That parser keeps only the first character of a chain id containing a dot.
    The dot rule is applied here and the comma rule by resolve_receptor, as on
    the database side.
    """
    out = {}
    for r in rows:
        pdb = (r.get("PDB") or "").upper()
        chain = r.get("ChainID") or ""
        if "." in chain:
            chain = chain[0]
        if pdb:
            out[pdb] = chain
    return out


def annotation_anchors(rows):
    """
    PDB -> [(HET, token, chain_res)] in file order, for the anchors Engine 1 serves.

    Mirrors schrodinger_import.is_in_scope: every row whose Name is a real
    chemical component, whatever its Type (the type can disagree between the
    annotation and the database). The (HET, token) pairs are the set
    check_map_covers compares against the database.
    """
    out = {}
    seen = set()
    for r in rows:
        pdb = (r.get("PDB") or "").upper()
        het = (r.get("Name") or "").strip().upper()
        if not pdb or not het or het in si.PLACEHOLDER_REFERENCES:
            continue
        chain_res = r.get("Residue_seq_id") or ""
        for tok in cm.split_tokens(chain_res) or [""]:
            if (pdb, het, tok) in seen:
                continue
            seen.add((pdb, het, tok))
            out.setdefault(pdb, []).append((het, tok, chain_res))
    return out


def write_chainmap(path, header, anchor_rows):
    buf = io.StringIO()
    for key, value in header:
        buf.write("# {}\t{}\n".format(key, _flat(value)))
    w = csv.DictWriter(
        buf, fieldnames=cm.ANCHOR_COLUMNS, delimiter="\t", lineterminator="\n"
    )
    w.writeheader()
    w.writerows(anchor_rows)
    # The out dir is usually the product tree the importer reads; never leave a
    # half written map behind if the run dies on the next structure.
    tmp = path + ".tmp"
    with open(tmp, "w", newline="") as fh:
        fh.write(buf.getvalue())
    os.replace(tmp, path)


class Command(BaseCommand):
    help = (
        "Build one per-PDB chainmap.tsv for the Schrodinger importer, from files only."
    )

    def add_arguments(self, parser):
        parser.add_argument(
            "--gpcrdb-data",
            default=None,
            help="gpcrdb_data checkout (uses structure_data/annotation and "
            "structure_data/pdbs). Defaults to DATA_DIR, so the maps are "
            "built from the same annotation the build will run on.",
        )
        parser.add_argument(
            "--data-dir",
            required=True,
            help="Product tree {data_dir}/{PDB}/{instance}/.",
        )
        parser.add_argument(
            "--index-dir",
            default=None,
            help="Where the coordinate indexes are, {index_dir}/{PDB}/{PDB}"
            + cm.INDEX_SUFFIX
            + "; defaults to --data-dir, where the "
            "producer delivers them.",
        )
        parser.add_argument(
            "--annotation-commit",
            default=None,
            help="gpcrdb_data commit the annotation was read at; recorded only. "
            "Defaults to the HEAD of --gpcrdb-data when git can read it "
            "there, 'unknown' otherwise.",
        )
        parser.add_argument(
            "--out-dir",
            default=None,
            help="Where to write {PDB}/chainmap.tsv; defaults to --data-dir, "
            "where the importer looks for it.",
        )
        parser.add_argument(
            "--allow-stray",
            action="store_true",
            help="Do not refuse product directories that are absent from "
            "structures.tsv. They still get no chainmap, and the "
            "importer fails on any of them the database still knows.",
        )
        parser.add_argument(
            "--pdb",
            action="append",
            default=[],
            help="Restrict to these PDB codes; default is every structure in "
            "structures.tsv.",
        )

    def handle(self, *args, **opt):
        gdata = opt["gpcrdb_data"] or settings.DATA_DIR
        data_dir = opt["data_dir"]
        index_dir = opt["index_dir"] or data_dir
        out_dir = opt["out_dir"] or data_dir
        commit = annotation_commit(opt["annotation_commit"], gdata)
        for label, path in (
            ("--gpcrdb-data", gdata),
            ("--index-dir", index_dir),
            ("--data-dir", data_dir),
        ):
            if not os.path.isdir(path):
                raise CommandError("{} {!r} is not a directory".format(label, path))
        # The out dir itself is created, its parent is not: a typo fails.
        out_parent = os.path.dirname(os.path.abspath(out_dir))
        if not os.path.isdir(out_parent):
            raise CommandError(
                "--out-dir {!r}: {} does not exist".format(out_dir, out_parent)
            )
        ann = os.path.join(gdata, "structure_data", "annotation")
        ligands_tsv = os.path.join(ann, "ligands.tsv")
        structures_tsv = os.path.join(ann, "structures.tsv")
        pdb_dir = os.path.join(gdata, "structure_data", "pdbs")
        for path in (ligands_tsv, structures_tsv):
            if not os.path.isfile(path):
                raise CommandError("{} is missing under --gpcrdb-data".format(path))
        if not os.path.isdir(pdb_dir):
            raise CommandError("{} is missing under --gpcrdb-data".format(pdb_dir))

        ligand_rows = read_tsv(ligands_tsv)
        labels = cm.annotation_labels(ligand_rows)
        anchors = annotation_anchors(ligand_rows)
        chains = preferred_chains(read_tsv(structures_tsv))

        # The corpus is structures.tsv, not the directories under --data-dir:
        # a structure whose run yielded nothing still needs its map, and a
        # stray directory must not get one.
        corpus = sorted(chains)
        if not corpus:
            raise CommandError("{} lists no structures".format(structures_tsv))
        if opt["pdb"]:
            pdbs = sorted({p.upper() for p in opt["pdb"]})
            unknown = [p for p in pdbs if p not in chains]
            if unknown:
                raise CommandError(
                    "not in {}: {}".format(structures_tsv, ", ".join(unknown))
                )
        else:
            pdbs = corpus
            # A directory holding products but absent from the corpus gets no
            # chainmap.tsv, and the importer fails on it if the database knows
            # the structure. Directories without an instance (.git, logs) are
            # only noted.
            dropped, housekeeping = [], []
            for name in sorted(os.listdir(data_dir)):
                if name.upper() in chains or not os.path.isdir(
                    os.path.join(data_dir, name)
                ):
                    continue
                (
                    dropped if si.instance_yaml_paths(data_dir, name) else housekeeping
                ).append(name)
            if housekeeping:
                self.stdout.write(
                    "directories that are not structures, ignored: {}".format(
                        _some(housekeeping)
                    )
                )
            if dropped and not opt["allow_stray"]:
                raise CommandError(
                    "these directories hold Engine 1 products but are not in {}, so they get "
                    "no chainmap: {}. Pass --allow-stray if that is intended, for "
                    "instance after a structure was retired from the annotation.".format(
                        structures_tsv, _some(dropped)
                    )
                )
            if dropped:
                self.stdout.write(
                    "product directories with no chainmap (--allow-stray): {}".format(
                        _some(dropped)
                    )
                )

        ligands_sha = _sha256_file(ligands_tsv)
        structures_sha = _sha256_file(structures_tsv)
        # What produced these rows: the algorithm module, this command, and what
        # it borrows from the importer module (scope, instance discovery and the
        # chainmap format constants), so an unrelated importer edit does not move
        # the stamp.
        builder_sha = _sha256_parts(
            [
                _sha256_file(cm.__file__),
                _sha256_file(os.path.abspath(__file__)),
                repr(sorted(si.PLACEHOLDER_REFERENCES)),
                si.INSTANCE_DIR_RE.pattern,
                repr(si.INSTANCE_DIR_RE.flags),
                inspect.getsource(si.instance_yaml_paths),
                si.CHAINMAP_SCHEMA,
                si.CHAINMAP_RECEPTOR_PREFIX,
                si.PRODUCT_SUMMARY_NAME,
                si.PRODUCT_SUMMARY_KEY,
            ]
        )

        counts, rstatus, written, no_products, checked, refused = {}, {}, 0, 0, 0, 0
        for pdb in pdbs:
            has_summary = has_product_summary(data_dir, pdb)
            summary_sha = product_input_sha256(data_dir, pdb)
            checked += bool(summary_sha)
            refused += summary_sha == ""
            rows, receptor, note = self.build_one(
                pdb,
                anchors.get(pdb, []),
                chains.get(pdb),
                labels,
                cm.index_path(index_dir, pdb),
                os.path.join(pdb_dir, pdb + ".pdb"),
                si.instance_yaml_paths(data_dir, pdb),
                has_summary,
                summary_sha,
            )
            no_products += receptor["product_instances_sha256"] == cm.instances_sha256(
                []
            )
            header = chainmap_header(
                pdb,
                dict(
                    annotation_commit=commit,
                    ligands_sha256=ligands_sha,
                    structures_sha256=structures_sha,
                    builder_sha256=builder_sha,
                    # The bytes on disk. receptor.gpcrdb_text_sha256 is what the
                    # importer compares against the database, hashed as decoded
                    # text; the two differ only on a CRLF file.
                    gpcrdb_pdb_sha256=note["gpcrdb_pdb_sha256"],
                    cif_sha256=note["cif_sha256"],
                ),
                receptor,
                has_summary,
            )
            target = os.path.join(out_dir, pdb)
            os.makedirs(target, exist_ok=True)
            write_chainmap(os.path.join(target, "chainmap.tsv"), header, rows)
            written += 1
            for r in rows:
                counts[(r["status"], r["source"])] = (
                    counts.get((r["status"], r["source"]), 0) + 1
                )
            key = (receptor["status"], receptor["method"])
            rstatus[key] = rstatus.get(key, 0) + 1

        self.stdout.write("out-dir {}".format(out_dir))
        self.stdout.write(
            "annotation_commit {} ligands_sha256 {} structures_sha256 {} "
            "builder_sha256 {}".format(commit, ligands_sha, structures_sha, builder_sha)
        )
        self.stdout.write(
            "chainmap.tsv written: {} ({} with no product instance)".format(
                written, no_products
            )
        )
        self.stdout.write(
            "summary.yaml with input_sha256: {} "
            "(no summary: {}; summary unreadable or without it, refused: {})".format(
                checked, written - checked - refused, refused
            )
        )
        self.stdout.write(
            "anchor rows {}: {}".format(sum(counts.values()), sorted(counts.items()))
        )
        self.stdout.write(
            "receptor rows {}: {}".format(
                sum(rstatus.values()), sorted(rstatus.items())
            )
        )

    def build_one(
        self,
        pdb,
        anchor_keys,
        preferred_chain,
        labels,
        index_path,
        gpcrdb_pdb_path,
        instances,
        has_summary,
        summary_sha=None,
    ):
        """
        (anchor_rows, receptor_row, provenance) for one structure.

        An unreadable input is not a reason to skip: every anchor is written as
        unresolved with the reason, and the importer refuses the structure. A
        silently missing structure would instead look like one with no anchors.
        """
        note = {
            "cif_sha256": "",
            "gpcrdb_pdb_sha256": "",
            "product_summary": has_summary,
        }
        if preferred_chain is None:
            preferred_chain = ""
        gtext = None

        def unresolved(exc):
            reason = (
                exc
                if isinstance(exc, str)
                else "input unreadable: {}: {}".format(type(exc).__name__, exc)
            )[:200]
            rows = [
                dict(
                    {c: "" for c in cm.ANCHOR_COLUMNS},
                    pdb=pdb,
                    het=het,
                    token=tok,
                    status="unresolved",
                    note=reason,
                )
                for het, tok, _ in anchor_keys
            ]
            receptor = dict(
                {c: "" for c in cm.RECEPTOR_COLUMNS},
                pdb=pdb,
                preferred_chain=preferred_chain,
                status="unresolved",
                note=reason,
                # Whenever the stored text was read, so the importer can tell
                # a map built from another dump.
                gpcrdb_text_sha256="" if gtext is None else cm.text_sha256(gtext),
                product_instances_sha256=cm.instances_sha256(instances),
            )
            return rows, receptor, note

        try:
            with open(index_path) as fh:
                index_text = fh.read()
            note["gpcrdb_pdb_sha256"] = _sha256_file(gpcrdb_pdb_path)
            with open(gpcrdb_pdb_path) as fh:
                gtext = fh.read()
        except (OSError, UnicodeDecodeError) as exc:
            return unresolved(exc)
        try:
            # Only the parsers get the wide clause: one malformed input costs
            # one structure, not the run; anything raised elsewhere is a bug.
            # The header records the sha256 of the mmCIF the index was read
            # from, the file the products were computed from.
            note["cif_sha256"], cif_atoms = cm.parse_structure_index(index_text)
            gatoms = cm.parse_gpcrdb_pdb(gtext)
        except (cm.ParseError, KeyError, ValueError) as exc:
            return unresolved(exc)
        mismatch = index_mismatch(summary_sha, note["cif_sha256"])
        if mismatch:
            return unresolved(mismatch)

        receptor = cm.resolve_receptor(pdb, preferred_chain, cif_atoms, gatoms)
        receptor["gpcrdb_text_sha256"] = cm.text_sha256(gtext)
        receptor["product_instances_sha256"] = cm.instances_sha256(instances)

        rows = []
        for het, tok, chain_res in anchor_keys:
            if tok:
                rows.append(
                    cm.resolve_anchor(
                        pdb,
                        het,
                        tok,
                        cif_atoms,
                        gatoms,
                        instances,
                        labels.get((pdb, het, tok)),
                    )
                )
                continue
            copies = sorted(n for n in instances if n.split("_", 1)[0].upper() == het)
            rows.append(
                dict(
                    {c: "" for c in cm.ANCHOR_COLUMNS},
                    pdb=pdb,
                    het=het,
                    token="",
                    instance=";".join(copies),
                    status="all_copies" if copies else "no_product",
                    note=(
                        "chain_res {!r} names no residue; every copy used".format(
                            chain_res
                        )
                        if copies
                        else "the product has no instance of {} (chain_res {!r})".format(
                            het, chain_res
                        )
                    ),
                )
            )
        return rows, receptor, note
