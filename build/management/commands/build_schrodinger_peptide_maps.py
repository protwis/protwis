"""Build one peptide_map.tsv per structure of an Engine 2 tree, from files only.

    python manage.py build_schrodinger_peptide_maps \\
        --data-dir <DATA_DIR>/structure_data/schrodinger/engine2 \\
        --index-dir <DATA_DIR>/structure_data/schrodinger/engine1

build_all runs it before the imports, so the maps always match the annotation
and structure text of that build. The author side comes from the coordinate
index the producer delivers with the Engine 1 products; no mmCIF is read here.
This relies on Engine 2 computing on Engine 1's prepared structure of the same
mmCIF: the Engine 2 products carry no record of their own input, so the index
is checked against Engine 1's summary.yaml only.

Where each input comes from:

    decision                              file
    ------------------------------------  ---------------------------------------
    which chains carry a "pep" ligand     structure_data/annotation/ligands.tsv
                                            Name / ChainID (Type and Title are
                                            only recorded)
    GPCRdb's receptor chain               structure_data/annotation/structures.tsv
                                            ChainID
    GPCRdb's coordinates                  structure_data/pdbs/<PDB>.pdb
    author chains and numbering           --index-dir/<PDB>/<PDB>_structure_index.tsv
    segments, work items, outcomes        --data-dir/<PDB>/plan.json and the
                                            item records beside it

The map lists every chain the annotation names for a "pep" ligand of the
structure, whatever its type, not only the one the database anchor chose; the
importer looks up which one that is. The annotation's Type column is recorded
only, never used for scope. The rules are in interaction.schrodinger_peptide.

The command issues no database query.
"""

import inspect
import os

from django.conf import settings
from django.core.management.base import BaseCommand, CommandError

from build.management.commands import build_schrodinger_chainmap_files as e1
from interaction import schrodinger_chain_map as cm
from interaction import schrodinger_peptide as sp


def peptide_chains(rows):
    """PDB -> {gpcrdb chain: (titles, types)} for every "pep" ligand of the annotation.

    No type filter: the importer serves every "pep" anchor whatever its type,
    and a map that left a chain out would make the importer refuse the
    structure.
    """
    out = {}
    for r in rows:
        pdb = (r.get("PDB") or "").upper()
        if not pdb or (r.get("Name") or "").strip().upper() != sp.PEPTIDE_REFERENCE:
            continue
        ltype = r.get("Type") or ""
        for chain in [
            c.strip() for c in (r.get("ChainID") or "").split(",") if c.strip()
        ]:
            titles, types = out.setdefault(pdb, {}).setdefault(chain, (set(), set()))
            titles.add(_flat(r.get("Title") or ""))
            types.add(ltype)
    return out


def _flat(value):
    text = str(value if value is not None else "")
    for bad in ("\t", "\r", "\n", sp.LIST_SEP):
        text = text.replace(bad, " ")
    return text


HEADER_KEYS = (
    "annotation_commit",
    "ligands_sha256",
    "structures_sha256",
    "gpcrdb_pdb_sha256",
    "cif_sha256",
    "plan_sha256",
    "builder_sha256",
)


class Command(BaseCommand):
    help = (
        "Build one peptide_map.tsv per structure of an Engine 2 tree, from files only."
    )

    def add_arguments(self, parser):
        parser.add_argument(
            "--gpcrdb-data",
            default=None,
            help="gpcrdb_data checkout (structure_data/annotation and "
            "structure_data/pdbs). Defaults to DATA_DIR.",
        )
        parser.add_argument(
            "--data-dir",
            required=True,
            help="Engine 2 tree {data_dir}/{PDB}/plan.json.",
        )
        parser.add_argument(
            "--index-dir",
            required=True,
            help="Where the coordinate indexes are, {index_dir}/{PDB}/{PDB}"
            + cm.INDEX_SUFFIX
            + ": the Engine 1 tree, where the producer "
            "delivers them.",
        )
        parser.add_argument(
            "--annotation-commit",
            default=None,
            help="gpcrdb_data commit the annotation was read at; recorded only. "
            "Defaults to the HEAD of --gpcrdb-data when git can read it "
            "there, 'unknown' otherwise.",
        )
        parser.add_argument(
            "--pdb", action="append", default=[], help="Restrict to these PDB codes."
        )

    def handle(self, *args, **opt):
        gdata = opt["gpcrdb_data"] or settings.DATA_DIR
        index_dir, data_dir = opt["index_dir"], opt["data_dir"]
        commit = e1.annotation_commit(opt["annotation_commit"], gdata)
        for label, path in (
            ("--gpcrdb-data", gdata),
            ("--index-dir", index_dir),
            ("--data-dir", data_dir),
        ):
            if not os.path.isdir(path):
                raise CommandError("{} {!r} is not a directory".format(label, path))
        ann = os.path.join(gdata, "structure_data", "annotation")
        ligands_tsv = os.path.join(ann, "ligands.tsv")
        structures_tsv = os.path.join(ann, "structures.tsv")
        pdb_dir = os.path.join(gdata, "structure_data", "pdbs")
        for path in (ligands_tsv, structures_tsv):
            if not os.path.isfile(path):
                raise CommandError("{} is missing under --gpcrdb-data".format(path))

        peptides = peptide_chains(e1.read_tsv(ligands_tsv))
        chains = e1.preferred_chains(e1.read_tsv(structures_tsv))
        # The corpus: structures the annotation lists and the producer planned.
        planned = {
            name.upper(): name
            for name in os.listdir(data_dir)
            if os.path.isfile(os.path.join(data_dir, name, sp.PLAN_NAME))
        }
        if opt["pdb"]:
            pdbs = sorted({p.upper() for p in opt["pdb"]})
            unknown = [p for p in pdbs if p not in planned or p not in chains]
            if unknown:
                raise CommandError(
                    "not planned or not in structures.tsv: {}".format(
                        ", ".join(unknown)
                    )
                )
        else:
            pdbs = sorted(p for p in planned if p in chains)
            missing = sorted(p for p in peptides if p in chains and p not in planned)
            if missing:
                self.stdout.write(
                    "structures with peptides and no Engine 2 plan "
                    "(the importer fails them if the database has the anchors): "
                    "{}".format(e1._some(missing))
                )

        builder_sha = e1._sha256_parts(
            [
                sp.sha256_file(sp.__file__),
                sp.sha256_file(cm.__file__),
                sp.sha256_file(os.path.abspath(__file__)),
                inspect.getsource(e1.read_tsv),
                inspect.getsource(e1.preferred_chains),
            ]
        )
        common = dict(
            annotation_commit=commit,
            ligands_sha256=sp.sha256_file(ligands_tsv),
            structures_sha256=sp.sha256_file(structures_tsv),
            builder_sha256=builder_sha,
        )

        counts, rstatus, written, checked, refused = {}, {}, 0, 0, 0
        for pdb in pdbs:
            name = planned[pdb]
            summary_sha = e1.product_input_sha256(index_dir, pdb)
            checked += bool(summary_sha)
            refused += summary_sha == ""
            receptor, rows, values = self.build_one(
                pdb,
                os.path.join(data_dir, name),
                chains.get(pdb) or "",
                peptides.get(pdb, {}),
                cm.index_path(index_dir, pdb),
                os.path.join(pdb_dir, pdb + ".pdb"),
                summary_sha,
            )
            header = [("pdb", pdb)] + [
                (k, dict(common, **values).get(k, "")) for k in HEADER_KEYS
            ]
            sp.write_peptide_map(
                os.path.join(data_dir, name, sp.MAP_NAME),
                header,
                # Header values may hold commas (segments is a list);
                # only what breaks a header line is replaced.
                {k: e1._flat(v) for k, v in receptor.items()},
                rows,
            )
            written += 1
            for r in rows:
                counts[r["status"]] = counts.get(r["status"], 0) + 1
            rstatus[receptor["status"]] = rstatus.get(receptor["status"], 0) + 1

        self.stdout.write(
            "annotation_commit {} builder_sha256 {}".format(commit, builder_sha)
        )
        self.stdout.write("peptide_map.tsv written: {}".format(written))
        self.stdout.write(
            "summary.yaml with input_sha256: {} "
            "(no summary: {}; summary unreadable or without it, refused: {})".format(
                checked, written - checked - refused, refused
            )
        )
        self.stdout.write(
            "peptide chain rows {}: {}".format(
                sum(counts.values()), sorted(counts.items())
            )
        )
        self.stdout.write("receptor {}".format(sorted(rstatus.items())))

    def build_one(
        self,
        pdb,
        tree,
        preferred_chain,
        chain_info,
        index_path,
        gpcrdb_pdb_path,
        summary_sha=None,
    ):
        """(receptor, rows, header values) for one structure.

        An unreadable input is not a reason to skip: the receptor is written
        unresolved with the reason, and the importer refuses the structure.
        """
        values = {"plan_sha256": "", "cif_sha256": "", "gpcrdb_pdb_sha256": ""}
        receptor = {k: "" for k in sp.RECEPTOR_KEYS}
        receptor["preferred_chain"] = preferred_chain.split(",")[0].strip()
        gtext = None

        def unresolved(reason, rows_status="chain_unresolved"):
            receptor.update(status="unresolved", note=reason[:200])
            # Whenever the stored text was read, so the importer can tell a
            # map built from another dump.
            if gtext is not None:
                receptor["gpcrdb_text_sha256"] = cm.text_sha256(gtext)
            rows = [
                self._row(pdb, chain, info, status=rows_status, note=reason[:200])
                for chain, info in sorted(chain_info.items())
            ]
            return receptor, rows, values

        try:
            values["plan_sha256"] = sp.sha256_file(os.path.join(tree, sp.PLAN_NAME))
            segments, items = sp.load_plan(
                os.path.dirname(tree), os.path.basename(tree)
            )
            # The header records the sha256 of the mmCIF the index was read from.
            with open(index_path) as fh:
                values["cif_sha256"], cif_atoms = cm.parse_structure_index(fh.read())
            mismatch = e1.index_mismatch(summary_sha, values["cif_sha256"])
            values["gpcrdb_pdb_sha256"] = sp.sha256_file(gpcrdb_pdb_path)
            with open(gpcrdb_pdb_path) as fh:
                gtext = fh.read()
            gatoms = cm.parse_gpcrdb_pdb(gtext)
        except (
            OSError,
            UnicodeDecodeError,
            sp.MalformedProduct,
            cm.ParseError,
            KeyError,
            ValueError,
        ) as exc:
            return unresolved(
                "input unreadable: {}: {}".format(type(exc).__name__, exc)
            )
        if mismatch:
            return unresolved(mismatch)

        res = cm.resolve_receptor(pdb, receptor["preferred_chain"], cif_atoms, gatoms)
        receptor.update(
            auth_chain=res["auth_chain"],
            status=res["status"],
            method=res["method"],
            n_ca_gpcrdb=res["n_ca_gpcrdb"],
            n_ca_matched=res["n_ca_matched"],
            note=res["note"],
            gpcrdb_text_sha256=cm.text_sha256(gtext),
        )
        rseg, rsegs = None, []
        if res["status"] == "ok":
            numbers = sp.receptor_ca_numbers(
                res["auth_chain"], receptor["preferred_chain"], cif_atoms, gatoms
            )
            rseg, covered, why = sp.receptor_segment(
                segments, res["auth_chain"], numbers
            )
            receptor.update(segment=rseg or "", n_covered=covered)
            if rseg is None:
                receptor.update(status="no_segment", note=why)
            else:
                rsegs = sp.receptor_chain_segments(segments, res["auth_chain"])
                receptor.update(segments=sp.LIST_SEP.join(rsegs))
        rows = []
        for chain, info in sorted(chain_info.items()):
            if rseg is None:
                rows.append(
                    self._row(
                        pdb,
                        chain,
                        info,
                        status="no_receptor_segment",
                        note=receptor["note"],
                    )
                )
                continue
            auth, method, note = sp.peptide_author_chain(pdb, chain, cif_atoms, gatoms)
            if not auth:
                rows.append(
                    self._row(pdb, chain, info, status="chain_unresolved", note=note)
                )
                continue
            found = sp.peptide_items(segments, items, auth, rsegs)
            if not found:
                rows.append(
                    self._row(
                        pdb,
                        chain,
                        info,
                        auth_chain=auth,
                        chain_method=method,
                        status="no_items",
                        note="no peptide-as-ligand item against chain {} "
                        "(peptide on the receptor chain?)".format(
                            receptor["auth_chain"]
                        ),
                    )
                )
                continue
            outcomes = []
            for _, key in found:
                record = sp.read_json(
                    sp.item_paths(os.path.dirname(tree), os.path.basename(tree), key)[0]
                )
                outcomes.append(str(record.get("outcome") or ""))
            rows.append(
                self._row(
                    pdb,
                    chain,
                    info,
                    auth_chain=auth,
                    chain_method=method,
                    status=sp.ROW_OK,
                    note=note,
                    segments=sp.LIST_SEP.join(sorted({s for s, _ in found})),
                    items=sp.LIST_SEP.join(k for _, k in found),
                    outcomes=sp.LIST_SEP.join(outcomes),
                )
            )
        return receptor, rows, values

    @staticmethod
    def _row(pdb, chain, info, **fields):
        titles, types = info
        row = {c: "" for c in sp.MAP_COLUMNS}
        row.update(
            pdb=pdb,
            gpcrdb_chain=chain,
            titles=" | ".join(sorted(titles)),
            types=" | ".join(sorted(types)),
        )
        row.update({k: _flat(v) if k == "note" else v for k, v in fields.items()})
        return row
