r"""
End-to-end tests of the map builders and the clean-up on a small synthetic delivery.

No database: each command is called with the options build/ligand_imports plans.

    python -c "import django; django.setup(); import unittest; \
        unittest.main(module='interaction.test_schrodinger_map_builders', argv=['x'])"
"""

import hashlib
import io
import json
import os
import shutil
import tempfile
import unittest
from unittest import mock

from django.core.management import call_command
from django.test import override_settings

from build import ligand_imports
from build.management.commands import build_schrodinger_chainmap_files as e1
from interaction import schrodinger_chain_map as cm
from interaction import schrodinger_import as si
from interaction import schrodinger_peptide as sp

PDB = "1ABC"
EMPTY = "9XYZ"  # in the annotation, nothing delivered
SHA = hashlib.sha256(b"the mmCIF the products were computed from").hexdigest()
REC_SEG, PEP_SEG = "R_REC_1_3", "P_1ABC_1_2"
K_PEP = "1ABC__P_1ABC_1_2__R_REC_1_3__aaaaaaaaaaaa"
K_REC = "1ABC__R_REC_1_3__P_1ABC_1_2__bbbbbbbbbbbb"

# (GPCRdb chain, author chain, label, group, resname, resnum, atom, x, y, z)
ATOMS = [
    ("A", "R", "A", "ATOM", "ALA", 1, "CA", 1.0, 0.0, 0.0),
    ("A", "R", "A", "ATOM", "ALA", 2, "CA", 2.0, 0.0, 0.0),
    ("A", "R", "A", "ATOM", "ALA", 3, "CA", 3.0, 0.0, 0.0),
    ("A", "R", "B", "HETATM", "LIG", 401, "C1", 5.0, 5.0, 5.0),
    ("A", "R", "B", "HETATM", "LIG", 401, "C2", 6.0, 5.0, 5.0),
    ("P", "P", "C", "ATOM", "GLY", 1, "CA", 10.0, 0.0, 0.0),
    ("P", "P", "C", "ATOM", "GLY", 2, "CA", 11.0, 0.0, 0.0),
]


def gpcrdb_text():
    out = []
    for n, (chain, _a, _l, group, resname, resnum, atom, x, y, z) in enumerate(
        ATOMS, start=1
    ):
        out.append(
            "%-6s%5d %-4s %3s %1s%4s    %8.3f%8.3f%8.3f  1.00 20.00          %2s\n"
            % (group, n, atom, resname, chain, resnum, x, y, z, "C")
        )
    return "".join(out) + "END\n"


def index_text(sha=SHA, coord="%.3f"):
    rows = [
        "\t".join(
            [
                label,
                auth,
                resname,
                str(resnum),
                "",
                atom,
                group,
                coord % x,
                coord % y,
                coord % z,
            ]
        )
        for _c, auth, label, group, resname, resnum, atom, x, y, z in ATOMS
    ]
    head = [
        "# schema\t" + cm.INDEX_SCHEMA,
        "# source\t%s.cif" % PDB,
        "# cif_sha256\t" + sha,
        "# atoms\t%d" % len(rows),
        "\t".join(cm.INDEX_COLUMNS),
    ]
    return "\n".join(head + rows) + "\n"


def write(path, text):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as fh:
        fh.write(text)


class Delivery(object):
    """A gpcrdb_data checkout with one structure delivered by both engines."""

    def __init__(self, root, summary_sha=SHA, index=True, coord="%.3f"):
        """Write the annotation, the structure texts and both deliveries under root."""
        self.gdata = os.path.join(root, "gpcrdb_data")
        self.engine1 = os.path.join(
            self.gdata, "structure_data", "schrodinger", "engine1"
        )
        self.engine2 = os.path.join(
            self.gdata, "structure_data", "schrodinger", "engine2"
        )
        ann = os.path.join(self.gdata, "structure_data", "annotation")
        write(
            os.path.join(ann, "ligands.tsv"),
            "PDB\tChainID\tName\tType\tTitle\tlabel_asym_id\tResidue_seq_id\n"
            "%s\tA\tLIG\tsmall-molecule\tA ligand\tB\tA:401\n"
            "%s\tP\tpep\tpeptide\tA peptide\t\t\n" % (PDB, PDB),
        )
        write(
            os.path.join(ann, "structures.tsv"),
            "PDB\tChainID\n%s\tA\n%s\tA\n" % (PDB, EMPTY),
        )
        write(
            os.path.join(self.gdata, "structure_data", "pdbs", PDB + ".pdb"),
            gpcrdb_text(),
        )
        write(
            os.path.join(self.gdata, "structure_data", "pdbs", EMPTY + ".pdb"),
            gpcrdb_text(),
        )
        summary = (
            "contract_version: %s\n" % si.PRODUCT_CONTRACT
            + "pdb_id: %s\n" % PDB
            + ("input_sha256: %s\n" % summary_sha if summary_sha else "")
        )
        write(
            os.path.join(self.engine1, PDB, "summary.yaml"),
            summary + "ligand_interaction_summary: []\n",
        )
        write(
            os.path.join(self.engine1, PDB, "LIG_R_401", "LIG_R_401.yaml"),
            "result:\n  interactions: []\n",
        )
        if index:
            write(cm.index_path(self.engine1, PDB), index_text(coord=coord))
        plan = {
            "contract_version": sp.PRODUCT_CONTRACT,
            "segments": [
                {"name": REC_SEG, "chain_id": "R", "ranges": [[1, 3]]},
                {"name": PEP_SEG, "chain_id": "P", "ranges": [[1, 2]]},
            ],
            "items": [
                {"key": K_PEP, "ligand_segment": PEP_SEG, "receptor_segment": REC_SEG},
                {"key": K_REC, "ligand_segment": REC_SEG, "receptor_segment": PEP_SEG},
            ],
        }
        write(os.path.join(self.engine2, PDB, sp.PLAN_NAME), json.dumps(plan))
        for key in (K_PEP, K_REC):
            write(
                os.path.join(self.engine2, PDB, key, key + ".json"),
                json.dumps({"work_item_key": key, "outcome": "done"}),
            )
            write(
                os.path.join(self.engine2, PDB, key, key + ".yaml"),
                "result:\n  interactions: []\n",
            )

    def options(self):
        return {
            "skip_ligand_import": False,
            "engine1_data_dir": self.engine1,
            "engine2_data_dir": self.engine2,
            "engine1_report_dir": "/r1",
            "engine2_report_dir": "/r2",
        }

    def run(self, *commands):
        """Run the planned steps named in ``commands``, as the build plans them."""
        with override_settings(DATA_DIR=self.gdata):
            for command, kwargs in ligand_imports.steps(self.options()):
                if command in commands:
                    call_command(command, stdout=io.StringIO(), **kwargs)


def read_map(path):
    header, rows, cols = {}, [], None
    with open(path) as fh:
        for line in fh.read().splitlines():
            if line.startswith("# "):
                key, _, value = line[2:].partition("\t")
                header[key] = value
            elif cols is None:
                cols = line.split("\t")
            else:
                rows.append(dict(zip(cols, line.split("\t"))))
    return header, rows


class MapBuilderTests(unittest.TestCase):
    def setUp(self):
        self.root = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.root)

    def chainmap(self, d):
        return read_map(os.path.join(d.engine1, PDB, si.CHAINMAP_NAME))

    def peptide_map(self, d):
        return read_map(os.path.join(d.engine2, PDB, sp.MAP_NAME))

    def test_the_build_steps_resolve_from_the_index(self):
        d = Delivery(self.root)
        d.run(*ligand_imports.MAP_COMMANDS)
        header, rows = self.chainmap(d)
        self.assertEqual(header["cif_sha256"], SHA)
        self.assertEqual(
            (header["receptor.auth_chain"], header["receptor.status"]), ("R", "ok")
        )
        self.assertEqual(
            [
                (r["het"], r["token"], r["instance"], r["status"], r["source"])
                for r in rows
            ],
            [("LIG", "A:401", "LIG_R_401", "ok", "coord+label")],
        )
        header, rows = self.peptide_map(d)
        self.assertEqual(header["cif_sha256"], SHA)
        self.assertEqual(
            (header["receptor.auth_chain"], header["receptor.segment"]), ("R", REC_SEG)
        )
        self.assertEqual(
            [
                (
                    r["gpcrdb_chain"],
                    r["auth_chain"],
                    r["status"],
                    r["items"],
                    r["outcomes"],
                )
                for r in rows
            ],
            [("P", "P", "ok", K_PEP, "done")],
        )

    def test_the_peptide_maps_read_the_index_from_the_engine1_tree(self):
        d = Delivery(self.root)
        self.assertFalse(os.path.exists(cm.index_path(d.engine2, PDB)))
        d.run(ligand_imports.MAP_COMMANDS[1])
        self.assertEqual(self.peptide_map(d)[0]["receptor.status"], "ok")

    def test_no_index_leaves_the_structure_unresolved(self):
        d = Delivery(self.root, index=False)
        d.run(*ligand_imports.MAP_COMMANDS)
        header, rows = self.chainmap(d)
        self.assertEqual(header["receptor.status"], "unresolved")
        self.assertEqual(header["receptor.gpcrdb_text_sha256"], "")
        self.assertEqual([r["status"] for r in rows], ["unresolved"])
        header = self.peptide_map(d)[0]
        self.assertEqual(header["receptor.status"], "unresolved")
        self.assertEqual(header["receptor.gpcrdb_text_sha256"], "")

    def test_a_reformatted_index_is_unresolved_not_matched_by_a_fallback(self):
        d = Delivery(self.root, coord="%.2f")
        d.run(*ligand_imports.MAP_COMMANDS)
        self.assertEqual(self.chainmap(d)[0]["receptor.status"], "unresolved")
        self.assertEqual(self.peptide_map(d)[0]["receptor.status"], "unresolved")

    def test_an_index_from_another_mmcif_than_the_products_is_refused(self):
        d = Delivery(self.root, summary_sha="ab" * 32)
        d.run(*ligand_imports.MAP_COMMANDS)
        header, rows = self.chainmap(d)
        self.assertEqual(header["receptor.status"], "unresolved")
        self.assertIn("another mmCIF", header["receptor.note"])
        text_sha = cm.text_sha256(gpcrdb_text())
        self.assertEqual(header["receptor.gpcrdb_text_sha256"], text_sha)
        instances = si.instance_yaml_paths(d.engine1, PDB)
        self.assertEqual(
            header["receptor.product_instances_sha256"], cm.instances_sha256(instances)
        )
        self.assertEqual([r["status"] for r in rows], ["unresolved"])
        header, _rows = self.peptide_map(d)
        self.assertEqual(header["receptor.status"], "unresolved")
        self.assertIn("another mmCIF", header["receptor.note"])
        self.assertEqual(header["receptor.gpcrdb_text_sha256"], text_sha)
        plan = os.path.join(d.engine2, PDB, sp.PLAN_NAME)
        self.assertEqual(header["plan_sha256"], sp.sha256_file(plan))

    def test_a_stored_text_that_does_not_parse_keeps_its_fingerprint(self):
        d = Delivery(self.root)
        bad = gpcrdb_text() + (
            "%-6s%5d %-4s %3s %1s%4s    %8s%8s%8s  1.00 20.00          %2s\n"
            % ("ATOM", 99, "CA", "ALA", "R", "999", "x", "y", "z", "C")
        )
        write(os.path.join(d.gdata, "structure_data", "pdbs", PDB + ".pdb"), bad)
        d.run(*ligand_imports.MAP_COMMANDS)
        for header in (self.chainmap(d)[0], self.peptide_map(d)[0]):
            self.assertEqual(header["receptor.status"], "unresolved")
            self.assertIn("input unreadable", header["receptor.note"])
            self.assertEqual(header["receptor.gpcrdb_text_sha256"], cm.text_sha256(bad))

    def test_a_summary_without_the_field_is_refused(self):
        d = Delivery(self.root, summary_sha=None)
        d.run(*ligand_imports.MAP_COMMANDS)
        self.assertEqual(self.chainmap(d)[0]["receptor.status"], "unresolved")
        self.assertEqual(self.peptide_map(d)[0]["receptor.status"], "unresolved")

    def test_the_commit_is_recorded_as_given_or_unknown(self):
        d = Delivery(self.root)
        d.run(*ligand_imports.MAP_COMMANDS)
        self.assertEqual(self.chainmap(d)[0]["annotation_commit"], "unknown")
        self.assertEqual(self.peptide_map(d)[0]["annotation_commit"], "unknown")
        self.assertEqual(e1.annotation_commit("2468ad4", d.gdata), "2468ad4")

    def test_the_clean_up_removes_the_maps_and_the_directories_the_build_made(self):
        d = Delivery(self.root)
        d.run(*ligand_imports.MAP_COMMANDS)
        self.assertTrue(
            os.path.isfile(os.path.join(d.engine1, EMPTY, si.CHAINMAP_NAME))
        )
        d.run(ligand_imports.CLEAN_COMMAND)
        self.assertFalse(os.path.exists(os.path.join(d.engine1, PDB, si.CHAINMAP_NAME)))
        self.assertFalse(os.path.exists(os.path.join(d.engine2, PDB, sp.MAP_NAME)))
        self.assertFalse(os.path.exists(os.path.join(d.engine1, EMPTY)))
        # The delivery is left as it came.
        self.assertTrue(os.path.isfile(cm.index_path(d.engine1, PDB)))
        self.assertTrue(os.path.isfile(os.path.join(d.engine1, PDB, "summary.yaml")))
        self.assertTrue(os.path.isfile(os.path.join(d.engine2, PDB, sp.PLAN_NAME)))

    def test_the_clean_up_touches_only_the_maps_and_the_directories_they_emptied(self):
        d = Delivery(self.root)
        os.makedirs(os.path.join(d.engine1, "EMPT"))  # empty, held no map
        write(
            os.path.join(d.engine1, "STRY", "notes.txt"), "x\n"
        )  # holds only a stray file
        outside = os.path.join(self.root, "outside")
        write(os.path.join(outside, si.CHAINMAP_NAME), "keep\n")
        os.symlink(outside, os.path.join(d.engine1, "LINK"))  # a symlinked directory
        d.run(*ligand_imports.MAP_COMMANDS)
        d.run(ligand_imports.CLEAN_COMMAND)
        self.assertTrue(os.path.isdir(os.path.join(d.engine1, "EMPT")))
        self.assertTrue(os.path.isfile(os.path.join(d.engine1, "STRY", "notes.txt")))
        self.assertTrue(os.path.isfile(os.path.join(outside, si.CHAINMAP_NAME)))
        self.assertFalse(os.path.exists(os.path.join(d.engine1, EMPTY)))
        d.run(ligand_imports.CLEAN_COMMAND)  # a second run is harmless
        self.assertTrue(os.path.isdir(os.path.join(d.engine1, "EMPT")))


class AnnotationCommitTests(unittest.TestCase):
    def test_a_given_commit_wins(self):
        self.assertEqual(e1.annotation_commit("abc1234", "/nowhere"), "abc1234")

    def test_no_git_binary_is_unknown(self):
        with mock.patch.object(
            e1.subprocess, "run", side_effect=FileNotFoundError("git")
        ):
            self.assertEqual(e1.annotation_commit(None, "/nowhere"), "unknown")

    def test_a_hanging_git_is_unknown(self):
        with mock.patch.object(
            e1.subprocess, "run", side_effect=e1.subprocess.TimeoutExpired("git", 30)
        ):
            self.assertEqual(e1.annotation_commit(None, "/nowhere"), "unknown")

    def test_only_the_top_level_of_a_repository_answers(self):
        def fake(cmd, **kw):
            out = {"--show-toplevel": "/repo", "--short=7": "1234567"}[cmd[4]]
            return e1.subprocess.CompletedProcess(cmd, 0, stdout=out + "\n")

        with mock.patch.object(e1.subprocess, "run", side_effect=fake):
            self.assertEqual(e1.annotation_commit(None, "/repo/inside"), "unknown")
            self.assertEqual(e1.annotation_commit(None, "/repo"), "1234567")


class ProductInputShaTests(unittest.TestCase):
    def setUp(self):
        self.tree = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.tree)

    def test_read_from_the_summary_none_without_one_empty_without_the_field(self):
        self.assertIsNone(e1.product_input_sha256(self.tree, PDB))
        write(os.path.join(self.tree, PDB, "summary.yaml"), "pdb_id: X\n")
        self.assertEqual(e1.product_input_sha256(self.tree, PDB), "")
        self.assertIn("no input_sha256", e1.index_mismatch("", SHA))
        self.assertIsNone(e1.index_mismatch(None, SHA))
        write(os.path.join(self.tree, PDB, "summary.yaml"), "- not a mapping\n")
        self.assertEqual(e1.product_input_sha256(self.tree, PDB), "")
        write(os.path.join(self.tree, PDB, "summary.yaml"), "a: [\n")
        self.assertEqual(e1.product_input_sha256(self.tree, PDB), "")
        write(os.path.join(self.tree, PDB, "summary.yaml"), "")
        self.assertEqual(e1.product_input_sha256(self.tree, PDB), "")
        with open(os.path.join(self.tree, PDB, "summary.yaml"), "wb") as fh:
            fh.write(b"input_sha256: \xff\xfe\n")
        self.assertEqual(e1.product_input_sha256(self.tree, PDB), "")
        write(os.path.join(self.tree, PDB, "summary.yaml"), "input_sha256: %s\n" % SHA)
        self.assertEqual(e1.product_input_sha256(self.tree, PDB), SHA)

    def test_mismatch_only_when_both_are_known_and_differ(self):
        self.assertIsNone(e1.index_mismatch(None, SHA))
        self.assertIsNone(e1.index_mismatch(SHA, SHA))
        self.assertIn("another mmCIF", e1.index_mismatch("ab" * 32, SHA))


if __name__ == "__main__":
    unittest.main()
