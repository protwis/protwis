"""
Unit tests for the database-free layer of build_schrodinger_chainmap_files.

They need Django configured (the command imports models through
interaction.schrodinger_import) but never touch the database, so they run with
plain unittest:

    python -c "import django; django.setup(); import unittest; \
        unittest.main(module='interaction.test_schrodinger_chainmap_files', argv=['x'])"
"""

import os
import shutil
import tempfile
import unittest

from build.management.commands import build_schrodinger_chainmap_files as b
from interaction import schrodinger_chain_map as cm
from interaction import schrodinger_import as si


def tsv(rows):
    return "".join("\t".join(r) + "\n" for r in rows)


class ReadTsvTests(unittest.TestCase):
    def setUp(self):
        self.root = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.root)

    def write(self, text):
        path = os.path.join(self.root, "x.tsv")
        with open(path, "w", newline="") as fh:
            fh.write(text)
        return path

    def test_values_are_stripped(self):
        path = self.write(tsv([["PDB", " Name "], [" 2rh1 ", " CAU "]]))
        self.assertEqual(b.read_tsv(path), [{"PDB": "2rh1", "Name": "CAU"}])

    def test_a_row_wider_than_the_header_is_refused(self):
        from django.core.management.base import CommandError

        path = self.write(tsv([["PDB", "Name"], ["2RH1", "CAU", "extra"]]))
        with self.assertRaises(CommandError):
            b.read_tsv(path)

    def test_a_short_row_is_padded_and_loses_its_token(self):
        """
        The safe direction: the anchor keeps its component but no residue.

        The map cannot cover the database's token, so the importer refuses the
        structure.
        """
        path = self.write(tsv([["PDB", "Name", "Type"], ["2RH1", "CAU"]]))
        self.assertEqual(b.read_tsv(path), [{"PDB": "2RH1", "Name": "CAU", "Type": ""}])
        self.assertEqual(
            b.annotation_anchors(b.read_tsv(path)), {"2RH1": [("CAU", "", "")]}
        )


class PreferredChainTests(unittest.TestCase):
    def test_the_parser_rule_is_mirrored(self):
        rows = [
            {"PDB": "2rh1", "ChainID": "A"},
            {"PDB": "8JCU", "ChainID": "1.0"},
            {"PDB": "9ZZZ", "ChainID": "A,B"},
        ]
        self.assertEqual(
            b.preferred_chains(rows), {"2RH1": "A", "8JCU": "1", "9ZZZ": "A,B"}
        )

    def test_a_row_without_a_pdb_is_dropped(self):
        self.assertEqual(b.preferred_chains([{"PDB": "", "ChainID": "A"}]), {})


class AnnotationAnchorTests(unittest.TestCase):
    @staticmethod
    def row(**kw):
        r = {
            "PDB": "2RH1",
            "Name": "CAU",
            "Type": "small-molecule",
            "Residue_seq_id": "A:408",
        }
        r.update(kw)
        return r

    def test_every_component_reference_whatever_its_type(self):
        rows = [
            self.row(),
            self.row(Name="OLA", Type="lipid", Residue_seq_id="A:1"),
            self.row(Name="PEP", Type="peptide", Residue_seq_id="B:1"),
            self.row(Name="D2U", Type="peptide", Residue_seq_id="A:1201"),
            self.row(Name="X", Type="protein", Residue_seq_id="B:2"),
            self.row(Name="Y", Type="None", Residue_seq_id="B:3"),
            self.row(Name="Z", Type="small molecule", Residue_seq_id="B:4"),
        ]
        got = b.annotation_anchors(rows)
        self.assertEqual(
            [(h, t) for h, t, _ in got["2RH1"]],
            [
                ("CAU", "A:408"),
                ("OLA", "A:1"),
                ("D2U", "A:1201"),
                ("X", "B:2"),
                ("Y", "B:3"),
                ("Z", "B:4"),
            ],
        )

    def test_the_map_selects_what_the_importer_serves(self):
        """
        annotation_anchors and schrodinger_import.is_in_scope agree.

        Checked on every combination of reference and type.
        """
        import types as pytypes

        references = ["CAU", "D2U", "PEP", "pep", "APO", " apo ", "", None]
        kinds = [
            "small-molecule",
            "lipid",
            "peptide",
            "protein",
            "none",
            "None",
            "small molecule",
            "",
        ]
        for ref in references:
            for kind in kinds:
                sli = pytypes.SimpleNamespace(
                    pdb_reference=ref,
                    ligand=pytypes.SimpleNamespace(
                        ligand_type=pytypes.SimpleNamespace(slug=kind)
                    ),
                )
                row = self.row(Name=ref if ref is not None else "", Type=kind)
                mapped = bool(b.annotation_anchors([row]))
                self.assertEqual(mapped, si.is_in_scope(sli), (ref, kind))

    def test_a_placeholder_name_is_not_an_anchor(self):
        rows = [self.row(Name=name) for name in sorted(si.PLACEHOLDER_REFERENCES)]
        rows += [
            self.row(Name=name.lower()) for name in sorted(si.PLACEHOLDER_REFERENCES)
        ]
        self.assertEqual(b.annotation_anchors(rows), {})

    def test_every_copy_is_kept_and_deduped(self):
        """The annotation splits a ligand into copies; all of them are written."""
        rows = [
            self.row(Residue_seq_id="A:977, B:978"),
            self.row(Residue_seq_id="A:977"),
        ]
        self.assertEqual(
            [(h, t) for h, t, _ in b.annotation_anchors(rows)["2RH1"]],
            [("CAU", "A:977"), ("CAU", "B:978")],
        )

    def test_an_unparsable_chain_res_becomes_the_empty_token(self):
        rows = [self.row(Residue_seq_id="A"), self.row(Residue_seq_id="")]
        got = b.annotation_anchors(rows)["2RH1"]
        self.assertEqual([(h, t) for h, t, _ in got], [("CAU", "")])
        self.assertEqual(got[0][2], "A")  # the raw text is carried for the note


class HeaderTests(unittest.TestCase):
    def setUp(self):
        self.root = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.root)

    def test_flat_replaces_only_what_would_break_the_format(self):
        self.assertEqual(
            b._flat("CA also matched on\tB:12,\nC:13  double"),
            "CA also matched on B:12, C:13  double",
        )
        self.assertEqual(b._flat(None), "")
        self.assertEqual(b._flat(0), "0")

    def test_some_says_what_it_left_out(self):
        self.assertEqual(b._some(["a", "b"]), "a, b")
        self.assertEqual(
            b._some([str(i) for i in range(25)], limit=20),
            ", ".join(str(i) for i in range(20)) + " (+5 more)",
        )

    def test_a_header_value_with_a_tab_still_round_trips(self):
        path = os.path.join(self.root, "chainmap.tsv")
        note = "CA also matched on\tB:12"
        b.write_chainmap(
            path,
            [
                ("schema", si.CHAINMAP_SCHEMA),
                ("pdb", "2RH1"),
                (si.CHAINMAP_RECEPTOR_PREFIX + "note", note),
            ],
            [],
        )
        header, fieldnames, rows = si._read_map(path)
        self.assertEqual(
            header[si.CHAINMAP_RECEPTOR_PREFIX + "note"], "CA also matched on B:12"
        )
        self.assertEqual(tuple(fieldnames), tuple(cm.ANCHOR_COLUMNS))
        self.assertEqual(rows, [])

    def test_the_write_is_atomic(self):
        path = os.path.join(self.root, "chainmap.tsv")
        with open(path + ".tmp", "w") as fh:
            fh.write("leftover from a killed run")
        b.write_chainmap(path, [("schema", si.CHAINMAP_SCHEMA)], [])
        self.assertEqual(sorted(os.listdir(self.root)), ["chainmap.tsv"])


class ProductSummaryTests(unittest.TestCase):
    """The one field that separates "ran and found nothing" from "never ran"."""

    def setUp(self):
        self.root = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.root)
        self.values = dict.fromkeys(b.HEADER_KEYS, "x")

    def test_it_is_read_from_the_product_tree(self):
        os.makedirs(os.path.join(self.root, "2RH1"))
        path = os.path.join(self.root, "2RH1", si.PRODUCT_SUMMARY_NAME)
        self.assertFalse(b.has_product_summary(self.root, "2RH1"))
        with open(path, "w") as fh:
            fh.write("pdb_id: 2RH1\n")
        self.assertTrue(b.has_product_summary(self.root, "2RH1"))
        self.assertFalse(b.has_product_summary(self.root, "6ZIN"))

    def test_an_empty_summary_is_a_truncated_copy_not_a_run(self):
        os.makedirs(os.path.join(self.root, "2RH1"))
        open(os.path.join(self.root, "2RH1", si.PRODUCT_SUMMARY_NAME), "w").close()
        self.assertFalse(b.has_product_summary(self.root, "2RH1"))

    def test_the_header_says_no_when_the_producer_left_nothing(self):
        for has_summary, expected in ((True, "yes"), (False, "no")):
            header = dict(
                b.chainmap_header("2RH1", self.values, {"auth_chain": "A"}, has_summary)
            )
            self.assertEqual(header[si.PRODUCT_SUMMARY_KEY], expected)
            self.assertEqual(header["schema"], si.CHAINMAP_SCHEMA)
            self.assertEqual(header["pdb"], "2RH1")
            self.assertEqual(header[si.CHAINMAP_RECEPTOR_PREFIX + "auth_chain"], "A")

    def test_a_never_run_structure_is_judged_not_run_by_the_reader(self):
        """The writer and the reader have to agree, or the field buys nothing."""
        receptor = dict.fromkeys(cm.RECEPTOR_COLUMNS, "")
        receptor["product_instances_sha256"] = cm.instances_sha256([])
        d = os.path.join(self.root, "6ZIN")
        os.makedirs(d)
        path = os.path.join(d, si.CHAINMAP_NAME)
        for has_summary, expected in ((False, ["6ZIN"]), (True, [])):
            b.write_chainmap(
                path, b.chainmap_header("6ZIN", self.values, receptor, has_summary), []
            )
            not_run = si.load_chainmap_dir(self.root, ["6ZIN"])[4]
            self.assertEqual(not_run, expected)


class BuilderStampTests(unittest.TestCase):
    def test_the_parts_cannot_be_recut(self):
        """One part is a function body and carries newlines."""
        self.assertNotEqual(b._sha256_parts(["ab", "c"]), b._sha256_parts(["a", "bc"]))
        self.assertNotEqual(b._sha256_parts(["a\nb"]), b._sha256_parts(["a", "b"]))
        self.assertEqual(b._sha256_parts(["a", "b"]), b._sha256_parts(["a", "b"]))
