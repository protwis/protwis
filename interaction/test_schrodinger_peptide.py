"""Unit tests for the database-free layer of interaction.schrodinger_peptide.

They need Django configured (the module imports models) but never touch the
database, so they run with plain unittest:

    python -c "import django; django.setup(); import unittest; \
        unittest.main(module='interaction.test_schrodinger_peptide', argv=['x'])"
"""

import json
import os
import shutil
import tempfile
import types
import unittest

from build.management.commands import build_schrodinger_peptide_maps as builder
from interaction import schrodinger_chain_map as cm
from interaction import schrodinger_import as si
from interaction import schrodinger_peptide as sp

# A peptide atom line as the producer writes it (no altloc column).
GLN = "HETATM   93  NE2GLN D  19      50.943 -11.924  33.028  1.00 83.04           N"


def producer_line(name, resname, chain, resnum, x=1.0, b=20.0, element="C", serial=1):
    """Build a line in the producer layout: one column short of standard PDB."""
    return "HETATM{:>5} {:<4}{:>3} {}{:>4}     {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}           {}".format(
        serial,
        (" " + name) if len(name) < 4 else name,
        resname,
        chain,
        resnum,
        x,
        2.0,
        3.0,
        1.0,
        b,
        element,
    )


def row(family, direction="", seq=147, aa="D", chain="R", atom="OD1", lig=""):
    return {
        "feature_family": family,
        "direction": direction,
        "receptor_atom_name": atom,
        "receptor_pdb_block": "ATOM",
        "ligand_pdb_block": lig,
        "receptor_residue": {
            "name_1_letter": aa,
            "pdb_residue_number": seq,
            "chain_id": chain,
            "insertion_code": "",
        },
    }


def fake_sli(reference, ltype):
    return types.SimpleNamespace(
        pdb_reference=reference,
        ligand=types.SimpleNamespace(ligand_type=types.SimpleNamespace(slug=ltype)),
    )


class ScopeTests(unittest.TestCase):
    def test_every_pep_anchor_whatever_its_type(self):
        self.assertTrue(sp.is_in_scope(fake_sli("pep", "peptide")))
        self.assertTrue(sp.is_in_scope(fake_sli("PEP", "small-molecule")))
        self.assertTrue(sp.is_in_scope(fake_sli("pep", "protein")))
        self.assertTrue(sp.is_in_scope(fake_sli(" Pep ", "protein")))
        for ltype in ("lipid", "na", "", "anything"):
            self.assertTrue(sp.is_in_scope(fake_sli("pep", ltype)), ltype)
        self.assertFalse(sp.is_in_scope(fake_sli("ZMA", "peptide")))
        self.assertFalse(sp.is_in_scope(fake_sli("", "peptide")))
        self.assertFalse(sp.is_in_scope(fake_sli(None, "peptide")))

    def test_no_anchor_is_served_by_both_imports(self):
        for ref in ("pep", "PEP", "ZMA", "apo", ""):
            for ltype in ("peptide", "small-molecule", "lipid", "protein"):
                sli = fake_sli(ref, ltype)
                self.assertFalse(
                    sp.is_in_scope(sli) and si.is_in_scope(sli), (ref, ltype)
                )


def seg(name, chain, ranges, accession=""):
    return {
        "name": name,
        "chain_id": chain,
        "ranges": ranges,
        "db_accession": accession,
    }


class SegmentTests(unittest.TestCase):
    SEGMENTS = {
        "A_OPRD_1_338": seg("A_OPRD_1_338", "A", [[1, 338]], "4RWA"),
        "A_BRIL_1001_1106": seg("A_BRIL_1001_1106", "A", [[1001, 1106]], "P0ABE7"),
        "H_PEP_1_5": seg("H_PEP_1_5", "H", [[1, 5]]),
        "H_seg_6_6": seg("H_seg_6_6", "H", [[6, 6]]),
    }

    def test_segment_residues(self):
        self.assertEqual(
            sp.segment_residues(seg("x", "A", [[3, 5], [9, 9]])), {3, 4, 5, 9}
        )

    def test_the_receptor_is_the_segment_covering_most_receptor_residues(self):
        name, covered, note = sp.receptor_segment(
            self.SEGMENTS, "A", set(range(40, 300)) | {1050}
        )
        self.assertEqual((name, covered, note), ("A_OPRD_1_338", 260, ""))

    def test_no_coverage_is_refused(self):
        one = {"A_R_1_300": seg("A_R_1_300", "A", [[1, 300]])}
        name, covered, note = sp.receptor_segment(one, "A", {2000})
        self.assertEqual((name, covered), (None, 0))
        self.assertIn("covers a receptor residue", note)
        self.assertIsNone(sp.receptor_segment(self.SEGMENTS, "B", {10})[0])

    def test_a_tie_is_refused(self):
        name, covered, note = sp.receptor_segment(self.SEGMENTS, "A", {10, 1010})
        self.assertEqual((name, covered), (None, 1))
        self.assertIn("equally", note)

    def test_the_receptor_side_is_every_segment_of_the_receptor_chain(self):
        self.assertEqual(
            sp.receptor_chain_segments(self.SEGMENTS, "A"),
            ["A_BRIL_1001_1106", "A_OPRD_1_338"],
        )
        reordered = dict(reversed(list(self.SEGMENTS.items())))
        self.assertEqual(
            sp.receptor_chain_segments(reordered, "A"),
            ["A_BRIL_1001_1106", "A_OPRD_1_338"],
        )

    def test_peptide_items_take_every_segment_of_the_chain_as_the_ligand_side(self):
        items = [
            {
                "key": "k1",
                "ligand_segment": "H_PEP_1_5",
                "receptor_segment": "A_OPRD_1_338",
            },
            {
                "key": "k2",
                "ligand_segment": "A_OPRD_1_338",
                "receptor_segment": "H_PEP_1_5",
            },
            {
                "key": "k3",
                "ligand_segment": "H_seg_6_6",
                "receptor_segment": "A_OPRD_1_338",
            },
            {
                "key": "k4",
                "ligand_segment": "H_PEP_1_5",
                "receptor_segment": "A_BRIL_1001_1106",
            },
        ]
        self.assertEqual(
            sp.peptide_items(self.SEGMENTS, items, "H", {"A_OPRD_1_338"}),
            [("H_PEP_1_5", "k1"), ("H_seg_6_6", "k3")],
        )
        # A receptor declared in two pieces: items against both.
        self.assertEqual(
            sp.peptide_items(
                self.SEGMENTS, items, "H", {"A_OPRD_1_338", "A_BRIL_1001_1106"}
            ),
            [("H_PEP_1_5", "k1"), ("H_PEP_1_5", "k4"), ("H_seg_6_6", "k3")],
        )

    def test_the_receptor_segment_is_never_its_own_peptide(self):
        segs = {
            "A_R_1_300": seg("A_R_1_300", "A", [[1, 300]]),
            "A_PEP_401_410": seg("A_PEP_401_410", "A", [[401, 410]]),
        }
        items = [
            {
                "key": "k",
                "ligand_segment": "A_PEP_401_410",
                "receptor_segment": "A_R_1_300",
            },
            {
                "key": "self",
                "ligand_segment": "A_R_1_300",
                "receptor_segment": "A_R_1_300",
            },
        ]
        self.assertEqual(
            sp.peptide_items(segs, items, "A", {"A_R_1_300"}), [("A_PEP_401_410", "k")]
        )
        # When every segment of the chain is receptor, the peptide has none.
        self.assertEqual(
            sp.peptide_items(segs, items, "A", {"A_R_1_300", "A_PEP_401_410"}), []
        )


def cif_atom(chain, seq, atom, key, group="ATOM", comp="ALA"):
    return {
        "label_asym": chain,
        "auth_asym": chain,
        "comp": comp,
        "auth_seq": str(seq),
        "icode": "",
        "atom": atom,
        "group": group,
        "key": key,
    }


def g_atom(chain, resnum, atom, key, group="ATOM", resname="ALA"):
    return {
        "chain": chain,
        "resnum": str(resnum),
        "icode": "",
        "resname": resname,
        "atom": atom,
        "group": group,
        "key": key,
    }


class PeptideChainTests(unittest.TestCase):
    def test_ca_coordinates_name_the_author_chain(self):
        cif = [cif_atom("P", 1, "CA", "1 1 1"), cif_atom("P", 2, "CA", "2 2 2")]
        gp = [g_atom("D", 1, "CA", "1 1 1"), g_atom("D", 2, "CA", "2 2 2")]
        self.assertEqual(sp.peptide_author_chain("X", "D", cif, gp)[:2], ("P", "ca"))

    def test_an_all_hetatm_peptide_falls_back_to_every_atom(self):
        cif = [
            cif_atom("C", 1, "CA", "1 1 1", group="HETATM", comp="DAL"),
            cif_atom("C", 1, "N", "1 1 2", group="HETATM", comp="DAL"),
        ]
        gp = [
            g_atom("C", 1, "CA", "1 1 1", group="HETATM", resname="DAL"),
            g_atom("C", 1, "N", "1 1 2", group="HETATM", resname="DAL"),
        ]
        auth, method, note = sp.peptide_author_chain("X", "C", cif, gp)
        self.assertEqual((auth, method, note), ("C", "any_atom", "2/2 atoms"))

    def test_a_weak_any_atom_match_is_refused(self):
        gp = [
            g_atom("C", 1, a, "k%d" % i, group="HETATM") for i, a in enumerate("ABCD")
        ]
        weak = [cif_atom("C", 1, "A", "k0", group="HETATM")]
        split = [
            cif_atom("C", 1, "A", "k0", group="HETATM"),
            cif_atom("C", 1, "B", "k1", group="HETATM"),
            cif_atom("E", 1, "C", "k2", group="HETATM"),
            cif_atom("E", 1, "D", "k3", group="HETATM"),
        ]
        self.assertEqual(sp.peptide_author_chain("X", "C", weak, gp)[0], "")
        self.assertEqual(sp.peptide_author_chain("X", "C", split, gp)[0], "")
        self.assertEqual(sp.peptide_author_chain("X", "Z", weak, gp)[0], "")

    def test_receptor_ca_numbers_are_author_numbers_of_shared_ca(self):
        cif = [
            cif_atom("R", 147, "CA", "a"),
            cif_atom("R", 148, "CA", "b"),
            cif_atom("R", 149, "CB", "c"),
        ]
        gp = [g_atom("A", 47, "CA", "a"), g_atom("A", 49, "CB", "c")]
        self.assertEqual(sp.receptor_ca_numbers("R", "A", cif, gp), {147})


class PeptideLineTests(unittest.TestCase):
    def test_the_producer_layout_is_read(self):
        a = sp.parse_peptide_line(GLN, "D")
        self.assertEqual(
            (a["name"], a["resname"], a["resnum"], a["icode"], a["element"]),
            ("NE2", "GLN", 19, "", "N"),
        )
        self.assertEqual(
            (a["x"], a["y"], a["z"], a["b"]), (50.943, -11.924, 33.028, 83.04)
        )

    def test_prepared_names_become_standard_ones(self):
        for prepared, standard in (
            ("HID", "HIS"),
            ("HIE", "HIS"),
            ("HIP", "HIS"),
            ("CYX", "CYS"),
            ("ASH", "ASP"),
            ("GLH", "GLU"),
            ("LYN", "LYS"),
            ("ARN", "ARG"),
        ):
            a = sp.parse_peptide_line(producer_line("CA", prepared, "D", 3), "D")
            self.assertEqual(a["resname"], standard)
        self.assertEqual(len(sp.PREPARED_NAMES), 8)

    def test_negative_numbers_and_long_names(self):
        a = sp.parse_peptide_line(producer_line("C1", "A1ABC", "D", -3), "D")
        self.assertEqual((a["resname"], a["resnum"]), ("A1ABC", -3))

    def test_another_chain_or_garbage_is_refused(self):
        with self.assertRaises(sp.MalformedPeptideLine):
            sp.parse_peptide_line(GLN, "E")
        with self.assertRaises(sp.MalformedPeptideLine):
            sp.parse_peptide_line("REMARK nothing here", "D")
        with self.assertRaises(sp.MalformedPeptideLine):
            sp.parse_peptide_line(GLN[:50], "D")

    def test_standard_columns_and_gpcrdb_chain(self):
        line, capped = sp.standard_peptide_line(sp.parse_peptide_line(GLN, "D"), "P")
        self.assertEqual(len(line), 78)
        self.assertEqual(
            (line[:6], line[12:16], line[17:20], line[21], line[22:26]),
            ("HETATM", " NE2", "GLN", "P", "  19"),
        )
        self.assertFalse(capped)

    def test_every_standard_column_with_a_long_name_and_a_long_chain(self):
        line = producer_line(
            "NZ", "A1D5B", "CCC", 1007, x=-12.345, b=33.5, element="N", serial=99999
        )
        a = sp.parse_peptide_line(line, "CCC")
        out, capped = sp.standard_peptide_line(a, "C")
        self.assertEqual(len(out), 78)
        self.assertEqual(
            (
                out[0:6],
                out[6:11],
                out[12:16],
                out[16],
                out[17:20],
                out[21],
                out[22:26],
                out[26],
            ),
            ("HETATM", "99999", " NZ ", " ", "A1D", "C", "1007", " "),
        )
        self.assertEqual(
            (out[30:38], out[38:46], out[46:54], out[54:60], out[60:66], out[76:78]),
            (" -12.345", "   2.000", "   3.000", "  1.00", " 33.50", " N"),
        )
        self.assertFalse(capped)

    def test_a_large_b_factor_is_capped(self):
        line, capped = sp.standard_peptide_line(
            sp.parse_peptide_line(producer_line("CA", "GLY", "D", 1, b=1500.0), "D"),
            "D",
        )
        self.assertTrue(capped)
        self.assertEqual(line[60:66], "999.99")

    def test_standardise_blocks_rewrites_in_place(self):
        rows = [row("HPhob", lig=GLN + "\n" + producer_line("CB", "GLN", "D", 19))]
        sp.standardise_blocks(rows, "D", "P")
        self.assertEqual(
            [line[21] for line in rows[0]["ligand_pdb_block"].splitlines()], ["P", "P"]
        )


class PeptideTypeTests(unittest.TestCase):
    def test_vocabulary_names_the_receptor_first(self):
        self.assertEqual(
            sp.peptide_type(row("Donor", "ligand-donor"))[:2],
            ("polar", "acceptor-donor"),
        )
        self.assertEqual(
            sp.peptide_type(row("Acceptor", "ligand-acceptor"))[:2],
            ("polar", "donor-acceptor"),
        )
        self.assertEqual(
            sp.peptide_type(row("NegCharge", "neg-pos"))[:2],
            ("ionic", "positive-negative"),
        )
        self.assertEqual(
            sp.peptide_type(row("PosCharge", "pos-neg"))[:2],
            ("ionic", "negative-positive"),
        )
        self.assertEqual(sp.peptide_type(row("HPhob"))[:2], ("hydrophobic", ""))
        self.assertEqual(sp.peptide_type(row("VdW"))[:2], ("van-der-waals", ""))

    def test_ring_sides(self):
        self.assertEqual(
            sp.peptide_type(row("PiCat", "ligand-cation")),
            ("aromatic", "pi-cation", True, False),
        )
        self.assertEqual(
            sp.peptide_type(row("PiCat", "receptor-cation")),
            ("aromatic", "cation-pi", False, True),
        )
        self.assertEqual(
            sp.peptide_type(row("Aromatic", "edge-to-face")),
            ("aromatic", "edge-to-face", True, True),
        )
        self.assertEqual(
            sp.peptide_type(row("Aromatic", "face-to-face"))[:2],
            ("aromatic", "face-to-face"),
        )

    def test_families_without_a_type_and_unknown_ones(self):
        for family in ("Accessible", "Covalent", "XBond", "Metal", "Wat-HBond"):
            self.assertIsNone(sp.peptide_type(row(family)))
        with self.assertRaises(si.UnroutableRow):
            sp.peptide_type(row("Donor", "sideways"))
        with self.assertRaises(si.UnroutableRow):
            sp.peptide_type(row("Mystery"))

    def test_every_rfi_family_is_known_to_the_peptide_vocabulary(self):
        # A family the RFI side routes must either have a peptide type or be
        # explicitly not in the peptide tables; nothing falls between.
        for family, direction in si.load_type_map():
            known = (
                family,
                direction,
            ) in sp.PEPTIDE_TYPES or family in sp.NOT_IN_PEPTIDE_TABLES
            self.assertTrue(known, (family, direction))


class PlanPeptidePairsTests(unittest.TestCase):
    def test_pairs_atoms_and_accounting(self):
        tyr_n = producer_line("N", "TYR", "D", 1, serial=1, element="N")
        tyr_ring = "\n".join(
            producer_line(n, "TYR", "D", 1, serial=i + 2)
            for i, n in enumerate(("CG", "CD1", "CD2", "CE1", "CE2", "CZ"))
        )
        rows = [
            row("PosCharge", "pos-neg", lig=tyr_n),  # salt bridge to D147
            row("Donor", "ligand-donor", lig=tyr_n),  # same pair, polar
            row("Donor", "ligand-donor", lig=tyr_n),  # duplicate
            row("Aromatic", "edge-to-face", seq=293, aa="W", atom="CG", lig=tyr_ring),
            row("Accessible", lig=tyr_n),  # not in the tables
            row("HPhob", aa="X", lig=tyr_n),  # nonstandard receptor
            row("HPhob", chain="B", lig=tyr_n),  # another chain
        ]
        pairs, counts = sp.plan_peptide_pairs(rows, "R", "D")
        self.assertEqual(
            counts["rows_in"],
            counts["not_in_peptide_tables"]
            + counts["nonstandard_residue"]
            + counts["other_chain"]
            + counts["used"],
        )
        self.assertEqual(
            (
                counts["not_in_peptide_tables"],
                counts["nonstandard_residue"],
                counts["other_chain"],
                counts["used"],
            ),
            (1, 1, 1, 4),
        )
        self.assertEqual(
            pairs[(1, "", "TYR", 147, "D")],
            [
                ("N", "OD1", "ionic", "negative-positive"),
                ("N", "OD1", "polar", "acceptor-donor"),
            ],
        )
        self.assertEqual(
            pairs[(1, "", "TYR", 293, "W")],
            [("RN1", "RN1", "aromatic", "edge-to-face")],
        )
        self.assertEqual((counts["pairs"], counts["interactions"]), (2, 3))

    def test_a_row_spanning_two_residues_makes_two_pairs(self):
        lig = (
            producer_line("CB", "ALA", "D", 2)
            + "\n"
            + producer_line("CB", "LEU", "D", 3)
        )
        pairs, _ = sp.plan_peptide_pairs(
            [row("HPhob", seq=100, aa="F", atom="CZ", lig=lig)], "R", "D"
        )
        self.assertEqual(
            sorted(pairs), [(2, "", "ALA", 100, "F"), (3, "", "LEU", 100, "F")]
        )

    def test_rings_and_cations_through_the_planner(self):
        ring = "\n".join(
            producer_line(n, "PHE", "D", 4, serial=i)
            for i, n in enumerate(("CG", "CD1", "CZ"))
        )
        nz = producer_line("NZ", "LYS", "D", 5, element="N")
        rows = [
            row("PiCat", "receptor-cation", seq=200, aa="R", atom="NH1", lig=ring),
            row("PiCat", "ligand-cation", seq=300, aa="W", atom="CD2", lig=nz),
            row("Aromatic", "face-to-face", seq=310, aa="F", atom="CG", lig=ring),
        ]
        pairs, _ = sp.plan_peptide_pairs(rows, "R", "D")
        self.assertEqual(
            pairs[(4, "", "PHE", 200, "R")], [("RN1", "NH1", "aromatic", "cation-pi")]
        )
        self.assertEqual(
            pairs[(5, "", "LYS", 300, "W")], [("NZ", "RN1", "aromatic", "pi-cation")]
        )
        self.assertEqual(
            pairs[(4, "", "PHE", 310, "F")],
            [("RN1", "RN1", "aromatic", "face-to-face")],
        )

    def test_a_residue_with_an_insertion_code_makes_no_pair(self):
        plain = producer_line("CZ2", "TRP", "D", 100)
        coded = plain[:25] + "C" + plain[26:]
        self.assertEqual(sp.parse_peptide_line(coded, "D")["icode"], "C")
        ser = producer_line("OG", "SER", "D", 100, element="O")
        pairs, counts = sp.plan_peptide_pairs(
            [row("HPhob", seq=185, aa="F", atom="CZ", lig=coded + "\n" + ser)], "R", "D"
        )
        self.assertEqual(sorted(pairs), [(100, "", "SER", 185, "F")])
        self.assertEqual(counts["insertion_code_atoms"], 1)
        pairs, counts = sp.plan_peptide_pairs(
            [row("HPhob", seq=185, aa="F", atom="CZ", lig=coded)], "R", "D"
        )
        self.assertEqual(
            (pairs, counts["used"], counts["insertion_code_atoms"]), ({}, 1, 1)
        )

    def test_the_level_is_the_normal_definition(self):
        self.assertEqual(sp.LEVEL, 0)

    def test_an_unroutable_row_raises_before_anything_is_planned(self):
        with self.assertRaises(si.UnroutableRow):
            sp.plan_peptide_pairs(
                [row("HPhob", lig=GLN), row("Donor", "odd", lig=GLN)], "R", "D"
            )


class MapFileTests(unittest.TestCase):
    def setUp(self):
        self.dir = tempfile.mkdtemp()
        self.path = os.path.join(self.dir, "peptide_map.tsv")
        self.receptor = {
            "preferred_chain": "R",
            "auth_chain": "R",
            "status": "ok",
            "method": "exact",
            "segment": "R_OPRM_1_300",
            "segments": "R_OPRM_1_300,R_seg_24_26",
            "n_ca_gpcrdb": 281,
            "n_ca_matched": 281,
            "n_covered": 281,
            "gpcrdb_text_sha256": "abc",
            "note": "",
        }
        self.rows = [
            dict(
                {c: "" for c in sp.MAP_COLUMNS},
                pdb="6DDF",
                gpcrdb_chain="D",
                auth_chain="D",
                chain_method="ca",
                status="ok",
                segments="D_X_1_5,D_seg_6_6",
                items="k1,k2",
                outcomes="done,selections_apart",
            ),
            dict(
                {c: "" for c in sp.MAP_COLUMNS},
                pdb="6DDF",
                gpcrdb_chain="E",
                status="chain_unresolved",
                note="no atom",
            ),
        ]

    def tearDown(self):
        shutil.rmtree(self.dir)

    def write(self, rows=None, header=None):
        sp.write_peptide_map(
            self.path,
            header or [("pdb", "6DDF"), ("annotation_commit", "abc")],
            self.receptor,
            self.rows if rows is None else rows,
        )

    def test_round_trip(self):
        self.write()
        pdb, receptor, table, prov, header = sp.load_peptide_map(self.path)
        self.assertEqual(pdb, "6DDF")
        self.assertEqual(receptor["segment"], "R_OPRM_1_300")
        self.assertEqual(receptor["segments_list"], ["R_OPRM_1_300", "R_seg_24_26"])
        self.assertEqual(header["annotation_commit"], "abc")
        self.assertEqual(table["D"]["items_list"], ["k1", "k2"])
        self.assertEqual(table["D"]["outcomes_list"], ["done", "selections_apart"])
        self.assertEqual(table["E"]["items_list"], [])
        self.assertEqual(prov["annotation_commit"], "abc")

    def test_refusals(self):
        bad_status = [dict(self.rows[0], status="maybe")]
        mismatched = [dict(self.rows[0], outcomes="done")]
        ok_without_items = [dict(self.rows[0], items="", outcomes="")]
        items_without_ok = [dict(self.rows[1], items="k", outcomes="done")]
        duplicate = [self.rows[0], dict(self.rows[0])]
        other_pdb = [dict(self.rows[0], pdb="1ABC")]
        for rows in (
            bad_status,
            mismatched,
            ok_without_items,
            items_without_ok,
            duplicate,
            other_pdb,
        ):
            self.write(rows)
            with self.assertRaises(sp.MapMismatch):
                sp.load_peptide_map(self.path)

    def test_a_header_without_pdb_is_refused(self):
        self.write(rows=[], header=[("annotation_commit", "abc")])
        with self.assertRaisesRegex(sp.MapMismatch, "names no pdb"):
            sp.load_peptide_map(self.path)

    def test_a_missing_receptor_key_is_refused(self):
        self.write()
        with open(self.path) as fh:
            text = fh.read().replace("# receptor.gpcrdb_text_sha256\tabc\n", "")
        with open(self.path, "w") as fh:
            fh.write(text)
        with self.assertRaisesRegex(sp.MapMismatch, "gpcrdb_text_sha256"):
            sp.load_peptide_map(self.path)

    def test_other_columns_are_refused(self):
        self.write()
        with open(self.path) as fh:
            text = fh.read().replace("\toutcomes\t", "\toutcome\t")
        with open(self.path, "w") as fh:
            fh.write(text)
        with self.assertRaisesRegex(sp.MapMismatch, "columns"):
            sp.load_peptide_map(self.path)

    def test_another_schema_is_refused(self):
        self.write()
        with open(self.path) as fh:
            text = fh.read().replace(sp.MAP_SCHEMA, "engine2-peptide-map/0")
        with open(self.path, "w") as fh:
            fh.write(text)
        with self.assertRaises(sp.MapMismatch):
            sp.load_peptide_map(self.path)


class TreeTests(unittest.TestCase):
    def setUp(self):
        self.dir = tempfile.mkdtemp()
        os.makedirs(os.path.join(self.dir, "6DDF"))

    def tearDown(self):
        shutil.rmtree(self.dir)

    def item(self, key, outcome, rows=None, record_key=None, contract="engine2/3.0"):
        d = os.path.join(self.dir, "6DDF", key)
        os.makedirs(d)
        with open(os.path.join(d, key + ".json"), "w") as fh:
            json.dump(
                {
                    "work_item_key": record_key or key,
                    "outcome": outcome,
                    "provenance": {"contract_version": contract},
                },
                fh,
            )
        if rows is not None:
            with open(os.path.join(d, key + ".yaml"), "w") as fh:
                fh.write("result:\n  interactions: {}\n".format(json.dumps(rows)))

    def map_row(self, keys, outcomes):
        return {"items_list": keys, "outcomes_list": outcomes}

    def test_rows_of_done_items_and_nothing_from_no_interface_answers(self):
        self.item("k1", "done", [row("HPhob", lig=GLN)])
        self.item("k2", "selections_apart")
        rows, failed = sp.anchor_rows(
            self.dir, "6DDF", self.map_row(["k1", "k2"], ["done", "selections_apart"])
        )
        self.assertEqual((len(rows), failed), (1, []))

    def test_a_record_that_disagrees_with_the_map_raises(self):
        self.item("k1", "selections_apart")
        with self.assertRaises(sp.MalformedProduct):
            sp.anchor_rows(self.dir, "6DDF", self.map_row(["k1"], ["done"]))
        self.item("k2", "done", [], record_key="other")
        with self.assertRaises(sp.MalformedProduct):
            sp.anchor_rows(self.dir, "6DDF", self.map_row(["k2"], ["done"]))

    def test_an_open_question_raises(self):
        for n, outcome in enumerate(("started", "eligible", "done_maybe", "")):
            key = "k%d" % n
            self.item(key, outcome)
            with self.assertRaisesRegex(
                sp.MalformedProduct, "leaves the question open"
            ):
                sp.anchor_rows(self.dir, "6DDF", self.map_row([key], [outcome]))

    def test_a_failed_item_is_reported_and_gives_no_row(self):
        self.assertEqual(
            sp.FAILED, {"preparation_failed", "compute_failed", "timed_out", "crashed"}
        )
        for n, outcome in enumerate(sorted(sp.FAILED)):
            key = "f%d" % n
            self.item(key, outcome)
            self.assertEqual(
                sp.anchor_rows(self.dir, "6DDF", self.map_row([key], [outcome])),
                ([], [key + ":" + outcome]),
            )

    def test_one_failed_item_drops_the_rows_of_the_done_ones(self):
        self.item("k1", "done", [row("HPhob", lig=GLN)])
        self.item("k2", "preparation_failed")
        self.item("k3", "selections_apart")
        self.assertEqual(
            sp.anchor_rows(
                self.dir,
                "6DDF",
                self.map_row(
                    ["k1", "k2", "k3"],
                    ["done", "preparation_failed", "selections_apart"],
                ),
            ),
            ([], ["k2:preparation_failed"]),
        )

    def test_a_failed_record_must_agree_with_the_map(self):
        self.item("k1", "done", [row("HPhob", lig=GLN)])
        with self.assertRaisesRegex(sp.MalformedProduct, "the map says"):
            sp.anchor_rows(self.dir, "6DDF", self.map_row(["k1"], ["compute_failed"]))
        self.item("k2", "compute_failed")
        with self.assertRaisesRegex(sp.MalformedProduct, "the map says"):
            sp.anchor_rows(self.dir, "6DDF", self.map_row(["k2"], ["done"]))

    def test_a_record_done_where_the_map_says_no_interface_raises(self):
        self.item("k1", "done", [row("HPhob", lig=GLN)])
        with self.assertRaisesRegex(sp.MalformedProduct, "the map says"):
            sp.anchor_rows(self.dir, "6DDF", self.map_row(["k1"], ["selections_apart"]))

    def test_fingerprints(self):
        plan = os.path.join(self.dir, "6DDF", "plan.json")
        with open(plan, "w") as fh:
            fh.write("{}")
        receptor = {"gpcrdb_text_sha256": cm.text_sha256("ATOM text\n")}
        header = {"plan_sha256": sp.sha256_file(plan)}
        sp.check_fingerprints("6DDF", self.dir, receptor, header, "ATOM text\n")
        with self.assertRaisesRegex(sp.MapMismatch, "structure text"):
            sp.check_fingerprints("6DDF", self.dir, receptor, header, "ATOM other\n")
        with self.assertRaisesRegex(sp.MapMismatch, "plan.json"):
            sp.check_fingerprints(
                "6DDF", self.dir, receptor, {"plan_sha256": "0" * 64}, "ATOM text\n"
            )

    def test_a_stale_map_is_reported_before_its_receptor_row(self):
        plan = os.path.join(self.dir, "6DDF", "plan.json")
        with open(plan, "w") as fh:
            fh.write("{}")
        header = {"plan_sha256": sp.sha256_file(plan)}
        good = {
            "status": "ok",
            "segment": "R_1",
            "segments_list": ["R_1"],
            "note": "",
            "gpcrdb_text_sha256": cm.text_sha256("ATOM text\n"),
        }
        check = sp.check_receptor_and_fingerprints
        check("6DDF", self.dir, good, header, "ATOM text\n")
        renumbered = dict(good, status="renumbered", note="numbering differs")
        with self.assertRaisesRegex(sp.MapMismatch, "structure text"):
            check("6DDF", self.dir, renumbered, header, "ATOM other\n")
        with self.assertRaisesRegex(sp.MapMismatch, "numbering differs"):
            check("6DDF", self.dir, renumbered, header, "ATOM text\n")
        gave_up = dict(good, status="unresolved", note="input unreadable")
        gave_up.update(gpcrdb_text_sha256="")
        with self.assertRaisesRegex(sp.MapMismatch, "input unreadable"):
            check("6DDF", self.dir, gave_up, header, "ATOM other\n")
        no_text = dict(good, gpcrdb_text_sha256="")
        with self.assertRaisesRegex(sp.MapMismatch, "structure text"):
            check("6DDF", self.dir, no_text, header, "ATOM text\n")
        other_plan = {"plan_sha256": "0" * 64}
        with self.assertRaisesRegex(sp.MapMismatch, "plan.json"):
            check("6DDF", self.dir, renumbered, other_plan, "ATOM text\n")

    def test_a_record_of_another_contract_version_raises(self):
        self.item("k1", "done", [row("HPhob", lig=GLN)], contract="engine2/2.0")
        with self.assertRaisesRegex(sp.MalformedProduct, "contract_version"):
            sp.anchor_rows(self.dir, "6DDF", self.map_row(["k1"], ["done"]))
        self.item("k2", "done", [row("HPhob", lig=GLN)], contract=None)
        with self.assertRaisesRegex(sp.MalformedProduct, "contract_version"):
            sp.anchor_rows(self.dir, "6DDF", self.map_row(["k2"], ["done"]))

    def test_a_done_item_without_its_yaml_raises(self):
        self.item("k1", "done")
        with self.assertRaises(si.MalformedProduct):
            sp.anchor_rows(self.dir, "6DDF", self.map_row(["k1"], ["done"]))

    def test_load_plan(self):
        plan = {
            "contract_version": "engine2/3.0",
            "segments": [seg("A", "A", [[1, 2]])],
            "items": [{"key": "k", "ligand_segment": "A", "receptor_segment": "B"}],
        }
        with open(os.path.join(self.dir, "6DDF", "plan.json"), "w") as fh:
            json.dump(plan, fh)
        segments, items = sp.load_plan(self.dir, "6DDF")
        self.assertEqual((list(segments), len(items)), (["A"], 1))
        v = {"contract_version": "engine2/3.0"}
        other = dict(plan, contract_version="engine2/2.0")
        unversioned = {k: plan[k] for k in ("segments", "items")}
        for broken in (
            dict(v, segments=[]),
            dict(v, segments=[seg("A", "A", [])] * 2, items=[]),
            dict(v, segments=[], items=[{"key": "k"}]),
            [1],
            other,
            unversioned,
        ):
            with open(os.path.join(self.dir, "6DDF", "plan.json"), "w") as fh:
                json.dump(broken, fh)
            with self.assertRaises(sp.MalformedProduct):
                sp.load_plan(self.dir, "6DDF")


class OutcomeLogTests(unittest.TestCase):
    def test_an_anchor_without_product_is_a_cleared_warning(self):
        from build.management.commands.import_schrodinger_peptides import Command

        logged = []

        class _Log:
            def log(self, *args, **kwargs):
                logged.append((args, kwargs))

        for mode in ("no_product", "cleared"):
            del logged[:]
            o = sp.AnchorOutcome(2471, "P")
            o.mode = mode
            o.rfi_deleted, o.pairs_deleted = 101, 36
            o.notes.append("items failed: k:preparation_failed")
            Command._log_outcome(_Log(), "9BUD", o)
            self.assertEqual(len(logged), 1, mode)
            args, kwargs = logged[0]
            self.assertEqual(
                args[:5], ("9BUD", "WARNING", "anchor_cleared", 2471, "P"), mode
            )
            self.assertIn("preparation_failed", kwargs["detail"])


class AnchorDecisionTests(unittest.TestCase):
    ROWS = {
        "D": {"status": "ok", "note": ""},
        "E": {"status": "chain_unresolved", "note": "no atom"},
        "F": {"status": "no_items", "note": "none"},
    }

    def test_an_anchor_without_a_chain_is_cleared(self):
        self.assertEqual(sp.anchor_action("X", 1, "", self.ROWS), ("clear", None))

    def test_an_ok_row_is_imported(self):
        self.assertEqual(
            sp.anchor_action("X", 1, "D", self.ROWS), ("import", self.ROWS["D"])
        )

    def test_a_row_that_is_not_ok_fails_the_structure(self):
        for chain in ("E", "F"):
            with self.assertRaises(sp.UnresolvedAnchor):
                sp.anchor_action("X", 1, chain, self.ROWS)
        self.assertTrue(issubclass(sp.UnresolvedAnchor, si.UnresolvedAnchor))

    def test_a_chain_the_map_does_not_list_fails_the_structure(self):
        with self.assertRaisesRegex(sp.MapMismatch, "no row"):
            sp.anchor_action("X", 1, "Z", self.ROWS)

    def test_two_anchors_cannot_share_a_peptide_structure(self):
        used = {}
        sp.claim_peptide_structure(used, 7, 101, "X")
        sp.claim_peptide_structure(used, 8, 102, "X")
        with self.assertRaisesRegex(sp.MissingPeptideStructure, "share"):
            sp.claim_peptide_structure(used, 7, 103, "X")

    def test_the_peptide_structure_of_the_anchor_chain(self):
        a, b = types.SimpleNamespace(chain="D"), types.SimpleNamespace(chain="E")
        self.assertIs(sp.choose_peptide_structure([a], "Q", "x"), a)
        self.assertIs(sp.choose_peptide_structure([a, b], "E", "x"), b)
        for candidates, chain in (
            ([], "D"),
            ([a, b], "Q"),
            ([a, types.SimpleNamespace(chain="D")], "D"),
        ):
            with self.assertRaises(sp.MissingPeptideStructure):
                sp.choose_peptide_structure(candidates, chain, "x")

    def test_the_receptor_must_be_resolved_to_a_listed_primary_segment(self):
        good = {
            "status": "ok",
            "segment": "R_1",
            "segments_list": ["R_1", "R_2"],
            "note": "",
        }
        sp.check_receptor("X", good)
        no_key = {k: v for k, v in good.items() if k != "segment"}
        blank_listed = dict(good, segment="", segments_list=["", "R_2"])
        for bad in (
            dict(good, status="no_segment"),
            dict(good, segment=""),
            dict(good, segments_list=["R_2"]),
            no_key,
            blank_listed,
        ):
            with self.assertRaises(sp.MapMismatch):
                sp.check_receptor("X", bad)


class BuilderTests(unittest.TestCase):
    def test_every_pep_chain_whatever_its_type(self):
        rows = [
            {
                "PDB": "6ddf",
                "Name": "pep",
                "Type": "peptide",
                "Title": "DAMGO",
                "ChainID": "D",
            },
            {
                "PDB": "8K3Z",
                "Name": "pep",
                "Type": "protein",
                "Title": "CXCL12",
                "ChainID": "D",
            },
            {
                "PDB": "10TM",
                "Name": "pep",
                "Type": "peptide",
                "Title": "DAMGO",
                "ChainID": "H, S",
            },
            {
                "PDB": "10TM",
                "Name": "pep",
                "Type": "small-molecule",
                "Title": "a, b",
                "ChainID": "H",
            },
            {
                "PDB": "2RH1",
                "Name": "CAU",
                "Type": "small-molecule",
                "Title": "carazolol",
                "ChainID": "A",
            },
        ]
        out = builder.peptide_chains(rows)
        self.assertEqual(sorted(out), ["10TM", "6DDF", "8K3Z"])
        self.assertEqual(sorted(out["10TM"]), ["H", "S"])
        self.assertEqual(
            out["10TM"]["H"], ({"DAMGO", "a  b"}, {"peptide", "small-molecule"})
        )
        self.assertEqual(out["8K3Z"]["D"][1], {"protein"})

    def test_flat_keeps_one_field(self):
        self.assertEqual(builder._flat("a\tb\nc,d"), "a b c d")


if __name__ == "__main__":
    unittest.main()
