r"""
Unit tests for interaction.stored_interactions (no database).

    python -c "import django; django.setup(); import unittest; \
        unittest.main(module='interaction.test_stored_interactions', argv=['x'])"
"""

import unittest

from interaction import stored_interactions as si
from interaction.views import regexaa


class BuildResultsTests(unittest.TestCase):
    ROWS = [
        (
            ("ZMA", False),
            "D",
            113,
            "polar_donor_protein",
            "polar (hydrogen bond)",
            "polar",
            "protein",
        ),
        (("ZMA", False), "F", 290, "hyd", "hydrophobic", "hydrophobic", ""),
        (("CLR", False), "W", 286, "hyd", "hydrophobic", "hydrophobic", None),
        (("ZMA", False), "X", 400, "hyd", "hydrophobic", "hydrophobic", ""),
        (("ZMA", False), "N", 253, "acc", "accessible", "hidden", ""),
        (("DAMGO", True), "D", 147, "polar_neg_protein", "polar", "polar", "protein"),
        (("DAMGO", True), "Y", 148, "hyd", "hydrophobic", "hydrophobic", ""),
        (("DAMGO", True), "W", 293, "hyd", "hydrophobic", "hydrophobic", ""),
    ]

    def test_shape_order_and_residue_names(self):
        out = si.build_results(self.ROWS, "A")
        # HET ligands first even when a chain has more rows; hidden rows not counted.
        self.assertEqual(list(out), ["ZMA", "CLR", "DAMGO"])
        self.assertEqual(
            (out["ZMA"]["score"], out["CLR"]["score"], out["DAMGO"]["score"]), (2, 1, 3)
        )
        self.assertEqual(len(out["ZMA"]["interactions"]), 3)
        self.assertEqual(
            out["ZMA"]["interactions"][0],
            [
                "ASP113A",
                "",
                "polar_donor_protein",
                "polar (hydrogen bond)",
                "polar",
                "protein",
            ],
        )
        self.assertEqual(
            out["CLR"]["interactions"],
            [["TRP286A", "", "hyd", "hydrophobic", "hydrophobic", ""]],
        )
        # The page splits the residue with regexaa.
        self.assertEqual(regexaa(out["ZMA"]["interactions"][0][0]), ("D", "113", "A"))

    def test_a_chain_only_structure_has_the_chain_as_main(self):
        out = si.build_results([r for r in self.ROWS if r[0][1]], "R")
        self.assertEqual(list(out), ["DAMGO"])

    def test_ties_keep_key_order_and_empty_is_empty(self):
        out = si.build_results(
            [
                (("B", False), "A", 1, "hyd", "h", "hydrophobic", ""),
                (("A", False), "A", 2, "hyd", "h", "hydrophobic", ""),
            ],
            "R",
        )
        self.assertEqual(list(out), ["A", "B"])
        self.assertEqual(si.build_results([], "R"), {})

    def test_two_copies_of_one_het_stay_apart(self):
        keys = si.anchor_keys(
            {
                (1, "CY8", "x", "A:1201"),
                (2, "CY8", "x", "A:1202"),
                (3, "ZMA", "y", "A:401"),
                (4, "pep", "DAMGO", "B"),
            }
        )
        self.assertEqual(
            keys,
            {
                1: ("CY8 A:1201", False),
                2: ("CY8 A:1202", False),
                3: ("ZMA", False),
                4: ("DAMGO", True),
            },
        )

    def test_ligand_key(self):
        self.assertEqual(si.ligand_key("zma", "ZM241385"), ("ZMA", False))
        self.assertEqual(si.ligand_key(" pep ", "DAMGO"), ("DAMGO", True))


if __name__ == "__main__":
    unittest.main()
