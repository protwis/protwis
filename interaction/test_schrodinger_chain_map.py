"""Unit tests for interaction.schrodinger_chain_map (no database).

Run with plain unittest after django.setup(), like test_schrodinger_import.
"""

import unittest

from interaction import schrodinger_chain_map as cm


SHA = "ab" * 32
INDEX_HEAD = (
    "# schema\tstructure-index/1\n# source\t2RH1.cif\n# cif_sha256\t" + SHA + "\n"
    "# atoms\t{n}\n"
    "label_asym\tauth_asym\tcomp\tauth_seq\ticode\tatom\tgroup\tx\ty\tz\n"
)


def index(*lines, head=INDEX_HEAD, n=None):
    return head.format(n=len(lines) if n is None else n) + "".join(lines)


def index_line(label, auth_asym, comp, auth_seq, icode, atom, group, x, y, z):
    return (
        "\t".join(
            [
                label,
                auth_asym,
                comp,
                auth_seq,
                icode,
                atom,
                group,
                "%.3f" % x,
                "%.3f" % y,
                "%.3f" % z,
            ]
        )
        + "\n"
    )


def pdb_line(group, n, atom, resname, chain, resnum, x, y, z, element):
    return "%-6s%5d %-4s %3s %1s%4s    %8.3f%8.3f%8.3f  1.00 20.00          %2s\n" % (
        group,
        n,
        atom,
        resname[:3],
        chain,
        resnum,
        x,
        y,
        z,
        element,
    )


class ParseTests(unittest.TestCase):
    def test_index_gives_the_sha256_and_the_atoms_in_file_order(self):
        text = index(
            index_line("A", "AAA", "ASP", "113", "", "CA", "ATOM", 1, 2, 3),
            index_line("B", "AAA", "Q6Q", "1000", "A", "C1", "HETATM", 7, 8, 9),
        )
        sha, atoms = cm.parse_structure_index(text)
        self.assertEqual(sha, SHA)
        self.assertEqual(
            [
                (
                    a["comp"],
                    a["auth_asym"],
                    a["auth_seq"],
                    a["label_asym"],
                    a["icode"],
                    a["atom"],
                    a["group"],
                )
                for a in atoms
            ],
            [
                ("ASP", "AAA", "113", "A", "", "CA", "ATOM"),
                ("Q6Q", "AAA", "1000", "B", "A", "C1", "HETATM"),
            ],
        )
        self.assertEqual(atoms[1]["key"], "7.000 8.000 9.000")
        self.assertEqual(atoms[0]["key"], cm.coord_key(1, 2, 3))

    def test_index_keeps_a_quote_inside_an_atom_name(self):
        _sha, atoms = cm.parse_structure_index(
            index(index_line("C", "A", "NAG", "901", "", "C1'", "HETATM", 1, 2, 3))
        )
        self.assertEqual(atoms[0]["atom"], "C1'")

    def test_index_of_another_schema_is_refused(self):
        with self.assertRaisesRegex(cm.ParseError, "schema"):
            cm.parse_structure_index(
                index(
                    index_line("A", "A", "ALA", "1", "", "CA", "ATOM", 1, 2, 3),
                    head=INDEX_HEAD.replace("structure-index/1", "structure-index/2"),
                )
            )

    def test_index_without_a_sha256_is_refused(self):
        for bad in ("", "abc", "AB" * 32):
            with self.assertRaisesRegex(cm.ParseError, "cif_sha256"):
                cm.parse_structure_index(
                    index(
                        index_line("A", "A", "ALA", "1", "", "CA", "ATOM", 1, 2, 3),
                        head=INDEX_HEAD.replace(SHA, bad),
                    )
                )

    def test_index_with_other_columns_is_refused(self):
        with self.assertRaisesRegex(cm.ParseError, "columns"):
            cm.parse_structure_index(
                index(
                    index_line("A", "A", "ALA", "1", "", "CA", "ATOM", 1, 2, 3),
                    head=INDEX_HEAD.replace("\tatom\tgroup", "\tgroup\tatom"),
                )
            )

    def test_index_row_shape_and_content_are_checked(self):
        with self.assertRaisesRegex(cm.ParseError, "fields"):
            cm.parse_structure_index(index("A\tA\tALA\t1\t\tCA\tATOM\t1.000\t2.000\n"))
        with self.assertRaisesRegex(cm.ParseError, "coordinate"):
            cm.parse_structure_index(
                index("A\tA\tALA\t1\t\tCA\tATOM\tx\t2.000\t3.000\n")
            )
        with self.assertRaisesRegex(cm.ParseError, "no atom"):
            cm.parse_structure_index(index())

    def test_index_coordinates_must_have_three_decimals(self):
        for bad in ("1.00", "1.0000", "1", "nan", "inf", "1e3", "+1.000"):
            for line in (
                "A\tA\tALA\t1\t\tCA\tATOM\t%s\t2.000\t3.000\n",
                "A\tA\tALA\t1\t\tCA\tATOM\t1.000\t%s\t3.000\n",
                "A\tA\tALA\t1\t\tCA\tATOM\t1.000\t2.000\t%s\n",
            ):
                with self.assertRaisesRegex(cm.ParseError, "coordinate"):
                    cm.parse_structure_index(index(line % bad))
        _sha, atoms = cm.parse_structure_index(
            index("A\tA\tALA\t1\t\tCA\tATOM\t-0.000\t2.000\t-31.500\n")
        )
        self.assertEqual(atoms[0]["key"], cm.coord_key("-0.000", 2, -31.5))

    def test_index_record_type_is_atom_or_hetatm(self):
        with self.assertRaisesRegex(cm.ParseError, "record type"):
            cm.parse_structure_index(
                index(index_line("A", "A", "ALA", "1", "", "CA", "atom", 1, 2, 3))
            )

    def test_index_atom_count_catches_a_truncated_file(self):
        line = index_line("A", "A", "ALA", "1", "", "CA", "ATOM", 1, 2, 3)
        with self.assertRaisesRegex(cm.ParseError, "truncated"):
            cm.parse_structure_index(index(line, n=2))
        with self.assertRaisesRegex(cm.ParseError, "atom count"):
            cm.parse_structure_index(
                index(line, head=INDEX_HEAD.replace("# atoms\t{n}\n", ""))
            )
        with self.assertRaisesRegex(cm.ParseError, "atom count"):
            cm.parse_structure_index(
                index(line, head=INDEX_HEAD.replace("# atoms\t{n}", "# atoms\t\u0661"))
            )

    def test_index_lives_beside_the_products(self):
        self.assertEqual(
            cm.index_path("/e1", "2RH1"), "/e1/2RH1/2RH1_structure_index.tsv"
        )

    def test_gpcrdb_text_reads_like_build_structures(self):
        text = (
            pdb_line("ATOM", 1, "CA", "ASP", "A", 113, 1, 2, 3, "C")
            + pdb_line("HETATM", 2, "C1", "A1C5S", "A", 1202, 7, 8, 9, "C")
            + pdb_line("HETATM", 3, "H1", "A1C5S", "A", 1202, 7, 8, 10, "H")
            + "ENDMDL\n"
            + pdb_line("ATOM", 4, "CA", "ASP", "A", 113, 9, 9, 9, "C")
        )
        atoms = cm.parse_gpcrdb_pdb(text)
        self.assertEqual(
            [(a["chain"], a["resnum"], a["resname"]) for a in atoms],
            [("A", "113", "ASP"), ("A", "1202", "A1C")],
        )


def lig(comp, label, auth_asym, auth_seq, coords):
    return [
        {
            "label_asym": label,
            "auth_asym": auth_asym,
            "comp": comp,
            "auth_seq": auth_seq,
            "icode": "",
            "atom": "C%d" % i,
            "group": "HETATM",
            "key": cm.coord_key(*xyz),
        }
        for i, xyz in enumerate(coords)
    ]


def glig(resname, chain, resnum, coords):
    return [
        {
            "chain": chain,
            "resnum": resnum,
            "icode": "",
            "resname": resname[:3],
            "atom": "C%d" % i,
            "group": "HETATM",
            "key": cm.coord_key(*xyz),
        }
        for i, xyz in enumerate(coords)
    ]


A_XYZ = [(1, 1, 1), (2, 2, 2), (3, 3, 3)]
B_XYZ = [(40, 1, 1), (41, 2, 2), (42, 3, 3)]


class ResolveAnchorTests(unittest.TestCase):
    def test_renamed_chain_resolved_by_coordinates_and_label(self):
        # 6ZIN: GPCRdb A:1000, product Q6Q_AAA_1000, annotation label B.
        cif = lig("Q6Q", "B", "AAA", "1000", A_XYZ)
        g = glig("Q6Q", "A", "1000", A_XYZ)
        row = cm.resolve_anchor("6ZIN", "Q6Q", "A:1000", cif, g, ["Q6Q_AAA_1000"], "B")
        self.assertEqual(
            (row["instance"], row["status"], row["source"], row["n_exact"]),
            ("Q6Q_AAA_1000", "ok", "coord+label", 3),
        )

    def test_swapped_annotation_is_errata_and_coordinates_win(self):
        # 6N51: annotation pairs A:903 with label G, but G is the B copy.
        cif = lig("QUS", "J", "A", "903", A_XYZ) + lig("QUS", "G", "B", "903", B_XYZ)
        g = glig("QUS", "A", "903", A_XYZ) + glig("QUS", "B", "903", B_XYZ)
        row = cm.resolve_anchor(
            "6N51", "QUS", "A:903", cif, g, ["QUS_A_903", "QUS_B_903"], "G"
        )
        self.assertEqual(
            (row["instance"], row["status"], row["source"]),
            ("QUS_A_903", "errata", "coord_exact"),
        )
        self.assertEqual(row["label_instance"], "QUS_B_903")

    def test_no_label_is_coord_only(self):
        cif = lig("RRY", "D", "C", "1", A_XYZ)
        row = cm.resolve_anchor(
            "9UTB", "RRY", "C:1", cif, glig("RRY", "C", "1", A_XYZ), ["RRY_C_1"], None
        )
        self.assertEqual((row["status"], row["source"]), ("ok", "coord_only"))

    def test_drifted_model_needs_label_agreement(self):
        cif = lig("2CU", "E", "R", "502", A_XYZ)
        g = glig("2CU", "R", "502", [(1.5, 1, 1), (2.5, 2, 2), (3.5, 3, 3)])
        ok = cm.resolve_anchor("7V68", "2CU", "R:502", cif, g, ["2CU_R_502"], "E")
        self.assertEqual(
            (ok["instance"], ok["status"], ok["source"]),
            ("2CU_R_502", "ok", "fallback+label"),
        )
        nolabel = cm.resolve_anchor("7V68", "2CU", "R:502", cif, g, ["2CU_R_502"], None)
        self.assertEqual((nolabel["instance"], nolabel["status"]), ("", "unresolved"))
        wrong = cm.resolve_anchor("7V68", "2CU", "R:502", cif, g, ["2CU_R_502"], "Z")
        self.assertEqual(wrong["status"], "unresolved")

    def test_no_product(self):
        row = cm.resolve_anchor("8E0G", "A1A7R", "A:54", [], [], ["CLR_A_403"], None)
        self.assertEqual(row["status"], "no_product")

    def test_five_char_het_matches_truncated_gpcrdb_name(self):
        cif = lig("A1C5S", "C", "A", "1202", A_XYZ)
        g = glig("A1C5S", "A", "1202", A_XYZ)
        row = cm.resolve_anchor(
            "10KT", "A1C5S", "A:1202", cif, g, ["A1C5S_A_1202"], "C"
        )
        # Decided by coordinates, not by the name fallback.
        self.assertEqual(
            (row["instance"], row["status"], row["source"], row["n_exact"]),
            ("A1C5S_A_1202", "ok", "coord+label", 3),
        )

    def test_identical_coordinates_on_two_instances_is_unresolved(self):
        cif = lig("XYZ", "C", "A", "1", A_XYZ) + lig("XYZ", "D", "B", "1", A_XYZ)
        row = cm.resolve_anchor(
            "1ABC",
            "XYZ",
            "A:1",
            cif,
            glig("XYZ", "A", "1", A_XYZ),
            ["XYZ_A_1", "XYZ_B_1"],
            "C",
        )
        self.assertEqual(row["status"], "unresolved")

    def test_split_tokens(self):
        self.assertEqual(cm.split_tokens("R:401, R:402"), ["R:401", "R:402"])
        self.assertEqual(cm.split_tokens("A:12B"), ["A:12B"])
        for value in ("", None, "L", "AAA:1", "A:1, L"):
            self.assertEqual(cm.split_tokens(value), [], value)

    def test_annotation_labels_copy_aligned(self):
        rows = [
            {
                "PDB": "6N51",
                "Name": "QUS",
                "Residue_seq_id": "A:903, B:903",
                "label_asym_id": "G, J",
            },
            {
                "PDB": "1ABC",
                "Name": "XYZ",
                "Residue_seq_id": "A:1, A:2",
                "label_asym_id": "C",
            },
        ]
        labels = cm.annotation_labels(rows)
        self.assertEqual(labels[("6N51", "QUS", "A:903")], "G")
        self.assertEqual(labels[("6N51", "QUS", "B:903")], "J")
        self.assertNotIn(("1ABC", "XYZ", "A:1"), labels)


def ca(chain, resnum, xyz, name="ASP", gpcrdb=False):
    if gpcrdb:
        return {
            "chain": chain,
            "resnum": resnum,
            "icode": "",
            "resname": name,
            "atom": "CA",
            "group": "ATOM",
            "key": cm.coord_key(*xyz),
        }
    return {
        "label_asym": "A",
        "auth_asym": chain,
        "comp": name,
        "auth_seq": resnum,
        "icode": "",
        "atom": "CA",
        "group": "ATOM",
        "key": cm.coord_key(*xyz),
    }


class ResolveReceptorTests(unittest.TestCase):
    def test_renamed_receptor_chain(self):
        cif = [ca("AAA", "113", (1, 1, 1)), ca("AAA", "114", (2, 2, 2))]
        g = [
            ca("A", "113", (1, 1, 1), gpcrdb=True),
            ca("A", "114", (2, 2, 2), gpcrdb=True),
        ]
        row = cm.resolve_receptor("6ZIN", "A", cif, g)
        self.assertEqual(
            (row["auth_chain"], row["status"], row["method"], row["n_ca_matched"]),
            ("AAA", "ok", "exact", 2),
        )

    def test_renumbered_is_refused(self):
        cif = [ca("R", "10", (1, 1, 1))]
        g = [ca("R", "11", (1, 1, 1), gpcrdb=True)]
        self.assertEqual(cm.resolve_receptor("X", "R", cif, g)["status"], "renumbered")

    def test_drift_identity_requires_every_residue(self):
        cif = [ca("R", "10", (1, 1, 1)), ca("R", "11", (2, 2, 2), name="GLY")]
        g_ok = [ca("R", "10", (1.5, 1, 1), gpcrdb=True)]
        row = cm.resolve_receptor("7V68", "R", cif, g_ok)
        self.assertEqual(
            (row["auth_chain"], row["status"], row["method"]),
            ("R", "ok", "identity_drift"),
        )
        g_bad = [ca("R", "12", (9, 9, 9), gpcrdb=True)]
        self.assertEqual(
            cm.resolve_receptor("7V68", "R", cif, g_bad)["status"], "unresolved"
        )

    def test_first_of_comma_preferred_chain(self):
        cif = [ca("A", "1", (1, 1, 1))]
        g = [ca("A", "1", (1, 1, 1), gpcrdb=True)]
        self.assertEqual(
            cm.resolve_receptor("X", "A,B", cif, g)["preferred_chain"], "A"
        )


if __name__ == "__main__":
    unittest.main()
