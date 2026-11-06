"""Unit tests for the database-free layer of interaction.schrodinger_import."""

import os
import shutil
import tempfile
import types
import unittest

from interaction import schrodinger_chain_map as cm
from interaction import schrodinger_import as si


def row(
    family, direction="", seq=100, aa="F", chain="A", atom="CB", block="ATOM", lig=""
):
    return {
        "feature_family": family,
        "direction": direction,
        "receptor_atom_name": atom,
        "receptor_pdb_block": block,
        "ligand_pdb_block": lig,
        "receptor_residue": {
            "name_1_letter": aa,
            "pdb_residue_number": seq,
            "chain_id": chain,
            "insertion_code": "",
        },
    }


ANCHOR_MAP = {
    ("6ZIN", "Q6Q", "A:1000"): {"status": "ok", "instance": "Q6Q_AAA_1000", "note": ""},
    ("6N51", "QUS", "A:903"): {
        "status": "errata",
        "instance": "QUS_A_903",
        "note": "label says B",
    },
    ("8E0G", "A1A7R", "A:54"): {
        "status": "no_product",
        "instance": "",
        "note": "no instance",
    },
    ("7V68", "2CU", "R:502"): {"status": "unresolved", "instance": "", "note": "drift"},
    ("9X9X", "U7D", "R:601"): {"status": "ok", "instance": "U7D_R_601", "note": ""},
    ("9X9X", "U7D", "R:602"): {"status": "no_product", "instance": "", "note": "gone"},
    ("7E2X", "CLR", ""): {
        "status": "all_copies",
        "instance": "CLR_A_1;CLR_A_2",
        "note": "",
    },
    ("6CMO", "RET", ""): {
        "status": "no_product",
        "instance": "",
        "note": "not modelled",
    },
}


class AnchorInstancesTests(unittest.TestCase):
    def test_renamed_chain(self):
        self.assertEqual(
            si.anchor_instances("6zin", "q6q", "A:1000", ANCHOR_MAP, []),
            (["Q6Q_AAA_1000"], "mapped", []),
        )

    def test_errata_still_imports_and_is_reported(self):
        names, mode, notes = si.anchor_instances("6N51", "QUS", "A:903", ANCHOR_MAP, [])
        self.assertEqual((names, mode), (["QUS_A_903"], "mapped"))
        self.assertTrue(notes and notes[0].startswith("errata: "))

    def test_no_product_and_partial(self):
        self.assertEqual(
            si.anchor_instances("8E0G", "A1A7R", "A:54", ANCHOR_MAP, [])[:2],
            ([], "no_product"),
        )
        self.assertEqual(
            si.anchor_instances("9X9X", "U7D", "R:601, R:602", ANCHOR_MAP, [])[:2],
            (["U7D_R_601"], "mapped_partial"),
        )

    def test_chain_res_without_residue_uses_all_copies_row(self):
        self.assertEqual(
            si.anchor_instances("7E2X", "CLR", None, ANCHOR_MAP, [])[:2],
            (["CLR_A_1", "CLR_A_2"], "all_copies"),
        )
        self.assertEqual(
            si.anchor_instances("6CMO", "RET", "", ANCHOR_MAP, [])[:2],
            ([], "no_product"),
        )

    def test_unresolved_and_missing_are_loud(self):
        with self.assertRaises(si.UnresolvedAnchor):
            si.anchor_instances("7V68", "2CU", "R:502", ANCHOR_MAP, [])
        with self.assertRaises(si.MapMismatch):
            si.anchor_instances("2RH1", "CAU", "A:408", ANCHOR_MAP, [])

    def test_receptor_chain(self):
        rmap = {
            "6ZIN": {"status": "ok", "auth_chain": "AAA", "note": ""},
            "7V68": {"status": "unresolved", "auth_chain": "", "note": "x"},
        }
        self.assertEqual(si.receptor_chain("6zin", rmap), "AAA")
        with self.assertRaises(si.UnresolvedAnchor):
            si.receptor_chain("7V68", rmap)
        with self.assertRaises(si.MapMismatch):
            si.receptor_chain("2RH1", rmap)


class MapGuardTests(unittest.TestCase):
    def test_no_product_with_a_copy_in_the_tree_is_a_map_mismatch(self):
        with self.assertRaises(si.MapMismatch):
            si.anchor_instances("8E0G", "A1A7R", "A:54", ANCHOR_MAP, ["A1A7R_A_54"])
        with self.assertRaises(si.MapMismatch):
            si.anchor_instances("6CMO", "RET", "", ANCHOR_MAP, ["RET_A_1"])
        # a copy of another HET does not count
        self.assertEqual(
            si.anchor_instances("8E0G", "A1A7R", "A:54", ANCHOR_MAP, ["CLR_A_403"])[:2],
            ([], "no_product"),
        )

    @staticmethod
    def sli(ref, chain_res):
        return types.SimpleNamespace(pdb_reference=ref, chain_res=chain_res)

    def test_map_must_cover_every_database_copy(self):
        amap = {
            ("9X9X", "U7D", "R:601"): {},
            ("9X9X", "U7D", "R:602"): {},
            ("9X9X", "CLR", ""): {},
        }
        si.check_map_covers(
            "9X9X", [self.sli("U7D", "R:601, R:602"), self.sli("CLR", None)], amap
        )
        with self.assertRaises(
            si.MapMismatch
        ):  # database has a copy the map does not list
            si.check_map_covers(
                "9X9X",
                [self.sli("U7D", "R:601, R:602, R:603"), self.sli("CLR", "")],
                amap,
            )

    def test_map_may_list_copies_the_database_cannot_hold(self):
        """The annotation splits copies by chain; SLI has no copy dimension."""
        amap = {
            ("9X9X", "U7D", "R:601"): {},
            ("9X9X", "U7D", "S:601"): {},
            ("9X9X", "CLR", ""): {},
        }
        self.assertEqual(
            si.check_map_covers(
                "9X9X", [self.sli("U7D", "R:601"), self.sli("CLR", "")], amap
            ),
            [("U7D", "S:601")],
        )
        # a whole HET the database does not have is also extra, not an error,
        # but it is still handed back to be reported
        self.assertEqual(
            si.check_map_covers("9X9X", [self.sli("U7D", "R:601")], amap),
            [("CLR", ""), ("U7D", "S:601")],
        )

    def test_map_covers_is_case_insensitive_in_the_pdb_code(self):
        amap = {("9X9X", "U7D", "R:601"): {}}
        self.assertEqual(
            si.check_map_covers("9x9x", [self.sli("U7D", "R:601")], amap), []
        )

    def test_map_rows_of_another_structure_do_not_count_as_coverage(self):
        amap = {("OTHR", "U7D", "R:601"): {}}
        with self.assertRaises(si.MapMismatch):
            si.check_map_covers("9X9X", [self.sli("U7D", "R:601")], amap)

    def _write(self, text):
        fh = tempfile.NamedTemporaryFile("w", suffix=".tsv", delete=False)
        fh.write(text)
        fh.close()
        self.addCleanup(os.unlink, fh.name)
        return fh.name

    def test_a_quoted_newline_in_a_row_cannot_forge_a_header(self):
        """A note may hold a newline; its continuation must stay in the body."""
        note = "unreadable:\n# schema\tengine1-chainmap/999"
        path = self._write(
            "# schema\tengine1-chainmap/1\n"
            "pdb\thet\ttoken\tnote\n"
            '6ZIN\tQ6Q\tA:1\t"%s"\n' % note
        )
        head, _, table = si._read_map(path)
        self.assertEqual(head, {"schema": "engine1-chainmap/1"})
        self.assertEqual(len(table), 1)
        self.assertEqual(table[0]["note"], note)

    def test_fingerprints(self):
        rmap = {
            "6ZIN": {
                "gpcrdb_text_sha256": cm.text_sha256("ATOM 1\n"),
                "product_instances_sha256": cm.instances_sha256(
                    ["Q6Q_AAA_1000", "CLR_AAA_1"]
                ),
            }
        }
        si.check_fingerprints("6zin", rmap, "ATOM 1\n", ["CLR_AAA_1", "Q6Q_AAA_1000"])
        with self.assertRaises(si.MapMismatch):
            si.check_fingerprints(
                "6ZIN", rmap, "ATOM 2\n", ["CLR_AAA_1", "Q6Q_AAA_1000"]
            )
        with self.assertRaises(si.MapMismatch):
            si.check_fingerprints("6ZIN", rmap, "ATOM 1\n", ["Q6Q_AAA_1000"])
        with self.assertRaises(si.MapMismatch):
            si.check_fingerprints("2RH1", rmap, "", [])

    def test_ok_receptor_row_without_chain_is_refused(self):
        with self.assertRaises(si.MapMismatch):
            si.receptor_chain(
                "X", {"X": {"status": "ok", "auth_chain": "", "note": ""}}
            )

    def test_a_stale_map_is_reported_before_its_receptor_row(self):
        old, new, names = "ATOM 1\n", "ATOM 2\n", ["Q6Q_AAA_1000"]
        row = {
            "auth_chain": "AAA",
            "status": "ok",
            "note": "",
            "gpcrdb_text_sha256": cm.text_sha256(old),
            "product_instances_sha256": cm.instances_sha256(names),
        }
        chain = si.checked_receptor_chain("6zin", {"6ZIN": row}, old, names)
        self.assertEqual(chain, "AAA")
        renumbered = {"6ZIN": dict(row, status="renumbered", note="numbering differs")}
        with self.assertRaisesRegex(si.MapMismatch, "rebuild"):
            si.checked_receptor_chain("6ZIN", renumbered, new, names)
        with self.assertRaisesRegex(si.UnresolvedAnchor, "numbering differs"):
            si.checked_receptor_chain("6ZIN", renumbered, old, names)
        gave_up = dict(row, status="unresolved", note="input unreadable")
        gave_up.update(auth_chain="", gpcrdb_text_sha256="")
        with self.assertRaisesRegex(si.UnresolvedAnchor, "input unreadable"):
            si.checked_receptor_chain("6ZIN", {"6ZIN": gave_up}, new, names)
        no_text = {"6ZIN": dict(row, gpcrdb_text_sha256="")}
        with self.assertRaisesRegex(si.MapMismatch, "rebuild"):
            si.checked_receptor_chain("6ZIN", no_text, old, names)
        with self.assertRaisesRegex(si.MapMismatch, "product instances"):
            si.checked_receptor_chain("6ZIN", renumbered, old, ["CLR_AAA_1"])


CHAINMAP_HEADER = "".join(
    "# {}\t{}\n".format(k, v)
    for k, v in [
        ("schema", "engine1-chainmap/1"),
        ("pdb", "6ZIN"),
        ("annotation_commit", "9fe1875"),
        ("ligands_sha256", "aa"),
        ("structures_sha256", "bb"),
        ("gpcrdb_pdb_sha256", "cc"),
        ("cif_sha256", "dd"),
        ("builder_sha256", "ee"),
        ("product_summary", "yes"),
        ("receptor.preferred_chain", "A"),
        ("receptor.auth_chain", "A"),
        ("receptor.status", "ok"),
        ("receptor.method", "exact"),
        ("receptor.n_ca_gpcrdb", "10"),
        ("receptor.n_ca_matched", "10"),
        ("receptor.renumbered", "0"),
        ("receptor.note", ""),
        ("receptor.gpcrdb_text_sha256", "ff"),
        ("receptor.product_instances_sha256", "gg"),
    ]
)
CHAINMAP_COLUMNS = "\t".join(cm.ANCHOR_COLUMNS) + "\n"


def chainmap_row(  # nosec B107
    pdb="6ZIN", het="q6q", token="A:1000", instance="Q6Q_A_1000", status="ok"
):
    row = dict.fromkeys(cm.ANCHOR_COLUMNS, "")
    row.update(pdb=pdb, het=het, token=token, instance=instance, status=status)
    return "\t".join(row[c] for c in cm.ANCHOR_COLUMNS) + "\n"


CHAINMAP_BODY = CHAINMAP_COLUMNS + chainmap_row()


class ChainmapDirTests(unittest.TestCase):
    """The per-PDB chain maps the build consumes without being told anything."""

    def setUp(self):
        self.root = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.root)

    def write(self, pdb, text):
        d = os.path.join(self.root, pdb)
        os.makedirs(d, exist_ok=True)
        with open(os.path.join(d, si.CHAINMAP_NAME), "w") as fh:
            fh.write(text)
        return os.path.join(d, si.CHAINMAP_NAME)

    def test_one_file_round_trips(self):
        path = self.write("6ZIN", CHAINMAP_HEADER + CHAINMAP_BODY)
        pdb, anchors, receptor, prov, ran = si.load_chainmap(path)
        self.assertEqual(pdb, "6ZIN")
        self.assertTrue(ran)
        self.assertEqual(anchors[("6ZIN", "Q6Q", "A:1000")]["instance"], "Q6Q_A_1000")
        self.assertEqual(receptor["auth_chain"], "A")
        self.assertEqual(receptor["product_instances_sha256"], "gg")
        self.assertEqual(prov["annotation_commit"], "9fe1875")

    def test_an_unknown_schema_is_refused(self):
        for schema in ("engine1-chainmap/2", "", "engine1-chainmap"):
            path = self.write(
                "6ZIN",
                CHAINMAP_HEADER.replace("engine1-chainmap/1", schema) + CHAINMAP_BODY,
            )
            with self.assertRaises(si.MapMismatch):
                si.load_chainmap(path)
        # a file with no header at all
        path = self.write("6ZIN", CHAINMAP_BODY)
        with self.assertRaises(si.MapMismatch):
            si.load_chainmap(path)

    def test_a_missing_receptor_field_is_refused(self):
        path = self.write(
            "6ZIN",
            CHAINMAP_HEADER.replace("# receptor.auth_chain\tA\n", "") + CHAINMAP_BODY,
        )
        with self.assertRaises(si.MapMismatch):
            si.load_chainmap(path)

    def test_a_row_of_another_structure_is_refused(self):
        path = self.write(
            "6ZIN",
            CHAINMAP_HEADER
            + CHAINMAP_BODY
            + chainmap_row(pdb="2RH1", het="CAU", token="A:408", instance="CAU_A_408"),
        )
        with self.assertRaises(si.MapMismatch):
            si.load_chainmap(path)

    def test_a_duplicate_anchor_key_is_refused(self):
        path = self.write("6ZIN", CHAINMAP_HEADER + CHAINMAP_BODY + chainmap_row())
        with self.assertRaises(si.MapMismatch):
            si.load_chainmap(path)

    def test_a_directory_without_a_map_is_named_not_skipped(self):
        self.write("6ZIN", CHAINMAP_HEADER + CHAINMAP_BODY)
        os.makedirs(os.path.join(self.root, "2RH1"))
        anchors, receptors, missing, _, not_run = si.load_chainmap_dir(
            self.root, ["6ZIN", "2RH1"]
        )
        self.assertEqual(missing, ["2RH1"])
        self.assertEqual(not_run, [])
        self.assertEqual(sorted(receptors), ["6ZIN"])
        self.assertEqual(len(anchors), 1)

    def test_the_tree_can_contradict_a_not_run_header(self):
        """A header from another run must not drop this run's products."""
        empty = CHAINMAP_HEADER.replace(
            "# product_summary\tyes", "# product_summary\tno"
        ).replace(
            "# receptor.product_instances_sha256\tgg",
            "# receptor.product_instances_sha256\t" + cm.instances_sha256([]),
        )
        self.write("6ZIN", empty + CHAINMAP_BODY)
        inst = os.path.join(self.root, "6ZIN", "Q6Q_A_1000")
        os.makedirs(inst)
        open(os.path.join(inst, "Q6Q_A_1000.yaml"), "w").close()
        _, receptors, _, _, not_run = si.load_chainmap_dir(self.root, ["6ZIN"])
        self.assertEqual(not_run, [])
        self.assertEqual(sorted(receptors), ["6ZIN"])

    def test_the_receptor_table_is_keyed_like_the_anchor_table(self):
        self.write("6zin", CHAINMAP_HEADER + CHAINMAP_BODY)
        anchors, receptors, _, _, _ = si.load_chainmap_dir(self.root, ["6zin"])
        self.assertEqual(sorted(receptors), ["6ZIN"])
        self.assertEqual(sorted({p for p, _, _ in anchors}), ["6ZIN"])

    def test_a_structure_nobody_ran_is_not_treated_as_one_with_no_ligand(self):
        """Clearing an anchor states that Engine 1 looked; it must have looked."""
        empty = CHAINMAP_HEADER.replace(
            "# product_summary\tyes", "# product_summary\tno"
        ).replace(
            "# receptor.product_instances_sha256\tgg",
            "# receptor.product_instances_sha256\t" + cm.instances_sha256([]),
        )
        self.write("6ZIN", empty + CHAINMAP_BODY)
        anchors, receptors, _, _, not_run = si.load_chainmap_dir(self.root, ["6ZIN"])
        self.assertEqual(not_run, ["6ZIN"])
        self.assertEqual(anchors, {})
        self.assertEqual(receptors, {})
        # the same structure WITH products is imported, summary or not
        both = CHAINMAP_HEADER.replace(
            "# product_summary\tyes", "# product_summary\tno"
        )
        self.write(
            "2RH1", both.replace("# pdb\t6ZIN", "# pdb\t2RH1") + CHAINMAP_COLUMNS
        )
        _, receptors, _, _, not_run = si.load_chainmap_dir(self.root, ["2RH1"])
        self.assertEqual(not_run, [])
        self.assertEqual(sorted(receptors), ["2RH1"])

    def test_a_missing_or_bogus_product_summary_is_refused(self):
        for value in ("", "maybe", "Yes"):
            path = self.write(
                "6ZIN",
                CHAINMAP_HEADER.replace(
                    "# product_summary\tyes", "# product_summary\t" + value
                )
                + CHAINMAP_BODY,
            )
            with self.assertRaises(si.MapMismatch):
                si.load_chainmap(path)
        path = self.write(
            "6ZIN",
            CHAINMAP_HEADER.replace("# product_summary\tyes\n", "") + CHAINMAP_BODY,
        )
        with self.assertRaises(si.MapMismatch):
            si.load_chainmap(path)

    def test_provenance_counts_a_merged_tree(self):
        self.write("6ZIN", CHAINMAP_HEADER + CHAINMAP_BODY)
        self.write(
            "2RH1",
            CHAINMAP_HEADER.replace("# pdb\t6ZIN", "# pdb\t2RH1").replace(
                "9fe1875", "deadbee"
            )
            + CHAINMAP_COLUMNS,
        )
        _, _, _, prov, _ = si.load_chainmap_dir(self.root, ["6ZIN", "2RH1"])
        self.assertEqual(prov["annotation_commit"], {"9fe1875": 1, "deadbee": 1})
        self.assertEqual(prov["ligands_sha256"], {"aa": 2})

    def test_a_file_in_the_wrong_directory_is_refused(self):
        d = os.path.join(self.root, "2RH1")
        os.makedirs(d)
        with open(os.path.join(d, si.CHAINMAP_NAME), "w") as fh:
            fh.write(CHAINMAP_HEADER + CHAINMAP_BODY)  # header says 6ZIN
        with self.assertRaises(si.MapMismatch):
            si.load_chainmap_dir(self.root, ["2RH1"])

    def test_product_pdb_codes_keeps_the_name_on_disk(self):
        """Upper-casing here would hide a lower-case delivery, or fold two into one."""
        self.write("6ZIN", CHAINMAP_HEADER + CHAINMAP_BODY)
        os.makedirs(os.path.join(self.root, "2rh1"))
        os.makedirs(os.path.join(self.root, ".rsync-partial"))
        with open(os.path.join(self.root, "_import_anomalies.csv"), "w") as fh:
            fh.write("x\n")
        self.assertEqual(si.product_pdb_codes(self.root), ["2rh1", "6ZIN"])

    def test_a_truncated_or_reshaped_file_is_refused(self):
        header_only = self.write("6ZIN", CHAINMAP_HEADER)
        with self.assertRaises(si.MapMismatch):
            si.load_chainmap(header_only)
        lost_a_column = self.write(
            "6ZIN", CHAINMAP_HEADER + CHAINMAP_BODY.replace("token\t", "", 1)
        )
        with self.assertRaises(si.MapMismatch):
            si.load_chainmap(lost_a_column)

    def test_crlf_does_not_leak_into_the_header_values(self):
        path = self.write(
            "6ZIN", (CHAINMAP_HEADER + CHAINMAP_BODY).replace("\n", "\r\n")
        )
        _, _, receptor, _, _ = si.load_chainmap(path)
        self.assertEqual(receptor["auth_chain"], "A")


class StructureVerdictTests(unittest.TestCase):
    """The pre-flight policy: what happens to a delivered structure, and why."""

    IMPORTABLE = dict(in_db=True, experimental=True, was_run=True, has_chainmap=True)

    def verdict(self, **kw):
        return si.structure_verdict(**dict(self.IMPORTABLE, **kw))

    def test_a_structure_we_serve_and_can_read_is_imported(self):
        self.assertIsNone(self.verdict())

    def test_a_delivered_directory_with_no_chainmap_fails_the_run(self):
        """Never a skip: its anchors would keep whatever legacy wrote."""
        self.assertEqual(
            self.verdict(has_chainmap=False), ("failed", "ERROR", "no_chainmap")
        )

    def test_a_structure_nobody_ran_is_left_alone(self):
        self.assertEqual(self.verdict(was_run=False), ("not_run", "INFO", "not_run"))

    def test_leaving_anchors_as_legacy_wrote_them_fails_like_any_other_mixing(self):
        """Same end state as a missing map, so the same answer."""
        self.assertEqual(
            self.verdict(was_run=False, anchors_at_risk=1),
            ("failed", "ERROR", "not_run"),
        )
        self.assertEqual(
            self.verdict(was_run=False, anchors_at_risk=1, allow_not_run=True),
            ("not_run", "WARNING", "not_run"),
        )
        # the escape hatch changes nothing when there is nothing at risk
        self.assertEqual(
            self.verdict(was_run=False, allow_not_run=True),
            ("not_run", "INFO", "not_run"),
        )

    def test_a_structure_this_import_does_not_serve_needs_no_chainmap(self):
        self.assertEqual(
            self.verdict(in_db=False, has_chainmap=False),
            ("structure_not_in_db", "WARNING", "structure_not_in_db"),
        )
        self.assertEqual(
            self.verdict(experimental=False, has_chainmap=False),
            ("not_experimental", "INFO", "not_experimental"),
        )

    def test_never_run_outranks_a_missing_chainmap(self):
        """Both are true of a directory holding nothing at all; say the useful one."""
        self.assertEqual(self.verdict(was_run=False, has_chainmap=False)[0], "not_run")


class StandardLigandLineTests(unittest.TestCase):
    # Lines as the producer writes them (no altloc column; widened for long names).
    CAU = (
        "HETATM    1  O17CAU A 408     -33.477  10.957   8.170  1.00 50.96           O"
    )
    FIVE = "HETATM 4817  C9 A1C5S A1202     -20.354 -13.384  28.940  1.00 22.97           C"
    AAA = "HETATM   10  N4 T0B AAA 601      -2.398  28.100  15.000  1.00 30.00           N"
    ABUT = "HETATM    7  C1 LIG A   1    -100.123-200.456-300.789  1.0033084.00           C"
    # A four-character atom name runs straight into the residue name, as the producer writes it.
    HNAME = (
        "HETATM   40 HN12CAU A 408     -30.000  11.000   9.000  1.00 50.96           H"
    )

    def check_standard(self, out, resname, chain, resnum, xyz):
        self.assertEqual(len(out), 78)
        self.assertEqual(out[16], " ")  # altloc column present and blank
        self.assertEqual(out[17:20].strip(), resname)
        self.assertEqual(out[21], chain)
        self.assertEqual(out[22:26].strip(), resnum)
        self.assertEqual((float(out[30:38]), float(out[38:46]), float(out[46:54])), xyz)

    def test_ordinary_line(self):
        out, capped = si.standard_ligand_line(self.CAU, "CAU", "A", "408", "", "A")
        self.check_standard(out, "CAU", "A", "408", (-33.477, 10.957, 8.17))
        self.assertEqual(out[12:16], " O17")
        self.assertFalse(capped)

    def test_five_char_code_is_cut_to_three(self):
        out, _ = si.standard_ligand_line(self.FIVE, "A1C5S", "A", "1202", "", "A")
        self.check_standard(out, "A1C", "A", "1202", (-20.354, -13.384, 28.94))

    def test_multichar_chain_becomes_gpcrdb_chain(self):
        out, _ = si.standard_ligand_line(self.AAA, "T0B", "AAA", "601", "", "A")
        self.check_standard(out, "T0B", "A", "601", (-2.398, 28.1, 15.0))

    def test_abutting_fields_and_capped_b_factor(self):
        out, capped = si.standard_ligand_line(self.ABUT, "LIG", "A", "1", "", "A")
        self.check_standard(out, "LIG", "A", "1", (-100.123, -200.456, -300.789))
        self.assertTrue(capped)
        self.assertEqual(out[60:66], "999.99")

    def test_four_char_hydrogen_name(self):
        out, _ = si.standard_ligand_line(self.HNAME, "CAU", "A", "408", "", "A")
        self.assertEqual(out[12:16], "HN12")
        self.assertEqual(out[76:78], " H")

    def test_mismatches_raise(self):
        for het, chain, resnum in (
            ("CAZ", "A", "408"),
            ("CAU", "B", "408"),
            ("CAU", "A", "409"),
        ):
            with self.assertRaises(si.MalformedLigandLine):
                si.standard_ligand_line(self.CAU, het, chain, resnum, "", "A")
        with self.assertRaises(si.MalformedLigandLine):
            si.standard_ligand_line("REMARK nothing here", "CAU", "A", "408", "", "A")
        with self.assertRaises(si.MalformedLigandLine):
            si.standard_ligand_line(self.CAU[:60], "CAU", "A", "408", "", "A")

    def test_block_and_instance_chains(self):
        text, _ = si.standard_ligand_block(
            self.CAU + "\n\n" + self.HNAME + "\n", "CAU_A_408", "A"
        )
        self.assertEqual(len(text.splitlines()), 2)
        amap = {
            ("6ZIN", "Q6Q", "A:1000"): {
                "status": "ok",
                "instance": "Q6Q_AAA_1000",
                "note": "",
            },
            ("7E2X", "CLR", ""): {
                "status": "all_copies",
                "instance": "CLR_R_602;CLR_R_603",
                "note": "",
            },
            ("9X9X", "LIG", ""): {
                "status": "all_copies",
                "instance": "LIG_AB_1",
                "note": "",
            },
        }
        self.assertEqual(
            si.instance_chains("6ZIN", "Q6Q", "A:1000", amap, ["Q6Q_AAA_1000"]),
            {"Q6Q_AAA_1000": "A"},
        )
        self.assertEqual(
            si.instance_chains("7E2X", "CLR", "", amap, ["CLR_R_602", "CLR_R_603"]),
            {"CLR_R_602": "R", "CLR_R_603": "R"},
        )
        with self.assertRaises(si.MapMismatch):
            si.instance_chains("9X9X", "LIG", "", amap, ["LIG_AB_1"])

    def test_capped_lines_are_returned(self):
        _, capped = si.standard_ligand_block(self.CAU, "CAU_A_408", "A")
        self.assertEqual(capped, [])
        _, capped = si.standard_ligand_block(self.ABUT, "LIG_A_1", "A")
        self.assertEqual(len(capped), 1)

    def test_fields_too_wide_raise(self):
        wide = (
            (
                "HETATM    1  C1 LIG A10000     -1.000   2.000   3.000  1.00 20.00           C",
                "10000",
                "",
            ),
            (
                "HETATM    1  C1 LIG A   1   -1000.500   2.000   3.000  1.00 20.00           C",
                "1",
                "",
            ),
            (
                "HETATM    1  C1 LIG A   1    10000.000   2.000   3.000  1.00 20.00           C",
                "1",
                "",
            ),
            (
                "HETATM    1  C1 LIG A   1      1.000   2.000   3.000  1.00-100.00           C",
                "1",
                "",
            ),
            (
                "HETATM    1  C1 LIG A   1      1.000   2.000   3.0001000.00 20.00           C",
                "1",
                "",
            ),
        )
        for line, resnum, icode in wide:
            with self.assertRaises(si.MalformedLigandLine, msg=line):
                si.standard_ligand_line(line, "LIG", "A", resnum, icode, "A")
        with self.assertRaises(si.MalformedLigandLine):
            si.standard_ligand_line(self.CAU, "CAU", "A", "408", "", "")


class CascadeGuardTests(unittest.TestCase):
    def test_only_expected_models_may_be_deleted(self):
        si._only_deleted({"structure.Fragment": 3}, {"structure.Fragment"})
        si._only_deleted(
            {"structure.Fragment": 3, "structure.Rotamer": 0}, {"structure.Fragment"}
        )
        with self.assertRaises(si.UnexpectedCascade):
            si._only_deleted(
                {"structure.PdbData": 1, "structure.Rotamer": 2}, {"structure.PdbData"}
            )

    def test_fragment_text(self):
        self.assertEqual(
            si.fragment_text(["HETATM 1", "HETATM 2"]), "HETATM 1\nHETATM 2\n"
        )
        self.assertEqual(si.fragment_text([]), "")


class RoutingTests(unittest.TestCase):
    # Every (family, direction) pair the producer emits.
    PRODUCTION_PAIRS = {
        ("VdW", ""): "Van der Waals",
        ("Accessible", ""): "acc",
        ("HPhob", ""): "hyd",
        ("Acceptor", "ligand-acceptor"): "polar_donor_protein",
        ("Donor", "ligand-donor"): "polar_acceptor_protein",
        ("Aromatic", "edge-to-face"): "aro_ef_protein",
        ("Aromatic", "face-to-face"): "aro_ff",
        ("NegCharge", "neg-pos"): "polar_double_pos_protein",
        ("PosCharge", "pos-neg"): "polar_double_neg_protein",
        ("Metal", ""): "metal_coordination_protein",
        ("PiCat", "ligand-cation"): "aro_ion_protein",
        ("PiCat", "receptor-cation"): "aro_ion_protein",
        ("Wat-HBond", ""): "water_bridge_protein",
        ("XBond", ""): "halogen_protein",
        ("Covalent", ""): "covalent",
    }

    def test_every_production_pair_routes(self):
        for (family, direction), slug in self.PRODUCTION_PAIRS.items():
            self.assertEqual(
                si.resolve_slug(family, direction), slug, (family, direction)
            )

    def test_the_covalent_slug_is_the_one_the_seed_migration_seeds_and_it_is_visible(
        self,
    ):
        # The map and the migration name the same slug, or the first Covalent
        # row finds no ResidueFragmentInteractionType. And the type must not be
        # "hidden": the pages leave hidden types out, which would import the
        # rows and show them to no one.
        import importlib

        mig = importlib.import_module(
            "interaction.migrations.0009_seed_schrodinger_interaction_types"
        )
        slug, name, type_, direction = next(
            row for row in mig.SEEDED_TYPES if row[0] == "covalent"
        )
        self.assertEqual(si.resolve_slug("Covalent", ""), slug)
        self.assertNotEqual(type_, "hidden")
        self.assertEqual(
            (slug, name, type_, direction),
            ("covalent", "covalent bond", "covalent", ""),
        )

    def test_the_seed_writes_a_visible_covalent_type(self):
        # What seed() hands the ORM, not only the constant.
        import importlib

        mig = importlib.import_module(
            "interaction.migrations.0009_seed_schrodinger_interaction_types"
        )
        seen = []

        class _Objects:
            def get_or_create(self, **kwargs):
                seen.append(kwargs)
                return object(), True

        class _Model:
            objects = _Objects()

        class _Apps:
            def get_model(self, app, name):
                assert (app, name) == ("interaction", "ResidueFragmentInteractionType")
                return _Model

        mig.seed(_Apps(), None)
        self.assertEqual(
            [k for k in seen if k["slug"] == "covalent"],
            [
                {
                    "slug": "covalent",
                    "defaults": {
                        "name": "covalent bond",
                        "type": "covalent",
                        "direction": "",
                    },
                }
            ],
        )

    def test_the_map_is_the_one_the_producer_pins(self):
        # The producer pins the same digest (see the header of the map).
        import hashlib

        path = os.path.join(os.path.dirname(si.__file__), "interaction_type_map.yaml")
        with open(path, "rb") as fh:
            digest = hashlib.sha256(fh.read()).hexdigest()
        self.assertEqual(
            digest, "7b3fb89dbac84b7463de7fe59f902673046668d3f54dd5b39a81dd1c49d277e8"
        )

    def test_none_direction_is_empty(self):
        self.assertEqual(si.resolve_slug("HPhob", None), "hyd")

    def test_unknown_pair_raises(self):
        with self.assertRaises(si.UnroutableRow):
            si.resolve_slug("Acceptor", "ligand-donor")

    def test_backbone_override(self):
        self.assertEqual(
            si.apply_backbone_override("polar_donor_protein", "N"), "polar_backbone"
        )
        self.assertEqual(
            si.apply_backbone_override("polar_acceptor_protein", " O "),
            "polar_backbone",
        )
        self.assertEqual(
            si.apply_backbone_override("polar_donor_protein", "OG"),
            "polar_donor_protein",
        )
        self.assertEqual(si.apply_backbone_override("hyd", "N"), "hyd")
        self.assertEqual(
            si.apply_backbone_override("polar_donor_protein", None),
            "polar_donor_protein",
        )

    def test_only_main_chain_n_and_o_promote(self):
        for atom in ("CA", "C", "CB", "OG", "OG1", "ND2", "NE2", "OH", "H", ""):
            for slug in ("polar_donor_protein", "polar_acceptor_protein"):
                self.assertEqual(
                    si.apply_backbone_override(slug, atom), slug, (slug, atom)
                )

    def test_required_slugs_exclude_water_bridge(self):
        self.assertEqual(
            si.required_slugs(),
            frozenset(
                {
                    "hyd",
                    "polar_donor_protein",
                    "polar_acceptor_protein",
                    "aro_ef_protein",
                    "aro_ff",
                    "polar_double_pos_protein",
                    "polar_double_neg_protein",
                    "metal_coordination_protein",
                    "aro_ion_protein",
                    "halogen_protein",
                    "polar_backbone",
                    "covalent",
                    "Van der Waals",
                    "acc",
                }
            ),
        )


class PlanRowsTests(unittest.TestCase):
    def test_accounting_identity_and_each_bucket(self):
        rows = [
            row("HPhob", seq=100),
            row("HPhob", seq=100),  # duplicate (100, hyd)
            row("Acceptor", "ligand-acceptor", seq=100, aa="F", atom="N"),  # backbone
            row("Wat-HBond", seq=101),  # excluded family
            row("HPhob", seq=102, aa="X"),  # non-standard residue
            row("HPhob", seq=103, chain="B"),  # other chain
            row("Metal", seq=104, aa="H"),
        ]
        records, counts, by_chain = si.plan_rows(rows, "A")
        self.assertEqual(counts["rows_in"], 7)
        self.assertEqual(counts["duplicate"], 1)
        self.assertEqual(counts["excluded_family"], 1)
        self.assertEqual(counts["nonstandard_residue"], 1)
        self.assertEqual(counts["other_chain"], 1)
        self.assertEqual(by_chain, {"B": 1})
        self.assertEqual(counts["planned"], len(records))
        self.assertEqual(
            counts["rows_in"],
            counts["excluded_family"]
            + counts["nonstandard_residue"]
            + counts["other_chain"]
            + counts["duplicate"]
            + counts["planned"],
        )
        self.assertEqual(
            [(r["sequence_number"], r["slug"]) for r in records],
            [
                (100, "hyd"),
                (100, "polar_backbone"),
                (104, "metal_coordination_protein"),
            ],
        )

    def test_distinct_slugs_on_one_residue_both_survive(self):
        rows = [
            row("PosCharge", "pos-neg", seq=113, aa="D"),
            row("Donor", "ligand-donor", seq=113, aa="D", atom="OD1"),
        ]
        records, counts, _ = si.plan_rows(rows, "A")
        self.assertEqual(
            {r["slug"] for r in records},
            {"polar_double_neg_protein", "polar_acceptor_protein"},
        )
        self.assertEqual(counts["duplicate"], 0)

    def test_empty_preferred_chain_keeps_every_chain(self):
        records, counts, _ = si.plan_rows([row("HPhob", chain="B")], "")
        self.assertEqual(len(records), 1)
        self.assertEqual(counts["other_chain"], 0)

    def test_unroutable_row_raises(self):
        with self.assertRaises(si.UnroutableRow):
            si.plan_rows([row("Bogus", "x")], "A")

    def test_receptor_chain_is_the_product_chain(self):
        # The receptor chain is the product (author) chain from the map.
        records, counts, _ = si.plan_rows([row("HPhob", chain="AAA")], "AAA")
        self.assertEqual((len(records), counts["other_chain"]), (1, 0))


class LigandLinesTests(unittest.TestCase):
    def test_collapsed_rows_merge_ligand_atoms_in_first_seen_order(self):
        rows = [
            row("HPhob", seq=100, lig="HETATM C1\nHETATM C2\n"),
            row("HPhob", seq=100, lig="HETATM C2\nHETATM C3\n"),
            row("Acceptor", "ligand-acceptor", seq=100, atom="OG", lig="HETATM O1\n"),
        ]
        records, counts, _ = si.plan_rows(rows, "A")
        self.assertEqual(counts["duplicate"], 1)
        self.assertEqual(
            [r["ligand_lines"] for r in records],
            [["HETATM C1", "HETATM C2", "HETATM C3"], ["HETATM O1"]],
        )


class ProductFilesTests(unittest.TestCase):
    def setUp(self):
        self.root = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.root)

    def _write(self, rel, text):
        path = os.path.join(self.root, rel)
        os.makedirs(os.path.dirname(path), exist_ok=True)
        with open(path, "w") as fh:
            fh.write(text)
        return path

    def test_instance_discovery(self):
        good = self._write(
            "2RH1/CAU_A_408/CAU_A_408.yaml", "result: {interactions: []}\n"
        )
        self._write("2RH1/summary.yaml", "x: 1\n")
        self._write("2RH1/not_an_instance/not_an_instance.yaml", "x: 1\n")
        os.makedirs(os.path.join(self.root, "2RH1", "CLR_A_1"))  # no YAML
        self.assertEqual(si.instance_yaml_paths(self.root, "2rh1"), {"CAU_A_408": good})
        self.assertEqual(si.instance_yaml_paths(self.root, "9ZZZ"), {})

    def test_read_rows(self):
        ok = self._write("a.yaml", "result:\n  interactions: []\n")
        self.assertEqual(si.read_instance_rows(ok), [])
        for text in (
            "result: {}\n",
            "- 1\n",
            "",
            "result: {interactions: 3}\n",
            "a: [\n",
        ):
            bad = self._write("b.yaml", text)
            with self.assertRaises(si.MalformedProduct, msg=text):
                si.read_instance_rows(bad)


class ScopeTests(unittest.TestCase):
    @staticmethod
    def sli(reference, ligand_type):
        return types.SimpleNamespace(
            pdb_reference=reference,
            ligand=types.SimpleNamespace(
                ligand_type=types.SimpleNamespace(slug=ligand_type)
            ),
        )

    def test_scope(self):
        self.assertTrue(si.is_in_scope(self.sli("CAU", "small-molecule")))
        self.assertTrue(si.is_in_scope(self.sli("CLR", "lipid")))
        self.assertFalse(si.is_in_scope(self.sli("pep", "small-molecule")))
        self.assertFalse(si.is_in_scope(self.sli("APO", "none")))
        # A component reference is Engine 1's whatever the database calls the
        # ligand (a peptide drug referenced by one component).
        self.assertTrue(si.is_in_scope(self.sli("D2U", "peptide")))
        self.assertTrue(si.is_in_scope(self.sli("XYZ", "protein")))
        self.assertFalse(si.is_in_scope(self.sli("PEP", "peptide")))
        self.assertFalse(si.is_in_scope(self.sli("PEP", "protein")))
        self.assertFalse(si.is_in_scope(self.sli(" apo ", "none")))
        self.assertFalse(si.is_in_scope(self.sli("", "small-molecule")))
        self.assertFalse(si.is_in_scope(self.sli(None, "small-molecule")))


class ProductContractTests(unittest.TestCase):
    def setUp(self):
        self.dir = tempfile.mkdtemp()
        os.makedirs(os.path.join(self.dir, "2RH1"))
        self.addCleanup(shutil.rmtree, self.dir)

    def summary(self, text):
        with open(os.path.join(self.dir, "2RH1", "summary.yaml"), "w") as fh:
            fh.write(text)

    def test_the_version_this_importer_reads_passes(self):
        self.assertEqual(si.PRODUCT_CONTRACT, "engine1/1.0")
        self.summary("contract_version: engine1/1.0\npdb_id: 2RH1\n")
        si.check_product_contract(self.dir, "2RH1")

    def test_another_version_no_version_no_summary_or_garbage_is_refused(self):
        for text in (
            "contract_version: engine1/2.0\npdb_id: 2RH1\n",
            "pdb_id: 2RH1\n",
            "- a list\n",
            "contract_version: [unclosed\n",
        ):
            self.summary(text)
            with self.assertRaises(si.MalformedProduct):
                si.check_product_contract(self.dir, "2RH1")
        with open(os.path.join(self.dir, "2RH1", "summary.yaml"), "wb") as fh:
            fh.write(b"contract_version: \xff\n")
        with self.assertRaisesRegex(si.MalformedProduct, "utf-8"):
            si.check_product_contract(self.dir, "2RH1")
        os.remove(os.path.join(self.dir, "2RH1", "summary.yaml"))
        with self.assertRaisesRegex(si.MalformedProduct, "no summary.yaml"):
            si.check_product_contract(self.dir, "2RH1")

    def test_an_instance_yaml_that_is_not_utf8_is_malformed(self):
        path = os.path.join(self.dir, "x.yaml")
        with open(path, "wb") as fh:
            fh.write(b"result: \xff\n")
        with self.assertRaisesRegex(si.MalformedProduct, "utf-8"):
            si.read_instance_rows(path)

    def test_import_structure_checks_the_contract_before_writing(self):
        # import_structure needs the database; this pins the call statically,
        # on the syntax tree, so a commented-out call does not count.
        import ast
        import inspect
        import textwrap

        tree = ast.parse(textwrap.dedent(inspect.getsource(si.import_structure)))
        calls = [(n.lineno, n.func) for n in ast.walk(tree) if isinstance(n, ast.Call)]
        checks = [
            line
            for line, f in calls
            if isinstance(f, ast.Name) and f.id == "check_product_contract"
        ]
        deletes = [
            line
            for line, f in calls
            if isinstance(f, ast.Attribute) and f.attr == "delete"
        ]
        self.assertEqual(len(checks), 1)
        self.assertTrue(deletes)
        self.assertLess(checks[0], min(deletes))


class SeedTests(unittest.TestCase):
    """The seed migration creates every slug the imports can write."""

    @staticmethod
    def _migration():
        import importlib

        return importlib.import_module(
            "interaction.migrations.0009_seed_schrodinger_interaction_types"
        )

    def test_the_seed_covers_every_slug_the_imports_write(self):
        slugs = [row[0] for row in self._migration().SEEDED_TYPES]
        self.assertEqual(len(slugs), len(set(slugs)), "a slug is seeded twice")
        self.assertEqual(set(slugs), set(si.required_slugs()))
        # Creation order sets the ids of the rows a legacy database lacks.
        self.assertEqual(
            slugs[:4],
            [
                "aro_ion_protein",
                "halogen_protein",
                "metal_coordination_protein",
                "covalent",
            ],
        )

    def test_the_seed_keeps_the_legacy_names(self):
        # Written out, not read back from the module: a page shows the name
        # and leaves "hidden" types out, so a changed value must fail here.
        self.assertEqual(
            sorted(self._migration().SEEDED_TYPES),
            sorted(
                [
                    ("acc", "accessible", "hidden", ""),
                    (
                        "aro_ef_protein",
                        "aromatic (edge-to-face)",
                        "aromatic",
                        "protein",
                    ),
                    ("aro_ff", "aromatic (face-to-face)", "aromatic", "none"),
                    ("aro_ion_protein", "aromatic (pi-cation)", "aromatic", "protein"),
                    ("covalent", "covalent bond", "covalent", ""),
                    ("halogen_protein", "halogen contact", "polar", ""),
                    ("hyd", "hydrophobic", "hydrophobic", ""),
                    ("metal_coordination_protein", "metal coordination", "polar", ""),
                    (
                        "polar_acceptor_protein",
                        "polar (hydrogen bond)",
                        "polar",
                        "protein",
                    ),
                    (
                        "polar_backbone",
                        "polar (hydrogen bond with backbone)",
                        "polar",
                        "protein",
                    ),
                    (
                        "polar_donor_protein",
                        "polar (hydrogen bond)",
                        "polar",
                        "protein",
                    ),
                    ("polar_double_neg_protein", "polar (charge-charge)", "polar", ""),
                    ("polar_double_pos_protein", "polar (charge-charge)", "polar", ""),
                    ("Van der Waals", "Van der Waals", "waals", ""),
                ]
            ),
        )

    def test_the_seed_renames_an_old_halogen_bond_row_and_nothing_else(self):
        mig = self._migration()
        saved = []

        class _Row:
            def __init__(self, slug, name):
                self.slug, self.name = slug, name

            def save(self, update_fields):
                saved.append((self.slug, self.name, update_fields))

        existing = {
            "halogen_protein": _Row("halogen_protein", "halogen bond"),
            "hyd": _Row("hyd", "something else"),
        }

        class _Objects:
            def get_or_create(self, slug, defaults):
                if slug in existing:
                    return existing[slug], False
                return _Row(slug, defaults["name"]), True

        class _Model:
            objects = _Objects()

        class _Apps:
            def get_model(self, app, name):
                return _Model

        mig.seed(_Apps(), None)
        self.assertEqual(saved, [("halogen_protein", "halogen contact", ["name"])])
        self.assertEqual(existing["hyd"].name, "something else")

    def test_the_seed_follows_0008_and_reverses_as_a_no_op(self):
        # The remedy _check_slugs prints (migrate interaction 0008, then
        # migrate interaction) needs the reverse to be a no-op.
        from django.db import migrations

        mig = self._migration().Migration
        self.assertEqual(mig.dependencies, [("interaction", "0008_auto_20260921_1803")])
        self.assertEqual(len(mig.operations), 1)
        self.assertIs(mig.operations[0].reverse_code, migrations.RunPython.noop)

    def test_the_seed_hands_every_row_to_the_orm(self):
        mig = self._migration()
        seen = []

        class _Objects:
            def get_or_create(self, **kwargs):
                seen.append(kwargs)
                return object(), True

        class _Model:
            objects = _Objects()

        class _Apps:
            def get_model(self, app, name):
                assert (app, name) == ("interaction", "ResidueFragmentInteractionType")
                return _Model

        mig.seed(_Apps(), None)
        self.assertEqual(
            [
                (
                    k["slug"],
                    k["defaults"]["name"],
                    k["defaults"]["type"],
                    k["defaults"]["direction"],
                )
                for k in seen
            ],
            list(mig.SEEDED_TYPES),
        )


if __name__ == "__main__":
    unittest.main()
