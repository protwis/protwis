"""Unit tests for the ligand import steps of build_all and build_all_interactions."""

import io
import os
import tempfile
import unittest
from unittest import mock

from django.core.management.base import CommandError

from build import ligand_imports
from build.management.commands import build_all


def options(**kw):
    out = {
        "skip_ligand_import": False,
        "engine1_data_dir": "/e1",
        "engine2_data_dir": "/e2",
        "engine1_report_dir": "/r1",
        "engine2_report_dir": "/r2",
        "phase": None,
    }
    out.update(kw)
    return out


def no_import(name, **kw):
    raise AssertionError("a real import command was reached: %s" % name)


class LigandImportStepsTests(unittest.TestCase):
    def test_maps_then_both_dry_runs_then_both_imports(self):
        steps = build_all.Command().ligand_import_steps(options())
        self.assertEqual(
            [(c, o.get("dry_run", False), o.get("data_dir")) for c, o in steps],
            [
                ("build_schrodinger_chainmap_files", False, "/e1"),
                ("build_schrodinger_peptide_maps", False, "/e2"),
                ("import_schrodinger_interactions", True, "/e1"),
                ("import_schrodinger_peptides", True, "/e2"),
                ("import_schrodinger_interactions", False, "/e1"),
                ("import_schrodinger_peptides", False, "/e2"),
                ("remove_schrodinger_maps", False, None),
            ],
        )
        self.assertEqual(
            [o["report_json"] for c, o in steps if c in ligand_imports.COMMANDS],
            [
                "/r1/report.dryrun.json",
                "/r2/report.dryrun.json",
                "/r1/report.json",
                "/r2/report.json",
            ],
        )

    def test_the_peptide_maps_read_the_indexes_from_the_engine1_tree(self):
        steps = dict(
            (c, o)
            for c, o in build_all.Command().ligand_import_steps(options())
            if c in ligand_imports.MAP_COMMANDS
        )
        self.assertEqual(
            steps["build_schrodinger_peptide_maps"],
            {"data_dir": "/e2", "index_dir": "/e1"},
        )
        # The annotation decides which structures get a map; a product directory
        # it does not list is not a reason to stop the build.
        self.assertEqual(
            steps["build_schrodinger_chainmap_files"],
            {"data_dir": "/e1", "allow_stray": True},
        )

    def test_split_runs_the_maps_with_the_dry_runs(self):
        before, imports = ligand_imports.split(ligand_imports.steps(options()))
        self.assertEqual(
            [(c, o.get("dry_run", False)) for c, o in before],
            [
                ("build_schrodinger_chainmap_files", False),
                ("build_schrodinger_peptide_maps", False),
                ("import_schrodinger_interactions", True),
                ("import_schrodinger_peptides", True),
            ],
        )
        self.assertEqual(
            [(c, o.get("dry_run", False)) for c, o in imports],
            [
                ("import_schrodinger_interactions", False),
                ("import_schrodinger_peptides", False),
                ("remove_schrodinger_maps", False),
            ],
        )

    def test_the_clean_up_names_both_trees(self):
        steps = build_all.Command().ligand_import_steps(options())
        self.assertEqual(
            steps[-1],
            ["remove_schrodinger_maps", {"engine1_dir": "/e1", "engine2_dir": "/e2"}],
        )

    def test_the_peptide_maps_alone_need_the_engine1_delivery(self):
        with tempfile.TemporaryDirectory() as d2:
            with self.assertRaisesRegex(CommandError, "Engine 1 products are missing"):
                ligand_imports.check_deliveries(
                    options(engine1_data_dir="/nonexistent-e1", engine2_data_dir=d2),
                    ["build_schrodinger_peptide_maps"],
                )
            ligand_imports.check_deliveries(
                options(engine1_data_dir="/nonexistent-e1", engine2_data_dir=d2),
                ["import_schrodinger_peptides"],
            )

    def test_the_tests_lead_any_plan_with_an_import_and_need_no_delivery(self):
        planned = ligand_imports.with_tests(ligand_imports.steps(options()))
        self.assertEqual(planned[0], ["run_ligand_import_tests", {}])
        self.assertEqual(planned[1:], ligand_imports.steps(options()))
        self.assertEqual(ligand_imports.with_tests([]), [])
        self.assertEqual(
            ligand_imports.with_tests([["build_common"]]), [["build_common"]]
        )
        ligand_imports.check_deliveries(
            options(
                engine1_data_dir="/nonexistent-e1", engine2_data_dir="/nonexistent-e2"
            ),
            ["run_ligand_import_tests"],
        )

    def test_a_failing_test_stops_the_steps_before_the_maps(self):
        calls = []

        def tests_fail(name, **kw):
            calls.append(name)
            if name == "run_ligand_import_tests":
                raise CommandError("1 of 220 ligand import tests failed")

        with mock.patch.object(ligand_imports, "call_command", tests_fail):
            with self.assertRaisesRegex(CommandError, "tests failed"):
                ligand_imports.run(
                    ligand_imports.with_tests(ligand_imports.steps(options()))
                )
        self.assertEqual(calls, ["run_ligand_import_tests"])

    def first_command_of_build_all(self, **kw):
        cmd = build_all.Command()
        opts = vars(cmd.create_parser("manage.py", "build_all").parse_args([]))
        calls = []

        class Stop(Exception):
            pass

        def first(name, *a, **k):
            calls.append(name)
            raise Stop

        with tempfile.TemporaryDirectory() as d1, tempfile.TemporaryDirectory() as d2:
            opts.update(options(engine1_data_dir=d1, engine2_data_dir=d2))
            opts.update(kw)
            with mock.patch.object(build_all, "call_command", first):
                with self.assertRaises(Stop):
                    cmd.handle(**opts)
        return calls[0]

    def test_build_all_runs_the_tests_before_anything_else(self):
        self.assertEqual(self.first_command_of_build_all(), "run_ligand_import_tests")
        self.assertEqual(
            self.first_command_of_build_all(phase=1), "run_ligand_import_tests"
        )

    def test_build_all_without_the_imports_runs_no_tests(self):
        self.assertEqual(
            self.first_command_of_build_all(phase=2), "build_structure_angles"
        )
        self.assertEqual(
            self.first_command_of_build_all(skip_ligand_import=True), "clear_cache"
        )

    def test_skipping_imports_nothing(self):
        self.assertEqual(
            build_all.Command().ligand_import_steps(options(skip_ligand_import=True)),
            [],
        )


class BuildAllInteractionsTests(unittest.TestCase):
    """
    build_all_interactions runs the contact network and the same imports.

    The imports run dry before the contacts and for real after them.
    """

    def test_it_takes_the_import_options_and_dry_runs_before_the_contacts(self):
        from tools.management.commands import build_all_interactions as bai

        cmd = bai.Command()
        parser = cmd.create_parser("manage.py", "build_all_interactions")
        opts = vars(parser.parse_args([]))
        for key in ("engine1_data_dir", "engine2_data_dir", "skip_ligand_import"):
            self.assertIn(key, opts)
        calls = []
        with tempfile.TemporaryDirectory() as d1, tempfile.TemporaryDirectory() as d2:
            opts.update(options(engine1_data_dir=d1, engine2_data_dir=d2))
            with mock.patch.object(
                bai.Command,
                "prepare_input",
                lambda self, proc, pdbs: calls.append("contacts"),
            ), mock.patch.object(
                ligand_imports,
                "call_command",
                lambda name, **kw: calls.append((name, kw.get("dry_run", False))),
            ):
                cmd.handle(**opts)
        self.assertEqual(
            calls,
            [
                ("run_ligand_import_tests", False),
                ("build_schrodinger_chainmap_files", False),
                ("build_schrodinger_peptide_maps", False),
                ("import_schrodinger_interactions", True),
                ("import_schrodinger_peptides", True),
                "contacts",
                ("import_schrodinger_interactions", False),
                ("import_schrodinger_peptides", False),
                ("remove_schrodinger_maps", False),
            ],
        )

    def test_a_failing_dry_run_stops_it_before_the_contacts_and_imports_nothing(self):
        from tools.management.commands import build_all_interactions as bai

        calls = []

        def dry_fails(name, **kw):
            calls.append((name, kw.get("dry_run", False)))
            if kw.get("dry_run") and name == "import_schrodinger_peptides":
                raise CommandError("1 structure(s) failed")

        with tempfile.TemporaryDirectory() as d1, tempfile.TemporaryDirectory() as d2:
            with mock.patch.object(
                bai.Command,
                "prepare_input",
                lambda self, proc, pdbs: calls.append("contacts"),
            ), mock.patch.object(ligand_imports, "call_command", dry_fails):
                with self.assertRaisesRegex(CommandError, "failed"):
                    bai.Command().handle(
                        **options(engine1_data_dir=d1, engine2_data_dir=d2, proc=1)
                    )
        self.assertEqual(
            calls,
            [
                ("run_ligand_import_tests", False),
                ("build_schrodinger_chainmap_files", False),
                ("build_schrodinger_peptide_maps", False),
                ("import_schrodinger_interactions", True),
                ("import_schrodinger_peptides", True),
            ],
        )

    def test_a_failing_import_after_the_contacts_makes_the_command_fail(self):
        from tools.management.commands import build_all_interactions as bai

        calls = []

        def real_fails(name, **kw):
            calls.append((name, kw.get("dry_run", False)))
            if name in ligand_imports.COMMANDS and not kw.get("dry_run"):
                raise CommandError("1 structure(s) failed")

        with tempfile.TemporaryDirectory() as d1, tempfile.TemporaryDirectory() as d2:
            with mock.patch.object(
                bai.Command,
                "prepare_input",
                lambda self, proc, pdbs: calls.append("contacts"),
            ), mock.patch.object(ligand_imports, "call_command", real_fails):
                with self.assertRaisesRegex(CommandError, "failed"):
                    bai.Command().handle(
                        **options(engine1_data_dir=d1, engine2_data_dir=d2, proc=1)
                    )
        self.assertEqual(
            calls,
            [
                ("run_ligand_import_tests", False),
                ("build_schrodinger_chainmap_files", False),
                ("build_schrodinger_peptide_maps", False),
                ("import_schrodinger_interactions", True),
                ("import_schrodinger_peptides", True),
                "contacts",
                ("import_schrodinger_interactions", False),
            ],
        )

    def test_a_contact_network_that_raises_does_not_stop_the_imports(self):
        # The contact network's failures are printed and logged, and the ligand
        # imports still run.
        from tools.management.commands import build_all_interactions as bai

        calls = []

        def contacts_fail(self, proc, pdbs):
            calls.append("contacts")
            raise RuntimeError("contact network failed")

        with tempfile.TemporaryDirectory() as d1, tempfile.TemporaryDirectory() as d2:
            with mock.patch.object(
                bai.Command, "prepare_input", contacts_fail
            ), mock.patch.object(
                ligand_imports,
                "call_command",
                lambda name, **kw: calls.append((name, kw.get("dry_run", False))),
            ):
                bai.Command().handle(
                    **options(engine1_data_dir=d1, engine2_data_dir=d2, proc=1)
                )
        self.assertEqual(
            calls,
            [
                ("run_ligand_import_tests", False),
                ("build_schrodinger_chainmap_files", False),
                ("build_schrodinger_peptide_maps", False),
                ("import_schrodinger_interactions", True),
                ("import_schrodinger_peptides", True),
                "contacts",
                ("import_schrodinger_interactions", False),
                ("import_schrodinger_peptides", False),
                ("remove_schrodinger_maps", False),
            ],
        )

    def test_each_import_keeps_its_own_accounting_directory_from_one_plan(self):
        from tools.management.commands import build_all_interactions as bai

        paths = []
        with tempfile.TemporaryDirectory() as d1, tempfile.TemporaryDirectory() as d2:
            with mock.patch.object(
                bai.Command, "prepare_input", lambda self, proc, pdbs: None
            ), mock.patch.object(
                ligand_imports, "steps", wraps=ligand_imports.steps
            ) as plan, mock.patch.object(
                ligand_imports,
                "call_command",
                lambda name, **kw: (
                    paths.append((name, kw["report_json"], kw["anomaly_csv"]))
                    if name in ligand_imports.COMMANDS
                    else None
                ),
            ):
                bai.Command().handle(
                    **options(
                        engine1_data_dir=d1,
                        engine2_data_dir=d2,
                        engine1_report_dir=None,
                        engine2_report_dir=None,
                        proc=1,
                    )
                )
        self.assertEqual(plan.call_count, 1)
        self.assertEqual(len(paths), 4)
        dirs = {}
        for name, report, anomalies in paths:
            # A run's report and anomaly list sit side by side.
            self.assertEqual(os.path.dirname(report), os.path.dirname(anomalies))
            dirs.setdefault(name, set()).add(os.path.dirname(report))
        self.assertEqual(
            sorted(dirs),
            ["import_schrodinger_interactions", "import_schrodinger_peptides"],
        )
        for name, found in dirs.items():
            self.assertEqual(len(found), 1, name)
        # The two imports do not share a directory, and no run overwrites another's file.
        self.assertNotEqual(
            dirs["import_schrodinger_interactions"], dirs["import_schrodinger_peptides"]
        )
        files = [p for _n, report, anomalies in paths for p in (report, anomalies)]
        self.assertEqual(len(set(files)), 8)

    def test_skip_is_off_by_default_and_skipping_needs_no_delivery(self):
        from tools.management.commands import build_all_interactions as bai

        cmd = bai.Command()
        opts = vars(
            cmd.create_parser("manage.py", "build_all_interactions").parse_args([])
        )
        self.assertIs(opts["skip_ligand_import"], False)
        calls = []
        opts.update(
            options(
                skip_ligand_import=True,
                engine1_data_dir="/nonexistent-e1",
                engine2_data_dir="/nonexistent-e2",
            )
        )
        with mock.patch.object(
            bai.Command,
            "prepare_input",
            lambda self, proc, pdbs: calls.append("contacts"),
        ), mock.patch.object(
            ligand_imports, "call_command", lambda name, **kw: calls.append(name)
        ):
            cmd.handle(**opts)
        self.assertEqual(calls, ["contacts"])

    def test_a_missing_engine2_delivery_alone_stops_it_before_the_contacts(self):
        from tools.management.commands import build_all_interactions as bai

        calls = []
        with tempfile.TemporaryDirectory() as d1:
            with mock.patch.object(
                bai.Command,
                "prepare_input",
                lambda self, proc, pdbs: calls.append("contacts"),
            ), mock.patch.object(ligand_imports, "call_command", no_import):
                with self.assertRaisesRegex(
                    CommandError, "Engine 2 products are missing"
                ):
                    bai.Command().handle(
                        **options(
                            engine1_data_dir=d1,
                            engine2_data_dir="/nonexistent-e2",
                            proc=1,
                        )
                    )
        self.assertEqual(calls, [])

    def test_a_missing_delivery_stops_it_before_the_contacts(self):
        from tools.management.commands import build_all_interactions as bai

        calls = []
        with mock.patch.object(
            bai.Command,
            "prepare_input",
            lambda self, proc, pdbs: calls.append("contacts"),
        ), mock.patch.object(ligand_imports, "call_command", no_import):
            with self.assertRaisesRegex(CommandError, "Engine 1 products are missing"):
                bai.Command().handle(
                    **options(
                        engine1_data_dir="/nonexistent-e1",
                        engine2_data_dir="/nonexistent-e2",
                        proc=1,
                    )
                )
        self.assertEqual(calls, [])


class RunLigandImportTestsTests(unittest.TestCase):
    """The command fails when the suite does, in the child process or in this one."""

    def command(self):
        from build.management.commands import run_ligand_import_tests

        return run_ligand_import_tests

    def test_a_failing_child_process_is_a_command_error(self):
        mod = self.command()
        for code, fails in ((0, False), (1, True)):
            done = mock.Mock(returncode=code)
            with mock.patch.object(mod.subprocess, "run", return_value=done) as run:
                if fails:
                    with self.assertRaisesRegex(CommandError, "tests failed"):
                        mod.Command().handle(in_process=False)
                else:
                    mod.Command().handle(in_process=False)
            self.assertIn("--in-process", run.call_args[0][0])
            self.assertEqual(run.call_args[1]["timeout"], mod.TIMEOUT_SECONDS)

    def test_a_failing_suite_is_a_command_error(self):
        mod = self.command()

        class Fails(unittest.TestCase):
            def runTest(self):
                self.fail("on purpose")

        class Passes(unittest.TestCase):
            def runTest(self):
                pass

        class Errors(unittest.TestCase):
            def runTest(self):
                raise ImportError("on purpose")

        for cases, fails in (
            ([Passes()], None),
            ([Fails()], "1 of 1"),
            ([Errors()], "1 of 1"),
            ([], "no ligand import test"),
        ):
            suite = unittest.TestSuite(cases)
            with mock.patch.object(
                mod.unittest.defaultTestLoader, "loadTestsFromNames", return_value=suite
            ) as load:
                cmd = mod.Command(stdout=io.StringIO())
                if fails:
                    with self.assertRaisesRegex(CommandError, fails):
                        cmd.handle(in_process=True)
                else:
                    cmd.handle(in_process=True)
            load.assert_called_once_with(mod.TEST_MODULES)

    def test_a_child_that_hangs_is_a_command_error(self):
        mod = self.command()
        hang = mod.subprocess.TimeoutExpired(cmd="x", timeout=1)
        with mock.patch.object(mod.subprocess, "run", side_effect=hang):
            with self.assertRaisesRegex(CommandError, "did not finish"):
                mod.Command().handle(in_process=False)

    def test_every_test_module_of_the_imports_is_listed(self):
        import glob

        here = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        uses = ("schrodinger", "ligand_imports", "stored_interactions")
        found = set()
        for app in ("build", "interaction"):
            for path in glob.glob(os.path.join(here, app, "test_*.py")):
                with open(path) as fh:
                    text = fh.read()
                if any(word in text for word in uses):
                    found.add(os.path.relpath(path, here)[:-3].replace(os.sep, "."))
        self.assertEqual(found, set(self.command().TEST_MODULES))


if __name__ == "__main__":
    unittest.main()
