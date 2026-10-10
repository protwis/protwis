"""Run the unit tests of the Schrodinger ligand imports; fail if any test fails.

build_all and build_all_interactions run this before anything else
(build.ligand_imports.with_tests), so a build whose import code fails its own
tests stops before anything is written. The tests need no database and no
delivery.

The suite runs in a child process: a test that patches a module and fails to
restore it cannot affect the imports that run later in the build's process.

    python manage.py run_ligand_import_tests
"""

import contextlib
import io
import os
import subprocess
import sys
import unittest

from django.conf import settings
from django.core.management.base import BaseCommand, CommandError

TEST_MODULES = (
    "build.test_build_all_ligands",
    "interaction.test_schrodinger_chain_map",
    "interaction.test_schrodinger_chainmap_files",
    "interaction.test_schrodinger_complex",
    "interaction.test_schrodinger_import",
    "interaction.test_schrodinger_map_builders",
    "interaction.test_schrodinger_peptide",
    "interaction.test_stored_interactions",
)

# A suite that has not finished by then is treated as failed.
TIMEOUT_SECONDS = 1800


class Command(BaseCommand):
    help = (
        "Run the unit tests of the Schrodinger ligand imports in a child process; "
        "fail if any test fails."
    )
    # The system checks load every URL module, whose views query the database;
    # the tests need none of it.
    requires_system_checks = []

    def add_arguments(self, parser):
        parser.add_argument(
            "--in-process",
            action="store_true",
            dest="in_process",
            help="Run the tests in this process (what the child process does).",
        )

    def handle(self, *args, **options):
        if options["in_process"]:
            self.run_suite()
            return
        manage = os.path.join(settings.BASE_DIR, "manage.py")
        # The child writes to the same stream; what this process printed comes first.
        sys.stdout.flush()
        sys.stderr.flush()
        try:
            done = subprocess.run(
                [sys.executable, manage, "run_ligand_import_tests", "--in-process"],
                env=dict(os.environ, PYTHONWARNINGS="ignore"),
                timeout=TIMEOUT_SECONDS,
            )
        except subprocess.TimeoutExpired:
            raise CommandError(
                "the ligand import tests did not finish in {} s".format(TIMEOUT_SECONDS)
            )
        if done.returncode:
            raise CommandError(
                "the ligand import tests failed (exit {}); see the output above".format(
                    done.returncode
                )
            )

    def run_suite(self):
        suite = unittest.defaultTestLoader.loadTestsFromNames(TEST_MODULES)
        out = io.StringIO()
        # What the tests print goes into the report, shown only on failure.
        with contextlib.redirect_stdout(out), contextlib.redirect_stderr(out):
            result = unittest.TextTestRunner(stream=out, verbosity=1).run(suite)
        if not result.wasSuccessful():
            self.stdout.write(out.getvalue())
            raise CommandError(
                "{} of {} ligand import tests failed".format(
                    len(result.failures) + len(result.errors), result.testsRun
                )
            )
        if not result.testsRun:
            raise CommandError("no ligand import test was found")
        self.stdout.write(
            "ligand import tests: {} run, all passed".format(result.testsRun)
        )
