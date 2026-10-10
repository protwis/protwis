"""The Schrodinger ligand imports, as build_all and build_all_interactions run them.

Ligand interactions come from two imports of the Schrodinger deliveries:
Engine 1 serves the anchors named by a HET code (import_schrodinger_interactions),
Engine 2 the "pep" chains (import_schrodinger_peptides). The steps:

1. build the maps that tell each import which product answers which database
   anchor, from this build's annotation and structure text
   (build_schrodinger_chainmap_files, build_schrodinger_peptide_maps);
2. dry-run both imports, so a structure that would fail stops the caller with
   nothing imported;
3. run both imports;
4. remove the maps (remove_schrodinger_maps); after a failed import they stay.

Before anything else a caller runs, with_tests() puts the unit tests of the
import code (run_ligand_import_tests): a build whose import code fails its own
tests stops with nothing written. build_all_interactions runs the tests and
steps 1-2 before its contact-network pass, and 3-4 after it.
"""

import datetime
import os

from django.conf import settings
from django.core.management import call_command
from django.core.management.base import CommandError

# Where the Engine 1 products are delivered, relative to DATA_DIR. The build
# writes the per-PDB chain maps there and removes them once both imports
# succeed.
ENGINE1_DIR = os.sep.join(["structure_data", "schrodinger", "engine1"])
# Where the Engine 2 products are delivered; the per-PDB peptide maps come and
# go the same way.
ENGINE2_DIR = os.sep.join(["structure_data", "schrodinger", "engine2"])

# Where each import run leaves its accounting, relative to BASE_DIR (not in
# DATA_DIR, which is shared input). One directory per run.
ENGINE1_RUN_DIR = os.sep.join(["logs", "engine1_import"])
ENGINE2_RUN_DIR = os.sep.join(["logs", "engine2_peptide_import"])

COMMANDS = ("import_schrodinger_interactions", "import_schrodinger_peptides")
MAP_COMMANDS = ("build_schrodinger_chainmap_files", "build_schrodinger_peptide_maps")
CLEAN_COMMAND = "remove_schrodinger_maps"
TEST_COMMAND = "run_ligand_import_tests"

# The deliveries each command reads. The peptide maps read Engine 1's tree too:
# the producer delivers the coordinate indexes there.
READS = {
    COMMANDS[0]: ("Engine 1",),
    COMMANDS[1]: ("Engine 2",),
    MAP_COMMANDS[0]: ("Engine 1",),
    MAP_COMMANDS[1]: ("Engine 2", "Engine 1"),
    CLEAN_COMMAND: ("Engine 1", "Engine 2"),
    TEST_COMMAND: (),
}


def _now():
    return datetime.datetime.strftime(datetime.datetime.now(), "%Y-%m-%d %H:%M:%S")


def add_arguments(parser):
    """The options both callers take."""
    parser.add_argument(
        "--engine1_data_dir",
        action="store",
        dest="engine1_data_dir",
        default=None,
        help="Engine 1 product tree; default DATA_DIR/" + ENGINE1_DIR,
    )
    parser.add_argument(
        "--engine1_report_dir",
        action="store",
        dest="engine1_report_dir",
        default=None,
        help="Where this run leaves its Engine 1 import accounting; "
        "default BASE_DIR/" + ENGINE1_RUN_DIR + "/<timestamp>",
    )
    parser.add_argument(
        "--engine2_data_dir",
        action="store",
        dest="engine2_data_dir",
        default=None,
        help="Engine 2 product tree; default DATA_DIR/" + ENGINE2_DIR,
    )
    parser.add_argument(
        "--engine2_report_dir",
        action="store",
        dest="engine2_report_dir",
        default=None,
        help="Where this run leaves its Engine 2 peptide import accounting; "
        "default BASE_DIR/" + ENGINE2_RUN_DIR + "/<timestamp>",
    )
    parser.add_argument(
        "--skip_ligand_import",
        action="store_true",
        dest="skip_ligand_import",
        default=False,
        help="Do not import the Schrodinger ligand interactions (Engine 1 "
        'small molecules, Engine 2 "pep" chains) in this run; nothing '
        "else computes them, so the ligand tables keep what they hold",
    )


def engine1_dir(options):
    return options["engine1_data_dir"] or os.sep.join([settings.DATA_DIR, ENGINE1_DIR])


def engine2_dir(options):
    return options["engine2_data_dir"] or os.sep.join([settings.DATA_DIR, ENGINE2_DIR])


def steps(options):
    """[[command, options]]: the two maps, both imports as dry runs, both for real,
    the clean-up. with_tests() puts the tests in front."""
    if options["skip_ligand_import"]:
        print(
            "{} SKIPPING the ligand imports: no ligand interactions are written".format(
                _now()
            )
        )
        return []
    stamp = datetime.datetime.utcnow().strftime("%Y%m%dT%H%M%SZ")
    imports = []
    for command, data_dir, report_dir, default_dir in (
        (
            COMMANDS[0],
            engine1_dir(options),
            options["engine1_report_dir"],
            ENGINE1_RUN_DIR,
        ),
        (
            COMMANDS[1],
            engine2_dir(options),
            options["engine2_report_dir"],
            ENGINE2_RUN_DIR,
        ),
    ):
        run_dir = report_dir or os.sep.join([settings.BASE_DIR, default_dir, stamp])
        print("{} {} accounting goes to {}".format(_now(), command, run_dir))
        imports.append((command, data_dir, run_dir))
    dry = [
        [
            command,
            {
                "data_dir": data_dir,
                "dry_run": True,
                "anomaly_csv": os.path.join(run_dir, "anomalies.dryrun.csv"),
                "report_json": os.path.join(run_dir, "report.dryrun.json"),
            },
        ]
        for command, data_dir, run_dir in imports
    ]
    real = [
        [
            command,
            {
                "data_dir": data_dir,
                "anomaly_csv": os.path.join(run_dir, "anomalies.csv"),
                "report_json": os.path.join(run_dir, "report.json"),
            },
        ]
        for command, data_dir, run_dir in imports
    ]
    # allow_stray: a product directory the annotation does not list gets no
    # map; the database, built from the same annotation, has no anchor there.
    maps = [
        [MAP_COMMANDS[0], {"data_dir": engine1_dir(options), "allow_stray": True}],
        [
            MAP_COMMANDS[1],
            {"data_dir": engine2_dir(options), "index_dir": engine1_dir(options)},
        ],
    ]
    clean = [
        [
            CLEAN_COMMAND,
            {"engine1_dir": engine1_dir(options), "engine2_dir": engine2_dir(options)},
        ]
    ]
    return maps + dry + real + clean


def with_tests(planned):
    """``planned`` [[command, ...]], led by the import tests when it holds an import."""
    if any(step[0] in COMMANDS for step in planned):
        return [[TEST_COMMAND, {}]] + planned
    return planned


def check_deliveries(options, command_names):
    """Refuse to start when a delivery a command in ``command_names`` reads is missing."""
    needed = {label for command in command_names for label in READS.get(command, ())}
    for label, data_dir in (
        ("Engine 1", engine1_dir(options)),
        ("Engine 2", engine2_dir(options)),
    ):
        if label in needed and not os.path.isdir(data_dir):
            raise CommandError(
                "{} products are missing: {} is not a directory. Deliver them, or pass "
                "--skip_ligand_import to go on without ligand interactions.".format(
                    label, data_dir
                )
            )


def split(planned):
    """(tests, maps and dry runs; real imports and the clean-up) of a planned list.

    Each part keeps its order. Plan once and split, so the dry runs and the
    imports they vouch for share one accounting directory per import, and read the
    same maps.
    """

    def after(step):
        return step[0] == CLEAN_COMMAND or (
            step[0] in COMMANDS and not step[1].get("dry_run")
        )

    return ([s for s in planned if not after(s)], [s for s in planned if after(s)])


def run(planned):
    """Run [[command, options]] in order; a failing import raises and stops the caller."""
    for command, kwargs in planned:
        print(
            "{} Running {}{}".format(
                _now(), command, " (dry run)" if kwargs.get("dry_run") else ""
            )
        )
        call_command(command, **kwargs)
