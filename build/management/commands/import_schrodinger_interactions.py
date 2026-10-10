"""Import Schrodinger Engine 1 interaction YAMLs, one transaction per structure.

Usage, the way a build calls it -- the tree under --data-dir says what to
import, and build_schrodinger_chainmap_files has written a chainmap.tsv into
each structure directory, so nothing else has to be passed::

    python manage.py import_schrodinger_interactions \\
        --data-dir DATA_DIR/structure_data/schrodinger/engine1 \\
        --anomaly-csv /runs/anomalies.csv

Every directory under --data-dir is one structure; whether Engine 1 ran it is
read from its chain map and checked against the tree (see not_run below).
--pdb / --pdb-list narrow that to a few.

Each chainmap.tsv carries its own provenance, so a tree merged from several
deliveries is valid; the command reports the distinct provenance it saw.

A directory with no chainmap.tsv fails that structure and the run exits
non-zero; it is never passed over, or its anchors would keep their old rows.
For the same reason, when the corpus comes from the tree, a database structure
with Engine 1 anchors and no directory at all stops the run before anything is
touched.

Each structure is imported in its own transaction. A structure that fails
(unreadable YAML, a summary.yaml without the contract version this importer
reads, a row the type map cannot route, or any unexpected error) is rolled
back and reported; the others are unaffected. The command exits
non-zero when any structure failed, after all structures have been attempted.

An anchor with no product instance loses its existing rows and is reported
as a WARNING (anchor_cleared) with the map's reason. That applies only to a
structure Engine 1 ran: one whose directory holds neither a product instance
nor the producer's summary is left untouched and reported as not_run. When
such a structure has anchors, they would keep rows this import did not write,
so it fails the run; --allow-not-run accepts it.

Point the anomaly CSV and the report JSON at a fresh directory each run
(build_all gives each build a timestamped one under logs/), not into the
delivered tree. The anomaly CSV is written outside the transactions and
flushed per row, so it survives any rollback; every row that was read but not
written is accounted for in it.
"""

import collections
import csv
import datetime
import json
import os

from django.core.management.base import BaseCommand, CommandError
from django.db import transaction

from interaction import schrodinger_import as si
from interaction.models import (
    ResidueFragmentInteractionType,
    StructureLigandInteraction,
)
from structure.models import Structure


class _Rollback(Exception):
    """Raised inside a dry-run transaction to undo it."""


def _make_room_for(path):
    """Create the directory an output file is about to be written into.

    Each run is given its own directory, so the caller should not have to make
    it first.
    """
    parent = os.path.dirname(os.path.abspath(path))
    if parent:
        os.makedirs(parent, exist_ok=True)


class AnomalyLog(object):
    """CSV sink for per-anchor accounting; independent of the database."""

    COLUMNS = [
        "timestamp",
        "pdb",
        "sli_id",
        "het",
        "level",
        "category",
        "count",
        "detail",
    ]

    def __init__(self, path):
        _make_room_for(path)
        self._fh = open(path, "w", newline="")
        self._writer = csv.writer(self._fh)
        self._writer.writerow(self.COLUMNS)
        self._fh.flush()
        self.levels = {"INFO": 0, "WARNING": 0, "ERROR": 0}

    def log(self, pdb, level, category, sli_id="", het="", count=1, detail=""):
        self.levels[level] += 1
        self._writer.writerow(
            [
                datetime.datetime.utcnow().isoformat(timespec="seconds"),
                pdb,
                sli_id,
                het,
                level,
                category,
                count,
                detail,
            ]
        )
        self._fh.flush()

    def close(self):
        self._fh.close()


class Command(BaseCommand):
    help = "Import Schrodinger Engine 1 interactions (one transaction per structure)."

    def add_arguments(self, parser):
        parser.add_argument(
            "--data-dir",
            required=True,
            help="Root of the product tree: {data_dir}/{PDB}/{instance}/.",
        )
        parser.add_argument(
            "--pdb", action="append", default=[], help="PDB code to import; repeatable."
        )
        parser.add_argument(
            "--pdb-list",
            default=None,
            help="File with one PDB code per line (# comments allowed). "
            "Without it and without --pdb, every structure directory "
            "under --data-dir is imported.",
        )
        parser.add_argument(
            "--anomaly-csv",
            required=True,
            help="Where to write the per-anchor accounting CSV.",
        )
        parser.add_argument(
            "--report-json",
            default=None,
            help="Optional path for a machine-readable per-anchor report.",
        )
        parser.add_argument(
            "--allow-not-run",
            action="store_true",
            help="Report, instead of failing, a structure that Engine 1 never "
            "ran but that has anchors: its rows stay as the legacy "
            "pipeline left them.",
        )
        parser.add_argument(
            "--dry-run",
            action="store_true",
            help="Run every structure and roll each one back.",
        )

    def _pdb_codes(self, options):
        """(codes, came_from_the_tree). The second value decides whether the
        database's in-scope structures have to be covered: a list file that
        turns out to hold no code leaves the tree as the corpus, and the check
        has to follow the corpus, not the flag."""
        codes = [c.strip().upper() for c in options["pdb"] if c.strip()]
        if options["pdb_list"]:
            with open(options["pdb_list"]) as fh:
                for line in fh:
                    line = line.split("#", 1)[0].strip()
                    if line:
                        codes.append(line.upper())
        if not codes:
            # Every structure directory the tree offers; which of them the
            # producer ran is decided per structure from its chainmap.
            codes = si.product_pdb_codes(options["data_dir"])
            if not codes:
                raise CommandError(
                    "--data-dir {!r} holds no structure directory".format(
                        options["data_dir"]
                    )
                )
            return list(dict.fromkeys(codes)), True
        return list(dict.fromkeys(codes)), False

    @staticmethod
    def _uncovered(codes):
        """In-scope structures the database has and the delivered tree does not.

        Only meaningful when the corpus came from the tree: such a structure is
        never visited, so its anchors would keep rows this import did not
        write.
        """
        have = {c.upper() for c in codes}
        out = []
        for sli in StructureLigandInteraction.objects.select_related(
            "structure__pdb_code", "structure__structure_type", "ligand__ligand_type"
        ).order_by("id"):
            if sli.structure is None or not si.is_in_scope(sli):
                continue
            if sli.structure.structure_type.origin != si.STRUCTURE_ORIGIN:
                continue
            pdb = sli.structure.pdb_code.index.upper()
            if pdb not in have:
                out.append(pdb)
        return sorted(set(out))

    def _check_slugs(self):
        present = set(
            ResidueFragmentInteractionType.objects.values_list("slug", flat=True)
        )
        missing = sorted(si.required_slugs() - present)
        if missing:
            # Every type the imports write is seeded by interaction migration
            # 0009; nothing in the build creates one. A plain migrate does not
            # run an applied migration again, and the seed has a no-op reverse,
            # so unapplying and reapplying it is safe.
            raise CommandError(
                "interaction types missing from the database: {} (they are seeded by "
                "interaction migration 0009: run migrate, or, if it is already applied "
                "and the rows were deleted, migrate interaction 0008 and then migrate "
                "interaction)".format(", ".join(missing))
            )

    def handle(self, *args, **options):
        if not os.path.isdir(options["data_dir"]):
            raise CommandError(
                "--data-dir {!r} is not a directory".format(options["data_dir"])
            )
        codes, whole_tree = self._pdb_codes(options)
        self._check_slugs()
        try:
            anchor_map, receptor_map, no_chainmap, provenance, not_run = (
                si.load_chainmap_dir(options["data_dir"], codes)
            )
        except (OSError, si.MapMismatch) as exc:
            raise CommandError("cannot read the chain maps: {}".format(exc))
        for key, values in sorted(provenance.items()):
            if len(values) > 1:
                # Maps built at different times are expected: reported, not refused.
                self.stdout.write(
                    "{}: {} distinct values across the tree ({})".format(
                        key,
                        len(values),
                        ", ".join(
                            "{}x {}".format(n, v or "(blank)")
                            for v, n in sorted(values.items(), key=lambda kv: -kv[1])[
                                :5
                            ]
                        ),
                    )
                )
        if whole_tree:
            uncovered = self._uncovered(codes)
            if uncovered:
                raise CommandError(
                    "{} structure(s) have Engine 1 anchors in the database and no directory "
                    "under --data-dir, so they would keep whatever the legacy pipeline wrote: "
                    "{}{}".format(
                        len(uncovered),
                        ", ".join(uncovered[:20]),
                        ""
                        if len(uncovered) <= 20
                        else " (+%d more)" % (len(uncovered) - 20),
                    )
                )
        log = AnomalyLog(options["anomaly_csv"])
        no_chainmap_set = set(no_chainmap)
        not_run_set = set(not_run)
        report = []
        failed = []
        totals = {
            "anchors": 0,
            "cleared": 0,
            "deleted": 0,
            "written": 0,
            "fragments_deleted": 0,
            "pdbdata_deleted": 0,
            "pdbdata_kept_referenced": 0,
        }
        files = collections.Counter()
        try:
            for pdb in codes:
                entry = {"pdb": pdb}
                report.append(entry)
                structure = (
                    Structure.objects.filter(pdb_code__index__iexact=pdb)
                    .select_related(
                        "structure_type", "pdb_code", "protein_conformation"
                    )
                    .first()
                )
                was_run = pdb not in not_run_set
                at_risk = 0
                if structure is not None and not was_run:
                    at_risk = len(
                        [
                            s
                            for s in StructureLigandInteraction.objects.filter(
                                structure=structure
                            ).select_related("ligand__ligand_type")
                            if si.is_in_scope(s)
                        ]
                    )
                verdict = si.structure_verdict(
                    structure is not None,
                    structure is not None
                    and structure.structure_type.origin == si.STRUCTURE_ORIGIN,
                    was_run,
                    pdb not in no_chainmap_set,
                    at_risk,
                    options["allow_not_run"],
                )
                if verdict is not None:
                    status, level, category = verdict
                    detail = ""
                    if category == "not_experimental":
                        detail = structure.structure_type.slug
                    elif category == "not_run":
                        detail = "no product instance and no {} under {}; {}".format(
                            si.PRODUCT_SUMMARY_NAME,
                            os.path.join(options["data_dir"], pdb),
                            "{} anchor(s) would keep what the legacy pipeline wrote".format(
                                at_risk
                            )
                            if at_risk
                            else "nothing here for Engine 1 anyway",
                        )
                        entry["not_run_anchors"] = at_risk
                        if status == "failed":
                            entry["error"] = detail
                            failed.append(pdb)
                    elif category == "no_chainmap":
                        detail = "no {} under {}".format(
                            si.CHAINMAP_NAME, os.path.join(options["data_dir"], pdb)
                        )
                        entry["error"] = detail
                        failed.append(pdb)
                    log.log(pdb, level, category, count=at_risk, detail=detail)
                    entry["status"] = status
                    continue
                if not os.path.isdir(os.path.join(options["data_dir"], pdb)):
                    log.log(
                        pdb,
                        "WARNING",
                        "no_product_dir",
                        detail="no product directory; the structure has no summary.yaml "
                        "with a contract version, so it fails if it has in-scope "
                        "anchors",
                    )
                try:
                    with transaction.atomic():
                        outcomes, out_of_scope, cleanup, unused = si.import_structure(
                            structure, options["data_dir"], anchor_map, receptor_map
                        )
                        if options["dry_run"]:
                            raise _Rollback()
                except _Rollback:
                    pass
                except Exception as exc:
                    # The structure's transaction has been rolled back; record
                    # the failure and go on with the next structure.
                    message = "{}: {}".format(type(exc).__name__, exc)[:300]
                    log.log(pdb, "ERROR", type(exc).__name__, detail=message)
                    entry["status"] = "failed"
                    entry["error"] = message
                    failed.append(pdb)
                    continue
                entry["status"] = "rolled_back" if options["dry_run"] else "imported"
                entry["out_of_scope_anchors"] = out_of_scope
                if unused:
                    # The map lists ligand copies this database has no anchor
                    # for. One row per HET, so the CSV can be grepped by ligand.
                    # A further copy of a ligand that IS anchored is the normal
                    # shape of the copy axis and only worth an INFO; a HET with
                    # no anchor at all means a ligand was computed and will
                    # never reach the database, which is a WARNING.
                    entry["map_copies_unused"] = [
                        "{} {}".format(h, t).strip() for h, t in unused
                    ]
                    anchored = {o.het for o in outcomes}
                    by_het = collections.OrderedDict()
                    for het, tok in unused:
                        by_het.setdefault(het, []).append(tok or "(no token)")
                    for het, tokens in by_het.items():
                        if het in anchored:
                            log.log(
                                pdb,
                                "INFO",
                                "map_copies_unused",
                                het=het,
                                count=len(tokens),
                                detail="further copies of an anchored ligand, the database "
                                "holds one row for all of them: " + ", ".join(tokens),
                            )
                        else:
                            log.log(
                                pdb,
                                "WARNING",
                                "map_het_not_imported",
                                het=het,
                                count=len(tokens),
                                detail="computed but never imported, the database has no "
                                "anchor for this ligand: " + ", ".join(tokens),
                            )
                entry["cleanup"] = dict(cleanup)
                for key in (
                    "fragments_deleted",
                    "pdbdata_deleted",
                    "pdbdata_kept_referenced",
                ):
                    totals[key] += cleanup[key]
                entry["anchors"] = []
                for o in outcomes:
                    self._log_outcome(log, pdb, o)
                    entry["anchors"].append(
                        {
                            "sli_id": o.sli_id,
                            "het": o.het,
                            "mode": o.mode,
                            "instances": o.instances,
                            "notes": o.notes,
                            "deleted": o.deleted,
                            "written": o.written,
                            "fragments_created": o.fragments_created,
                            "counts": dict(o.counts),
                            "other_chain_by_chain": o.other_chain_by_chain,
                            "dropped": dict(o.dropped),
                            "complex_file": o.complex_file,
                        }
                    )
                    files[o.complex_file] += 1
                    totals["anchors"] += 1
                    totals["cleared"] += o.mode == "no_product"
                    totals["deleted"] += o.deleted
                    totals["written"] += o.written
        finally:
            log.close()
            if options["report_json"]:
                _make_room_for(options["report_json"])
                with open(options["report_json"], "w") as fh:
                    json.dump(
                        {
                            "dry_run": options["dry_run"],
                            "totals": totals,
                            "complex_files": dict(files),
                            "provenance": provenance,
                            "not_run": not_run,
                            "failed": failed,
                            "structures": report,
                        },
                        fh,
                        indent=1,
                    )

        self.stdout.write(
            "{} structures, {} in-scope anchors ({} cleared, no product), {} RFI rows deleted, "
            "{} written; 3D files {}; orphan fragments deleted {}, PdbData deleted {} (kept, "
            "referenced elsewhere: {}){}; anomalies INFO={} WARNING={} ERROR={}".format(
                len(codes),
                totals["anchors"],
                totals["cleared"],
                totals["deleted"],
                totals["written"],
                dict(sorted(files.items())),
                totals["fragments_deleted"],
                totals["pdbdata_deleted"],
                totals["pdbdata_kept_referenced"],
                " (dry run: rolled back)" if options["dry_run"] else "",
                log.levels["INFO"],
                log.levels["WARNING"],
                log.levels["ERROR"],
            )
        )
        if failed:
            raise CommandError(
                "{} structure(s) failed and were rolled back: {}".format(
                    len(failed), ", ".join(failed)
                )
            )

    @staticmethod
    def _log_outcome(log, pdb, o):
        c = o.counts
        if o.mode == "no_product":
            log.log(
                pdb,
                "WARNING",
                "anchor_cleared",
                o.sli_id,
                o.het,
                o.deleted,
                detail="no product instance; {} existing rows deleted; {}".format(
                    o.deleted, "; ".join(o.notes)
                ),
            )
            return
        if o.complex_file == "ligand_not_found":
            log.log(
                pdb,
                "WARNING",
                "complex_file_missing",
                o.sli_id,
                o.het,
                detail="rows written but the ligand was not found in the stored structure "
                "text; the anchor has no 3D file",
            )
        if o.mode == "mapped_partial":
            log.log(
                pdb,
                "WARNING",
                "anchor_instances_partial",
                o.sli_id,
                o.het,
                detail="; ".join(o.notes),
            )
        for note in o.notes:
            if note.startswith("errata: "):
                log.log(pdb, "INFO", "annotation_errata", o.sli_id, o.het, detail=note)
        if c["rows_in"] == 0:
            log.log(
                pdb,
                "INFO",
                "product_has_zero_rows",
                o.sli_id,
                o.het,
                detail=",".join(o.instances),
            )
        for category, level in (
            ("excluded_family", "INFO"),
            ("nonstandard_residue", "INFO"),
            ("duplicate", "INFO"),
            ("other_chain", "WARNING"),
        ):
            if c[category]:
                detail = ""
                if category == "other_chain":
                    detail = ",".join(
                        "{}:{}".format(k, v)
                        for k, v in sorted(o.other_chain_by_chain.items())
                    )
                log.log(pdb, level, category, o.sli_id, o.het, c[category], detail)
        for category, n in sorted(o.dropped.items()):
            log.log(pdb, "WARNING", category, o.sli_id, o.het, n)
