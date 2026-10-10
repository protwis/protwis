"""Import Schrodinger Engine 2 receptor x peptide interactions, one transaction per structure.

    python manage.py import_schrodinger_peptides \\
        --data-dir DATA_DIR/structure_data/schrodinger/engine2 \\
        --anomaly-csv /runs/peptides/anomalies.csv --report-json /runs/peptides/report.json

The corpus is the database: every experimental structure with an anchor this
import serves (schrodinger_peptide.is_in_scope). Each must have a directory with
a peptide_map.tsv under --data-dir (build_schrodinger_peptide_maps); one that
has none fails rather than keep its old rows. --pdb / --pdb-list narrow the
corpus.

Each structure is imported in its own transaction and replaces, for each of
its "pep" anchors, the anchor's RFI rows and its peptide pairs. An anchor whose
items the run failed (preparation_failed, compute_failed, timed_out, crashed)
gets no rows and is reported as a WARNING (anchor_cleared); the structure and
the run go on. A structure that fails (its map could not be built, an item
record is unreadable or carries another contract version, or any unexpected
error) is rolled back and reported; the others are unaffected, and the command
exits non-zero after all were attempted. --dry-run runs every structure and
rolls each one back.

The anomaly CSV is written outside the transactions and flushed per row, so it
survives any rollback. Point it, and the report, at a fresh directory per run.
"""

import collections
import json
import os

from django.core.management.base import BaseCommand, CommandError
from django.db import transaction

from build.management.commands.import_schrodinger_interactions import (
    AnomalyLog,
    _Rollback,
    _make_room_for,
)
from interaction import schrodinger_import as si
from interaction import schrodinger_peptide as sp
from interaction.models import (
    ResidueFragmentInteractionType,
    StructureLigandInteraction,
)
from structure.models import Structure


class Command(BaseCommand):
    help = "Import Schrodinger Engine 2 receptor x peptide interactions (one transaction per structure)."

    def add_arguments(self, parser):
        parser.add_argument(
            "--data-dir",
            required=True,
            help="Engine 2 tree: {data_dir}/{PDB}/plan.json and peptide_map.tsv.",
        )
        parser.add_argument(
            "--pdb", action="append", default=[], help="PDB code to import; repeatable."
        )
        parser.add_argument(
            "--pdb-list",
            default=None,
            help="File with one PDB code per line (# comments allowed).",
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
            "--dry-run",
            action="store_true",
            help="Run every structure and roll each one back.",
        )

    @staticmethod
    def corpus():
        """PDB codes of the experimental structures with an in-scope anchor."""
        out = set()
        for sli in StructureLigandInteraction.objects.select_related(
            "structure__pdb_code", "structure__structure_type", "ligand__ligand_type"
        ):
            if sli.structure is None or not sp.is_in_scope(sli):
                continue
            if sli.structure.structure_type.origin != si.STRUCTURE_ORIGIN:
                continue
            out.add(sli.structure.pdb_code.index.upper())
        return sorted(out)

    def _codes(self, options):
        codes = [c.strip().upper() for c in options["pdb"] if c.strip()]
        if options["pdb_list"]:
            with open(options["pdb_list"]) as fh:
                for line in fh:
                    line = line.split("#", 1)[0].strip()
                    if line:
                        codes.append(line.upper())
        return list(dict.fromkeys(codes)) or self.corpus()

    def handle(self, *args, **options):
        data_dir = options["data_dir"]
        if not os.path.isdir(data_dir):
            raise CommandError("--data-dir {!r} is not a directory".format(data_dir))
        present = set(
            ResidueFragmentInteractionType.objects.values_list("slug", flat=True)
        )
        missing = sorted(si.required_slugs() - present)
        if missing:
            raise CommandError(
                "interaction types missing from the database: {} (see "
                "import_schrodinger_interactions for where each comes "
                "from)".format(", ".join(missing))
            )
        dirs = {
            name.upper(): name
            for name in os.listdir(data_dir)
            if os.path.isdir(os.path.join(data_dir, name))
        }
        codes = self._codes(options)
        log = AnomalyLog(options["anomaly_csv"])
        report, failed = [], []
        provenance = collections.defaultdict(collections.Counter)
        totals = collections.Counter()
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
                if structure is None:
                    log.log(pdb, "WARNING", "structure_not_in_db")
                    entry["status"] = "structure_not_in_db"
                    continue
                path = os.path.join(data_dir, dirs.get(pdb, pdb), sp.MAP_NAME)
                try:
                    if dirs.get(pdb, pdb) != pdb:
                        # The item paths are built from the upper-case code.
                        raise sp.MapMismatch(
                            "directory {!r} under --data-dir is not the upper-case "
                            "PDB code {}".format(dirs[pdb], pdb)
                        )
                    if not os.path.isfile(path):
                        raise sp.MapMismatch(
                            "no {} under {}".format(
                                sp.MAP_NAME, os.path.join(data_dir, dirs.get(pdb, pdb))
                            )
                        )
                    named, receptor, rows, prov, header = sp.load_peptide_map(path)
                    if named != pdb:
                        raise sp.MapMismatch("{} names {}".format(path, named))
                    for key, value in prov.items():
                        provenance[key][value] += 1
                    with transaction.atomic():
                        outcomes, cleanup = sp.import_structure(
                            structure, data_dir, receptor, rows, header
                        )
                        if options["dry_run"]:
                            raise _Rollback()
                except _Rollback:
                    pass
                except Exception as exc:
                    message = "{}: {}".format(type(exc).__name__, exc)[:300]
                    log.log(pdb, "ERROR", type(exc).__name__, detail=message)
                    entry["status"] = "failed"
                    entry["error"] = message
                    failed.append(pdb)
                    continue
                entry["status"] = "rolled_back" if options["dry_run"] else "imported"
                entry["receptor_segment"] = receptor["segment"]
                entry["cleanup"] = dict(cleanup)
                totals.update(cleanup)
                entry["anchors"] = []
                for o in outcomes:
                    self._log_outcome(log, pdb, o)
                    entry["anchors"].append(
                        {
                            "sli_id": o.sli_id,
                            "chain": o.chain,
                            "mode": o.mode,
                            "items": o.items,
                            "notes": o.notes,
                            "rfi_deleted": o.rfi_deleted,
                            "rfi_written": o.rfi_written,
                            "rfi_counts": dict(o.rfi_counts),
                            "rfi_dropped": dict(o.rfi_dropped),
                            "fragments_created": o.fragments_created,
                            "pairs_deleted": o.pairs_deleted,
                            "interactions_deleted": o.interactions_deleted,
                            "pairs_written": o.pairs_written,
                            "interactions_written": o.interactions_written,
                            "pair_counts": dict(o.pair_counts),
                            "pair_dropped": dict(o.pair_dropped),
                            "complex_file": o.complex_file,
                        }
                    )
                    totals["complex_file_" + o.complex_file] += 1
                    totals["anchors"] += 1
                    totals["cleared"] += o.mode in ("cleared", "no_product")
                    for key in (
                        "rfi_deleted",
                        "rfi_written",
                        "pairs_deleted",
                        "pairs_written",
                        "interactions_deleted",
                        "interactions_written",
                    ):
                        totals[key] += getattr(o, key)
        finally:
            log.close()
            if options["report_json"]:
                _make_room_for(options["report_json"])
                with open(options["report_json"], "w") as fh:
                    json.dump(
                        {
                            "dry_run": options["dry_run"],
                            "totals": dict(totals),
                            "provenance": {k: dict(v) for k, v in provenance.items()},
                            "failed": failed,
                            "structures": report,
                        },
                        fh,
                        indent=1,
                    )
        for key, values in sorted(provenance.items()):
            if len(values) > 1:
                self.stdout.write(
                    "{}: {} distinct values across the tree".format(key, len(values))
                )
        self.stdout.write(
            "{} structures, {} anchors ({} cleared); RFI deleted {} written {}; peptide pairs "
            "deleted {} written {}; peptide interactions deleted {} written {}; 3D files {}; "
            "orphan fragments deleted {}{}; anomalies INFO={} WARNING={} ERROR={}".format(
                len(codes),
                totals["anchors"],
                totals["cleared"],
                totals["rfi_deleted"],
                totals["rfi_written"],
                totals["pairs_deleted"],
                totals["pairs_written"],
                totals["interactions_deleted"],
                totals["interactions_written"],
                {
                    k[len("complex_file_") :]: v
                    for k, v in sorted(totals.items())
                    if k.startswith("complex_file_")
                },
                totals["fragments_deleted"],
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
        if o.mode in ("cleared", "no_product"):
            log.log(
                pdb,
                "WARNING",
                "anchor_cleared",
                o.sli_id,
                o.chain,
                o.rfi_deleted + o.pairs_deleted,
                detail="no product rows; {} RFI rows and {} peptide pairs deleted; {}".format(
                    o.rfi_deleted, o.pairs_deleted, "; ".join(o.notes)
                ),
            )
            return
        if o.complex_file == "ligand_not_found":
            log.log(
                pdb,
                "WARNING",
                "complex_file_missing",
                o.sli_id,
                o.chain,
                detail="rows written but the chain was not found in the stored structure "
                "text; the anchor has no 3D file",
            )
        if o.rfi_counts["rows_in"] == 0:
            log.log(
                pdb,
                "INFO",
                "product_has_zero_rows",
                o.sli_id,
                o.chain,
                detail=",".join(o.items),
            )
        for category, level in (
            ("excluded_family", "INFO"),
            ("nonstandard_residue", "INFO"),
            ("duplicate", "INFO"),
            ("other_chain", "WARNING"),
        ):
            if o.rfi_counts[category]:
                log.log(
                    pdb,
                    level,
                    "rfi_" + category,
                    o.sli_id,
                    o.chain,
                    o.rfi_counts[category],
                )
        for category in (
            "not_in_peptide_tables",
            "nonstandard_residue",
            "other_chain",
            "insertion_code_atoms",
        ):
            if o.pair_counts[category]:
                log.log(
                    pdb,
                    "INFO",
                    "pairs_" + category,
                    o.sli_id,
                    o.chain,
                    o.pair_counts[category],
                )
        for category, n in sorted(o.rfi_dropped.items()):
            log.log(pdb, "WARNING", "rfi_" + category, o.sli_id, o.chain, n)
        for category, n in sorted(o.pair_dropped.items()):
            log.log(pdb, "WARNING", "pairs_" + category, o.sli_id, o.chain, n)
        for note in o.notes:
            log.log(pdb, "INFO", "note", o.sli_id, o.chain, detail=note)
