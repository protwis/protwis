from build.management.commands.base_build import Command as BaseBuild

from django.core.management.base import CommandError
from django.db import connection, transaction
from django.db.models import Q

import logging
import math
import time


class Command(BaseBuild):
    help = (
        "Build structural distance similarity tables for representative human GPCR structures "
        "(active/inactive), using the same distance definition as the structure clustering page."
    )

    logger = logging.getLogger(__name__)

    @staticmethod
    def _format_elapsed(seconds):
        seconds = int(round(seconds))
        h = seconds // 3600
        m = (seconds % 3600) // 60
        s = seconds % 60
        return f"{h} hours {m} mins {s} secs"

    def add_arguments(self, parser):
        super(Command, self).add_arguments(parser=parser)
        parser.add_argument(
            '--state',
            choices=['active', 'inactive', 'both'],
            default='both',
            help='Which structure state(s) to build.',
        )
        parser.add_argument(
            '--batch-size',
            type=int,
            default=5000,
            help='Bulk insert batch size.',
        )
        parser.add_argument(
            '--verbose',
            action='store_true',
            default=False,
            help='Print progress to stdout.',
        )
        parser.add_argument(
            '--dry-run',
            action='store_true',
            default=False,
            help='Compute counts but do not write to the database.',
        )

    @staticmethod
    def _truncate_table(model):
        table = model._meta.db_table
        with connection.cursor() as cursor:
            if connection.vendor == 'postgresql':
                cursor.execute(f'TRUNCATE TABLE {connection.ops.quote_name(table)} RESTART IDENTITY')
            else:
                model.objects.all().delete()

    @staticmethod
    def _delete_state(model, state_slug):
        """
        Delete rows for a single state ('active'/'inactive') without affecting the other.
        """
        model.objects.filter(state__slug=state_slug).delete()

    @staticmethod
    def _expected_pairs(n):
        return (n * (n - 1)) // 2

    def _select_structures(self, state_slug, limit=None):
        """
        Representative, annotated, non-AF structures for human canonical receptors only.
        """
        from structure.models import Structure

        qs = (
            Structure.objects
            .exclude(structure_type__slug__startswith='af-')
            .filter(representative=True, annotated=True, state__slug=state_slug)
            .select_related(
                'pdb_code',
                'protein_conformation__protein',
                'protein_conformation__protein__parent',
                'protein_conformation__protein__parent__species',
                'protein_conformation__protein__parent__family',
            )
        )

        # Canonical receptor protein (parent) must be human
        qs = qs.filter(protein_conformation__protein__parent__isnull=False)
        qs = qs.filter(protein_conformation__protein__parent__species_id=1)
        qs = qs.filter(protein_conformation__protein__parent__parent_id__isnull=True)
        qs = qs.filter(protein_conformation__protein__parent__accession__isnull=False)

        # GPCR families only (exclude unclassified/artificial buckets)
        qs = qs.filter(protein_conformation__protein__parent__family__slug__startswith='0')

        if limit:
            qs = qs[:limit]

        return list(qs)

    @staticmethod
    def _canonical_protein_id(structure_obj):
        prot = structure_obj.protein_conformation.protein
        parent = getattr(prot, 'parent', None)
        return (parent.id if parent is not None else prot.id)

    def _build_one_state(self, state_slug, model, *, batch_size, verbose, dry_run, test):
        from contactnetwork.distances import Distances
        from protein.models import ProteinState

        t0 = time.time()

        limit = 10 if test else None
        structures = self._select_structures(state_slug, limit=limit)

        pdb_to_ids = {}
        for s in structures:
            if not s.pdb_code or not s.pdb_code.index:
                continue
            pdb = s.pdb_code.index.upper()
            if pdb in pdb_to_ids:
                # Should not happen, but keep first and warn.
                continue
            pdb_to_ids[pdb] = {
                'structure_id': s.id,
                'protein_id': self._canonical_protein_id(s),
            }

        pdbs = sorted(pdb_to_ids.keys())
        n_input = len(pdbs)

        if verbose:
            print(f"[{state_slug}] selected structures: {n_input}")

        if n_input < 2:
            self.logger.warning("[%s] Not enough structures (%s) to build distances.", state_slug, n_input)
            return 0

        # Compute distance matrix (ordering may change inside Distances)
        dis = Distances()
        dis.load_pdbs(pdbs)
        dm_raw = dis.get_distance_matrix(normalize=False)
        dm_norm = dis.get_distance_matrix(normalize=True)
        ordered_pdbs_all = [p.upper() for p in dis.pdbs]

        # Safety: if Distances reorders/drops any pdbs, keep matrix aligned.
        keep_indices = [idx for idx, pdb in enumerate(ordered_pdbs_all) if pdb in pdb_to_ids]
        if len(keep_indices) != len(ordered_pdbs_all):
            import numpy as np
            dm_raw = dm_raw[np.ix_(keep_indices, keep_indices)]
            dm_norm = dm_norm[np.ix_(keep_indices, keep_indices)]
        ordered_pdbs = [ordered_pdbs_all[idx] for idx in keep_indices]
        n = len(ordered_pdbs)

        expected = self._expected_pairs(n)
        if verbose:
            print(f"[{state_slug}] matrix size: {n}x{n} (pairs={expected}) (raw+normalized)")

        if dry_run:
            return expected

        # Reset existing rows for this state only
        self._delete_state(model, state_slug)

        state_obj = ProteinState.objects.only('id').get(slug=state_slug)

        buffer = []
        total = 0

        def _flush():
            nonlocal total
            if not buffer:
                return
            with transaction.atomic():
                model.objects.bulk_create(buffer, batch_size=batch_size)
            total += len(buffer)
            buffer.clear()
            if verbose:
                print(f"[{state_slug}] inserted rows: {total}")

        # Unfold upper triangle -> canonical (structure_ref_id < structure_target_id)
        for i in range(n):
            pdb_i = ordered_pdbs[i]
            a = pdb_to_ids[pdb_i]
            for j in range(i + 1, n):
                pdb_j = ordered_pdbs[j]
                b = pdb_to_ids[pdb_j]

                # Distances are symmetric; store a canonical direction by structure id
                dist_raw = float(dm_raw[i, j])
                dist_norm = float(dm_norm[i, j])
                if (
                    math.isnan(dist_raw) or math.isinf(dist_raw) or
                    math.isnan(dist_norm) or math.isinf(dist_norm)
                ):
                    continue

                if a['structure_id'] < b['structure_id']:
                    ref, tgt = a, b
                else:
                    ref, tgt = b, a

                buffer.append(
                    model(
                        state_id=state_obj.id,
                        structure_ref_id=ref['structure_id'],
                        structure_target_id=tgt['structure_id'],
                        protein_ref_id=ref['protein_id'],
                        protein_target_id=tgt['protein_id'],
                        distance=dist_raw,
                        distance_normalized=dist_norm,
                    )
                )

                if len(buffer) >= batch_size:
                    _flush()

        _flush()

        t1 = time.time()
        self.logger.info("[%s] Inserted %s rows in %s", state_slug, total, self._format_elapsed(t1 - t0))
        return total

    def handle(self, *args, **options):
        try:
            from classification.models import StructureSimilarity
        except ImportError as e:
            raise CommandError(
                "StructureSimilarity models not available. Did you migrate the classification app?"
            ) from e

        state_opt = options['state']
        batch_size = int(options['batch_size'])
        verbose = bool(options['verbose'])
        dry_run = bool(options['dry_run'])
        test = bool(options.get('test'))

        start = time.time()

        targets = []
        if state_opt in ('active', 'both'):
            targets.append('active')
        if state_opt in ('inactive', 'both'):
            targets.append('inactive')

        if not targets:
            raise CommandError("No states selected.")

        # If rebuilding both, it is safe and faster to truncate once.
        if not dry_run and state_opt == 'both':
            self._truncate_table(StructureSimilarity)

        for state_slug in targets:
            self._build_one_state(
                state_slug,
                StructureSimilarity,
                batch_size=batch_size,
                verbose=verbose,
                dry_run=dry_run,
                test=test,
            )

        end = time.time()
        if verbose:
            print("Execution time:", self._format_elapsed(end - start))

