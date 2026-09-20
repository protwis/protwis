from build.management.commands.base_build import Command as BaseBuild

from django.core.cache import cache, caches
from django.core.management.base import CommandError
from django.db.models import Q
from django.db import transaction

import logging
import math
import time

import numpy as np
from sklearn.manifold import TSNE


try:
    cache_alignment = caches["alignments"]
except Exception:
    cache_alignment = cache


class Command(BaseBuild):
    help = (
        "Build persisted 2D coordinates for StructureSim plots and store them in "
        "classification_clustercoord (global + per-class sequence/structure datasets)."
    )

    logger = logging.getLogger(__name__)

    GLOBAL_GROUP_KEY = "global"
    CLASS_SLUG_BY_KEY = {
        "A": "001",
        "B1": "002",
        "B2": "003",
        "C": "004",
        "F": "006",
        "O1": "007",
        "O2": "008",
        "T2": "009",
        "V": "010",
        "U": "011",
    }

    @staticmethod
    def _format_elapsed(seconds):
        seconds = int(round(seconds))
        h = seconds // 3600
        m = (seconds % 3600) // 60
        s = seconds % 60
        return f"{h} hours {m} mins {s} secs"

    @classmethod
    def _class_group_key(cls, class_key):
        return f"class:{class_key}"

    @classmethod
    def _class_family_q(cls, class_key, prefix=""):
        slug = cls.CLASS_SLUG_BY_KEY.get(class_key)
        if not slug:
            return Q()
        return Q(**{prefix + "slug": slug})

    def add_arguments(self, parser):
        super(Command, self).add_arguments(parser=parser)
        parser.add_argument(
            "--state",
            choices=["active", "inactive", "both"],
            default="both",
            help="Which structure state(s) to build.",
        )
        parser.add_argument(
            "--batch-size",
            type=int,
            default=5000,
            help="Bulk insert batch size.",
        )
        parser.add_argument(
            "--verbose",
            action="store_true",
            default=False,
            help="Print progress to stdout.",
        )
        parser.add_argument(
            "--dry-run",
            action="store_true",
            default=False,
            help="Compute counts but do not write to the database.",
        )

    @staticmethod
    def _resolve_perplexity(n, perplexity=None):
        if n < 2:
            return 1.0
        if perplexity is None:
            perplexity = min(40.0, max(1.0, (n - 1) / 3.0))
        else:
            perplexity = float(perplexity)
        if perplexity >= (n - 1):
            perplexity = float(max(1, n - 2))
        return perplexity

    @classmethod
    def _compute_tsne(cls, D, perplexity=None):
        n = D.shape[0]
        if n < 2:
            return np.zeros((n, 2), dtype=float)
        perplexity = cls._resolve_perplexity(n, perplexity=perplexity)

        base_kwargs = dict(
            n_components=2,
            metric="precomputed",
            perplexity=perplexity,
            random_state=42,
            init="random",
        )
        try:
            tsne = TSNE(learning_rate="auto", square_distances=True, **base_kwargs)
        except TypeError:
            tsne = TSNE(learning_rate=200.0, **base_kwargs)
        return tsne.fit_transform(D)

    def _persist_coords(
        self,
        *,
        ClusterCoord,
        proteins,
        coords,
        dataset_type,
        group_key,
        batch_size,
        verbose,
        label,
        dry_run,
    ):
        if dry_run:
            return len(proteins)

        buffer = []
        total = 0

        def _flush():
            nonlocal total
            if not buffer:
                return
            with transaction.atomic():
                ClusterCoord.objects.bulk_create(buffer, batch_size=batch_size)
            total += len(buffer)
            buffer.clear()
            if verbose:
                print(f"[{label}] inserted rows: {total}")

        for i, protein in enumerate(proteins):
            buffer.append(
                ClusterCoord(
                    protein_id=protein.id,
                    dataset_type=dataset_type,
                    plot_type=ClusterCoord.PLOT_TSNE,
                    group_key=group_key,
                    x=float(coords[i, 0]),
                    y=float(coords[i, 1]),
                )
            )
            if len(buffer) >= batch_size:
                _flush()
        _flush()
        return total

    def _build_similarity_dataset_from_qs(
        self,
        *,
        label,
        dataset_type,
        qs,
        ClusterCoord,
        Protein,
        group_key,
        batch_size,
        verbose,
        dry_run,
    ):
        prot_ids = set()
        for ref_id, tgt_id, _sim in qs.iterator():
            prot_ids.add(ref_id)
            prot_ids.add(tgt_id)

        proteins = list(
            Protein.objects.filter(id__in=prot_ids, species_id=1)
            .order_by("entry_name")
            .only("id", "entry_name")
        )
        n = len(proteins)
        if verbose:
            print(f"[{label}] proteins: {n}")
        if n == 0:
            return 0

        idx = {protein.id: i for i, protein in enumerate(proteins)}
        D = np.full((n, n), 100.0, dtype=float)
        np.fill_diagonal(D, 0.0)

        for ref_id, tgt_id, sim in qs.iterator():
            i = idx.get(ref_id)
            j = idx.get(tgt_id)
            if i is None or j is None or i == j:
                continue
            try:
                dist = 100.0 - float(sim)
            except Exception:
                continue
            D[i, j] = dist
            D[j, i] = dist

        coords_tsne = self._compute_tsne(D)
        return self._persist_coords(
            ClusterCoord=ClusterCoord,
            proteins=proteins,
            coords=coords_tsne,
            dataset_type=dataset_type,
            group_key=group_key,
            batch_size=batch_size,
            verbose=verbose,
            label=label,
            dry_run=dry_run,
        )

    def _build_sequence(
        self,
        *,
        ClusterCoord,
        ReceptorSimilarity,
        Protein,
        ProteinFamily,
        batch_size,
        verbose,
        dry_run,
    ):
        counts = {}
        global_qs = (
            ReceptorSimilarity.objects.filter(
                protein_ref__species_id=1,
                protein_target__species_id=1,
            ).values_list("protein_ref_id", "protein_target_id", "similarity")
        )
        counts[self.GLOBAL_GROUP_KEY] = self._build_similarity_dataset_from_qs(
            label="sequence/global",
            dataset_type=ClusterCoord.DATASET_SEQUENCE,
            qs=global_qs,
            ClusterCoord=ClusterCoord,
            Protein=Protein,
            group_key=self.GLOBAL_GROUP_KEY,
            batch_size=batch_size,
            verbose=verbose,
            dry_run=dry_run,
        )

        for class_key in self.CLASS_SLUG_BY_KEY.keys():
            family_ids = list(
                ProteinFamily.objects
                .filter(self._class_family_q(class_key))
                .values_list("id", flat=True)
            )
            if not family_ids:
                counts[self._class_group_key(class_key)] = 0
                continue
            class_qs = (
                ReceptorSimilarity.objects.filter(
                    protein_ref__species_id=1,
                    protein_target__species_id=1,
                    ref_class_id__in=family_ids,
                    target_class_id__in=family_ids,
                ).values_list("protein_ref_id", "protein_target_id", "similarity")
            )
            group_key = self._class_group_key(class_key)
            counts[group_key] = self._build_similarity_dataset_from_qs(
                label=f"sequence/{class_key}",
                dataset_type=ClusterCoord.DATASET_SEQUENCE,
                qs=class_qs,
                ClusterCoord=ClusterCoord,
                Protein=Protein,
                group_key=group_key,
                batch_size=batch_size,
                verbose=verbose,
                dry_run=dry_run,
            )
        return counts

    def _build_structure_state(
        self,
        state_slug,
        *,
        ClusterCoord,
        StructureSimilarity,
        ProteinState,
        Protein,
        batch_size,
        verbose,
        dry_run,
        class_key=None,
    ):
        state_obj = ProteinState.objects.only("id").get(slug=state_slug)
        qs = StructureSimilarity.objects.filter(
            state_id=state_obj.id,
            protein_ref__species_id=1,
            protein_target__species_id=1,
        )
        if class_key:
            ref_class_q = self._class_family_q(
                class_key,
                prefix="protein_ref__family__parent__parent__parent__",
            )
            target_class_q = self._class_family_q(
                class_key,
                prefix="protein_target__family__parent__parent__parent__",
            )
            if not ref_class_q or not target_class_q:
                return 0
            qs = qs.filter(ref_class_q, target_class_q)

        value_qs = qs.values_list("protein_ref_id", "protein_target_id", "distance")

        prot_ids = set()
        max_dist = 0.0
        pairs = {}

        for a_id, b_id, distance in value_qs.iterator():
            if a_id == b_id:
                continue
            prot_ids.add(a_id)
            prot_ids.add(b_id)
            try:
                dist = float(distance)
            except Exception:
                continue
            if math.isnan(dist) or math.isinf(dist):
                continue
            if dist > max_dist:
                max_dist = dist
            key = (a_id, b_id) if a_id < b_id else (b_id, a_id)
            prev = pairs.get(key)
            pairs[key] = dist if (prev is None or dist < prev) else prev

        prot_qs = (
            Protein.objects.filter(id__in=prot_ids, species_id=1)
            .order_by("entry_name")
            .only("id", "entry_name")
        )
        if state_slug == "inactive":
            prot_qs = prot_qs.exclude(entry_name="ccr9_human")

        proteins = list(prot_qs)
        n = len(proteins)
        label_suffix = class_key or "global"
        if verbose:
            print(f"[{state_slug}/{label_suffix}] proteins: {n} unique pairs: {len(pairs)}")
        if n == 0:
            return 0

        if max_dist <= 0:
            max_dist = 1.0

        idx = {protein.id: i for i, protein in enumerate(proteins)}
        D = np.full((n, n), max_dist, dtype=float)
        np.fill_diagonal(D, 0.0)

        for (a_id, b_id), dist in pairs.items():
            i = idx.get(a_id)
            j = idx.get(b_id)
            if i is None or j is None or i == j:
                continue
            D[i, j] = dist
            D[j, i] = dist

        coords_tsne = self._compute_tsne(D)
        dataset_type = (
            ClusterCoord.DATASET_STRUCTURE_ACTIVE
            if state_slug == "active"
            else ClusterCoord.DATASET_STRUCTURE_INACTIVE
        )
        group_key = self._class_group_key(class_key) if class_key else self.GLOBAL_GROUP_KEY
        return self._persist_coords(
            ClusterCoord=ClusterCoord,
            proteins=proteins,
            coords=coords_tsne,
            dataset_type=dataset_type,
            group_key=group_key,
            batch_size=batch_size,
            verbose=verbose,
            label=f"{state_slug}/{label_suffix}",
            dry_run=dry_run,
        )

    def handle(self, *args, **options):
        try:
            from classification.models import ClusterCoord, ReceptorSimilarity, StructureSimilarity
        except ImportError as e:
            raise CommandError("Classification models not available. Did you migrate the classification app?") from e

        from protein.models import Protein, ProteinFamily, ProteinState

        state_opt = options["state"]
        batch_size = int(options["batch_size"])
        verbose = bool(options["verbose"])
        dry_run = bool(options["dry_run"])
        test = bool(options.get("test"))
        if test:
            raise CommandError("This command does not support --test; it writes summary coordinate tables.")

        t0 = time.time()

        if not dry_run:
            if verbose:
                print("[clustercoord] clearing classification_clustercoord ...")
            ClusterCoord.objects.all().delete()

        seq_counts = self._build_sequence(
            ClusterCoord=ClusterCoord,
            ReceptorSimilarity=ReceptorSimilarity,
            Protein=Protein,
            ProteinFamily=ProteinFamily,
            batch_size=batch_size,
            verbose=verbose,
            dry_run=dry_run,
        )

        targets = []
        if state_opt in ("active", "both"):
            targets.append("active")
        if state_opt in ("inactive", "both"):
            targets.append("inactive")

        struct_counts = {}
        for state_slug in targets:
            if not ProteinState.objects.filter(slug=state_slug).exists():
                message = (
                    f"[{state_slug}] ProteinState '{state_slug}' not found — skipping structure "
                    "build for this state (structure data likely not built yet)."
                )
                if verbose:
                    print(message)
                self.logger.warning(message)
                struct_counts[state_slug] = {}
                continue

            state_counts = {}
            state_counts[self.GLOBAL_GROUP_KEY] = self._build_structure_state(
                state_slug,
                ClusterCoord=ClusterCoord,
                StructureSimilarity=StructureSimilarity,
                ProteinState=ProteinState,
                Protein=Protein,
                batch_size=batch_size,
                verbose=verbose,
                dry_run=dry_run,
            )
            for class_key in self.CLASS_SLUG_BY_KEY.keys():
                group_key = self._class_group_key(class_key)
                state_counts[group_key] = self._build_structure_state(
                    state_slug,
                    ClusterCoord=ClusterCoord,
                    StructureSimilarity=StructureSimilarity,
                    ProteinState=ProteinState,
                    Protein=Protein,
                    batch_size=batch_size,
                    verbose=verbose,
                    dry_run=dry_run,
                    class_key=class_key,
                )
            struct_counts[state_slug] = state_counts

        if not dry_run:
            cache_alignment.delete("structuresim:payload:db:v5:clustercoord")

        t1 = time.time()
        self.logger.info(
            "Built ClusterCoord: seq(global)=%s active(global)=%s inactive(global)=%s in %s",
            seq_counts.get(self.GLOBAL_GROUP_KEY, 0),
            struct_counts.get("active", {}).get(self.GLOBAL_GROUP_KEY, 0),
            struct_counts.get("inactive", {}).get(self.GLOBAL_GROUP_KEY, 0),
            self._format_elapsed(t1 - t0),
        )

