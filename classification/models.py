from django.contrib.postgres.fields import JSONField
from django.db import models, connection


# Add models here as you copy/implement classification functionality.


class CustomReceptorSimilarityManager(models.Manager):
    def truncate_table(self):
        with connection.cursor() as cursor:
            cursor.execute(f'TRUNCATE TABLE "{self.model._meta.db_table}" CASCADE')


class ReceptorSimilarity(models.Model):
    protein_ref = models.ForeignKey(
        'protein.Protein',
        on_delete=models.CASCADE,
        related_name='receptor_similarity_as_ref',
        db_column='ref',
        db_index=True,
    )
    protein_target = models.ForeignKey(
        'protein.Protein',
        on_delete=models.CASCADE,
        related_name='receptor_similarity_as_target',
        db_column='target',
        db_index=True,
    )
    identity = models.PositiveSmallIntegerField()
    similarity = models.PositiveSmallIntegerField()

    # Top-level class FKs (nullable for backfill)
    ref_class = models.ForeignKey(
        'protein.ProteinFamily',
        null=True,
        blank=True,
        on_delete=models.CASCADE,
        related_name='sim_as_ref_class',
        db_column='ref_class',
        db_index=True,
    )
    target_class = models.ForeignKey(
        'protein.ProteinFamily',
        null=True,
        blank=True,
        on_delete=models.CASCADE,
        related_name='sim_as_target_class',
        db_column='target_class',
        db_index=True,
    )

    objects = models.Manager()
    custom_objects = CustomReceptorSimilarityManager()

    class Meta:
        db_table = 'classification_receptorsimilarity'
        constraints = [
            models.UniqueConstraint(
                fields=['protein_ref', 'protein_target'],
                name='crs_uniq_pair',
            ),
            models.CheckConstraint(
                check=~models.Q(protein_ref=models.F('protein_target')),
                name='crs_ref_ne_target',
            ),
        ]
        indexes = [
            models.Index(fields=['ref_class', 'target_class', 'identity'], name='crs_cls_id_idx'),
            models.Index(fields=['ref_class', 'target_class', 'similarity'], name='crs_cls_sim_idx'),
            models.Index(fields=['target_class', 'ref_class', 'identity'], name='crs_cls_id_rev_idx'),
            models.Index(fields=['target_class', 'ref_class', 'similarity'], name='crs_cls_sim_rev_idx'),
        ]

    def __str__(self):
        return f'{self.protein_ref} vs {self.protein_target}: sim={self.similarity} id={self.identity}'


class _StructureSimilarityBase(models.Model):
    """
    Pairwise structural distance between two representative structures.

    Stored in one canonical direction only (ref != target) and intended to be
    populated by a build management command.
    """

    structure_ref = models.ForeignKey(
        'structure.Structure',
        on_delete=models.CASCADE,
        related_name='+',
        db_index=True,
    )
    structure_target = models.ForeignKey(
        'structure.Structure',
        on_delete=models.CASCADE,
        related_name='+',
        db_index=True,
    )

    # Raw (unnormalized) distance and normalized distance (as used by clustering default)
    distance = models.FloatField()
    distance_normalized = models.FloatField(default=0.0)

    # Canonical (parent) receptor protein IDs for convenience
    protein_ref = models.ForeignKey(
        'protein.Protein',
        on_delete=models.CASCADE,
        related_name='+',
        db_index=True,
    )
    protein_target = models.ForeignKey(
        'protein.Protein',
        on_delete=models.CASCADE,
        related_name='+',
        db_index=True,
    )

    class Meta:
        abstract = True

    def __str__(self):
        return f'{self.structure_ref_id} vs {self.structure_target_id}: d={self.distance}'


class StructureSimilarity(_StructureSimilarityBase):
    """
    Pairwise structural distance between representative structures, for a given state.

    `state` indicates whether the structure pair belongs to the active or inactive ensemble.
    """

    # Use the same ProteinState as Structure.state (slug 'active'/'inactive')
    state = models.ForeignKey(
        'protein.ProteinState',
        on_delete=models.CASCADE,
        related_name='+',
        db_index=True,
    )

    class Meta(_StructureSimilarityBase.Meta):
        db_table = 'classification_structuresimilarity'
        constraints = [
            models.UniqueConstraint(
                fields=['state', 'structure_ref', 'structure_target'],
                name='css_uniq_state_pair',
            ),
            models.CheckConstraint(
                check=~models.Q(structure_ref=models.F('structure_target')),
                name='css_ref_ne_target',
            ),
        ]
        indexes = [
            models.Index(fields=['state', 'protein_ref', 'protein_target'], name='css_state_prot_pair_idx'),
            models.Index(fields=['state', 'structure_ref', 'structure_target'], name='css_state_struct_pair_idx'),
            models.Index(fields=['distance'], name='css_dist_idx'),
            models.Index(fields=['distance_normalized'], name='css_distn_idx'),
        ]


class ClusterCoord(models.Model):
    """
    Persisted 2D coordinates for clustering/embedding plots (e.g. StructureSim).

    One row per (protein, dataset_type, plot_type, group_key).
    """

    DATASET_SEQUENCE = 'sequence'
    DATASET_STRUCTURE_ACTIVE = 'structure_active'
    DATASET_STRUCTURE_INACTIVE = 'structure_inactive'

    PLOT_TSNE = 'tsne'

    DATASET_CHOICES = (
        (DATASET_SEQUENCE, 'Sequence'),
        (DATASET_STRUCTURE_ACTIVE, 'Structure (active)'),
        (DATASET_STRUCTURE_INACTIVE, 'Structure (inactive)'),
    )

    PLOT_CHOICES = (
        (PLOT_TSNE, 't-SNE'),
    )

    protein = models.ForeignKey(
        'protein.Protein',
        on_delete=models.CASCADE,
        related_name='+',
        db_index=True,
    )
    dataset_type = models.CharField(max_length=32, choices=DATASET_CHOICES, db_index=True)
    plot_type = models.CharField(max_length=16, choices=PLOT_CHOICES, default=PLOT_TSNE, db_index=True)
    group_key = models.CharField(max_length=32, default='global', db_index=True)

    x = models.FloatField()
    y = models.FloatField()

    class Meta:
        db_table = 'classification_clustercoord'
        constraints = [
            models.UniqueConstraint(
                fields=['protein', 'dataset_type', 'plot_type', 'group_key'],
                name='ccoord_uniq_prot_ds_plot_grp',
            ),
        ]
        indexes = [
            models.Index(fields=['dataset_type', 'plot_type', 'group_key'], name='ccoord_ds_plot_grp_idx'),
            models.Index(fields=['protein', 'dataset_type', 'group_key'], name='ccoord_prot_ds_grp_idx'),
        ]

    def __str__(self):
        return (
            f'{self.protein_id} {self.dataset_type}/{self.plot_type}/'
            f'{self.group_key}: ({self.x:.3f}, {self.y:.3f})'
        )


class TreeNetwork(models.Model):
    """
    Persisted render-ready family tree payloads.

    One row per receptor-family visualization group, including synthetic orphan
    splits such as Class A orphans.
    """

    group_key = models.CharField(max_length=64, unique=True)
    family = models.ForeignKey(
        'protein.ProteinFamily',
        null=True,
        blank=True,
        on_delete=models.CASCADE,
        related_name='+',
        db_index=True,
    )
    class_family = models.ForeignKey(
        'protein.ProteinFamily',
        null=True,
        blank=True,
        on_delete=models.CASCADE,
        related_name='+',
        db_index=True,
    )
    display_name = models.CharField(max_length=200)
    protein_count = models.PositiveSmallIntegerField(default=0)
    tree_method = models.CharField(max_length=64, default='neighbor_joining')
    segment_source = models.CharField(max_length=64, default='generic_conserved')
    bootstrap = models.PositiveSmallIntegerField(default=0)
    branch_mode = models.CharField(max_length=32, default='regular')
    tree_newick = models.TextField(blank=True, default='')
    payload = JSONField(default=dict)
    build_version = models.CharField(max_length=32, default='v1')
    source_hash = models.CharField(max_length=40, blank=True, default='')
    created_at = models.DateTimeField(auto_now_add=True)
    updated_at = models.DateTimeField(auto_now=True)

    class Meta:
        db_table = 'classification_treenetwork'
        indexes = [
            models.Index(fields=['class_family', 'group_key'], name='ctn_class_group_idx'),
            models.Index(fields=['family', 'group_key'], name='ctn_family_group_idx'),
        ]

    def __str__(self):
        return f'{self.group_key}: {self.display_name} ({self.protein_count})'
