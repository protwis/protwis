#Django Framework imports
from django.db.models import Count, Q, IntegerField, F, Value, Subquery, OuterRef,  Min, DecimalField
from django.db.models.functions import Coalesce, Replace

#GPCRdb imports
from structure.models import Structure, StructureExtraProteins
from interaction.models import StructureLigandInteraction
from protein.models import Protein, Gene


######################################
# GPCR Structure Coverage Statistics #
######################################

STIMULATORY_ROLES  = ['Agonist', 'Partial agonist', 'PAM']
INHIBITORY_ROLES   = ['Allosteric antagonist', 'Antagonist', 'Inverse agonist', 'NAM']

# ── subquery: count structures per primary/parent protein entity ──────────
def _struct_count(state_slug=None, extra_q=None):
    qs = Structure.objects.filter(
        ~Q(structure_type__slug__startswith='af-'),
        protein_conformation__protein__parent=OuterRef('pk'),
    )
    if state_slug:
        qs = qs.filter(state__slug=state_slug)
    if extra_q:
        qs = qs.filter(extra_q)
    return (
        qs.values('protein_conformation__protein__parent')
            .annotate(c=Count('id'))
            .values('c')
    )

# ── subquery: best (min) resolution per state ─────────────────────
def _best_resolution(state_slug=None):
    qs = Structure.objects.filter(
        ~Q(structure_type__slug__startswith='af-'),
        protein_conformation__protein__parent=OuterRef('pk'),
        resolution__isnull=False,
    )
    if state_slug:
        qs = qs.filter(state__slug=state_slug)
    return (
        qs.values('protein_conformation__protein__parent')
            .annotate(best=Min('resolution'))
            .values('best')
    )

# ── subquery: count structures with a given extra-protein type ────
def _transducer_count(category, family_slug_prefix=None):
    qs = StructureExtraProteins.objects.filter(
        ~Q(structure__structure_type__slug__startswith='af-'),
        structure__protein_conformation__protein__parent=OuterRef('pk'),
        category=category,
    )
    if family_slug_prefix:
        qs = qs.filter(wt_protein__family__slug__startswith=family_slug_prefix)
    return (
        qs.values('structure__protein_conformation__protein__parent')
            .annotate(c=Count('structure', distinct=True))
            .values('c')
    )

# ── subquery: count structures with GRK stabilising agent ─────────
def _grk_count():
    return (
        Structure.objects.filter(
            ~Q(structure_type__slug__startswith='af-'),
            protein_conformation__protein__parent=OuterRef('pk'),
            stabilizing_agents__name__icontains='GRK',
        )
            .values('protein_conformation__protein__parent')
            .annotate(c=Count('id', distinct=True))
            .values('c')
    )

# ── subquery: count annotated ligands of given types ──────────────
def _ligand_type_count(type_slugs):
    return (
        StructureLigandInteraction.objects.filter(
            ~Q(structure__structure_type__slug__startswith='af-'),
            structure__protein_conformation__protein__parent=OuterRef('pk'),
            annotated=True,
            ligand__ligand_type__slug__in=type_slugs,
        )
            .values('structure__protein_conformation__protein__parent')
            .annotate(c=Count('id'))
            .values('c')
    )

# ── subquery: count annotated ligand roles by name ────────────────
def _role_count(role_names, invert=False):

    qry = StructureLigandInteraction.objects.filter(
        structure__protein_conformation__protein__parent=OuterRef('pk'),
        annotated=True)

    if not invert:
        qry = qry.filter(ligand_role__name__in=role_names)
    else:
        qry = qry.filter(~Q(ligand_role__name__in=role_names))

    return (qry
            .filter(~Q(structure__structure_type__slug__startswith='af-'))
            .values('structure__protein_conformation__protein__parent')
            .annotate(c=Count('id'))
            .values('c')
    )

def get_primary_gene_name_subquery():
    return Subquery(
        Gene.objects.filter(proteins=OuterRef('pk'), position=0).order_by('entrez_id').values('name')[:1]
    )

def get_primary_gene_entrezid_subquery():
    return Subquery(
        Gene.objects.filter(proteins=OuterRef('pk'), position=0).order_by('entrez_id').values('entrez_id')[:1]
    )

#------------#
# Main Query #
#------------#
def GpcrStructureCoverageStatisticsQuery():
    """Query to generate a summary of structure coverage statistics for GPCRs

    Generates summary statistics for each receptor including counts of structures
    in different states, best resolution, transducer counts, ligand types, and modality counts.
    Returns a queryset of Protein objects annotated with these statistics.
    """
    queryset = Protein.objects.filter(
        #Q(source__name='SWISSPROT') | Q(source__name='TREMBL'),
        parent__isnull=True,
        accession__isnull=False
    ).select_related(
        'family__parent__parent__parent',
        'species',
    ).annotate(
        ## Basic protein identity
        uniprot_entry_name = F('entry_name'),
        uniprot_accession = F('accession'),
        protein_name = Replace(F('name'), Value(' receptor'), Value('')),
        receptor_family = F('family__parent__name'),
        receptor_class = F('family__parent__parent__parent__name'),
        gene_name = get_primary_gene_name_subquery(),
        gene_entrez_id = get_primary_gene_entrezid_subquery(),
        species_common_name = F('species__common_name'),

        # ── structure counts by state ──────────────────────────────
        st_inact_c=Coalesce(Subquery(
            _struct_count('inactive'), output_field=IntegerField()), 0),
        st_interm_c=Coalesce(Subquery(
            _struct_count('intermediate'), output_field=IntegerField()), 0),
        st_act_c=Coalesce(Subquery(
            _struct_count('active'), output_field=IntegerField()), 0),

        # ── best (lowest) resolution by state ─────────────────────
        best_res_inact=Subquery(
            _best_resolution('inactive'), output_field=DecimalField()),
        best_res_interm=Subquery(
            _best_resolution('intermediate'), output_field=DecimalField()),
        best_res_act=Subquery(
            _best_resolution('active'), output_field=DecimalField()),

        # ── G-protein family slug prefixes ──────────────────────────────
        # 100_001_001 = Gs, 100_001_002 = Gi/o,
        # 100_001_003 = Gq/11, 100_001_004 = G12/13
        # ── transducer complex counts (structural) ─────────────────
        transd_gio_c=Coalesce(Subquery(
            _transducer_count('G alpha', '100_001_002'), output_field=IntegerField()), 0),
        transd_gq11_c=Coalesce(Subquery(
            _transducer_count('G alpha', '100_001_003'), output_field=IntegerField()), 0),
        transd_g1213_c=Coalesce(Subquery(
            _transducer_count('G alpha', '100_001_004'), output_field=IntegerField()), 0),
        transd_arrest_c=Coalesce(Subquery(
            _transducer_count('Arrestin'), output_field=IntegerField()), 0),
        transd_grks_c=Coalesce(Subquery(
            _grk_count(), output_field=IntegerField()), 0),

        # ── ligand type counts ─────────────────────────────────────
        ligand_types_smmol_c=Coalesce(Subquery(
            _ligand_type_count(['small-molecule']), output_field=IntegerField()), 0),
        ligand_types_pep_c=Coalesce(Subquery(
            _ligand_type_count(['peptide', 'protein']), output_field=IntegerField()), 0),

        # ── modality group counts ──────────────────────────────────
        limodgrp_stim_c=Coalesce(Subquery(
            _role_count(STIMULATORY_ROLES), output_field=IntegerField()), 0),
        limodgrp_inhib_c=Coalesce(Subquery(
            _role_count(INHIBITORY_ROLES), output_field=IntegerField()), 0),
        limodgrp_unk_c=Coalesce(Subquery(
            _role_count(['Unknown']), output_field=IntegerField()), 0),
        limodgrp_other_c=Coalesce(Subquery(
            _role_count(STIMULATORY_ROLES + INHIBITORY_ROLES + ['Unknown'], invert=True), output_field=IntegerField()), 0),

        # ── individual modality counts ─────────────────────────────
        limod_apo_c=Coalesce(Subquery(
            _role_count(['Apo']), output_field=IntegerField()), 0),
        limod_ag_c=Coalesce(Subquery(
            _role_count(['Agonist']), output_field=IntegerField()), 0),
        limod_ag_partial_c=Coalesce(Subquery(
            _role_count(['Partial agonist']), output_field=IntegerField()), 0),
        limod_allo_antag_c=Coalesce(Subquery(
            _role_count(['Allosteric antagonist']), output_field=IntegerField()), 0),
        limod_antag_c=Coalesce(Subquery(
            _role_count(['Antagonist']), output_field=IntegerField()), 0),
        limod_inv_ag_c=Coalesce(Subquery(
            _role_count(['Inverse agonist']), output_field=IntegerField()), 0),
        limod_nam_c=Coalesce(Subquery(
            _role_count(['NAM']), output_field=IntegerField()), 0),
        limod_pam_c=Coalesce(Subquery(
            _role_count(['PAM']), output_field=IntegerField()), 0),
        limod_unk_c=Coalesce(Subquery(
            _role_count(['Unknown']), output_field=IntegerField()), 0),
    ).values(
            "uniprot_entry_name",
            "uniprot_accession",
            "protein_name",
            "receptor_family",
            "receptor_class",
            "gene_name",
            "gene_entrez_id",
            "species_common_name",
            "st_inact_c",
            "st_interm_c",
            "st_act_c",
            "best_res_inact",
            "best_res_interm",
            "best_res_act",
            "transd_gio_c",
            "transd_gq11_c",
            "transd_g1213_c",
            "transd_arrest_c",
            "transd_grks_c",
            "ligand_types_smmol_c",
            "ligand_types_pep_c",
            "limodgrp_stim_c",
            "limodgrp_inhib_c",
            "limodgrp_unk_c",
            "limod_apo_c",
            "limod_ag_c",
            "limod_ag_partial_c",
            "limod_allo_antag_c",
            "limod_antag_c",
            "limod_inv_ag_c",
            "limod_nam_c",
            "limod_pam_c",
            "limod_unk_c",
            ).filter(
            # Filter to only include proteins with at least one structure or ligand interaction
            Q(best_res_inact__isnull = False) |
            Q(best_res_interm__isnull = False) |
            Q(best_res_act__isnull = False) |
            Q(st_inact_c__gt=0) |
            Q(st_interm_c__gt=0) |
            Q(st_act_c__gt=0) |
            Q(transd_gio_c__gt=0) |
            Q(transd_gq11_c__gt=0) |
            Q(transd_g1213_c__gt=0) |
            Q(transd_arrest_c__gt=0) |
            Q(transd_grks_c__gt=0) |
            Q(ligand_types_smmol_c__gt=0) |
            Q(ligand_types_pep_c__gt=0) |
            Q(limodgrp_stim_c__gt=0) |
            Q(limodgrp_inhib_c__gt=0) |
            Q(limodgrp_unk_c__gt=0) |
            Q(limod_apo_c__gt=0) |
            Q(limod_ag_c__gt=0) |
            Q(limod_ag_partial_c__gt=0) |
            Q(limod_allo_antag_c__gt=0) |
            Q(limod_antag_c__gt=0) |
            Q(limod_inv_ag_c__gt=0) |
            Q(limod_nam_c__gt=0) |
            Q(limod_pam_c__gt=0) |
            Q(limod_unk_c__gt=0)
        )

    return queryset