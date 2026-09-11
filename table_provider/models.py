from django.db import models
from django.contrib.postgres.fields import JSONField

class GpcrStructureStatisticsTable(models.Model):
    """Model representing the structure statistics table for GPCRs generated from the data retrieved by GpcrStructureCoverageStatisticsQuery"""

    uniprot_entry_name = models.CharField(max_length=20, null=False)
    uniprot_accession = models.CharField(max_length=20, null=False)
    protein_name = models.CharField(max_length=100, null=False)
    receptor_family = models.CharField(max_length=100, null=True)
    receptor_class = models.CharField(max_length=100, null=True)
    gene_name = models.CharField(max_length=100, null=True)
    gene_entrez_id = models.CharField(max_length=100, null=True)
    species_common_name = models.CharField(max_length=100, null=True)

    # Structure counts by state
    st_inact_c = models.IntegerField(null=False)
    st_interm_c = models.IntegerField(null=False)
    st_act_c = models.IntegerField(null=False)

    # Best resolution by state
    best_res_inact = models.FloatField(null=True)
    best_res_interm = models.FloatField(null=True)
    best_res_act = models.FloatField(null=True)

    # Transducer / extra-protein counts
    transd_gio_c = models.IntegerField(null=False)
    transd_gq11_c = models.IntegerField(null=False)
    transd_g1213_c = models.IntegerField(null=False)
    transd_arrest_c = models.IntegerField(null=False)
    transd_grks_c = models.IntegerField(null=False)

    # Ligand type counts
    ligand_types_smmol_c = models.IntegerField(null=False)
    ligand_types_pep_c = models.IntegerField(null=False)

    # Modality group counts
    limodgrp_stim_c = models.IntegerField(null=False)
    limodgrp_inhib_c = models.IntegerField(null=False)
    limodgrp_unk_c = models.IntegerField(null=False)

    # Individual modality counts
    limod_apo_c = models.IntegerField(null=False)
    limod_ag_c = models.IntegerField(null=False)
    limod_ag_partial_c = models.IntegerField(null=False)
    limod_allo_antag_c = models.IntegerField(null=False)
    limod_antag_c = models.IntegerField(null=False)
    limod_inv_ag_c = models.IntegerField(null=False)
    limod_nam_c = models.IntegerField(null=False)
    limod_pam_c = models.IntegerField(null=False)
    limod_unk_c = models.IntegerField(null=False)

    def __str__(self):
        """Creates a string representation of the GpcrStructureStatisticsTable instance, including all relevant fields and their values.

        Returns:
            str: A string representation of the GpcrStructureStatisticsTable instance.
        """
        str = f'''
            "uniprot_entry_name" : {self.uniprot_entry_name},
            "protein_name" : {self.protein_name},
            "receptor_family" : {self.receptor_family},
            "receptor_class" : {self.receptor_class},
            "gene_name" : {self.gene_name},
            "gene_entrez_id" : {self.gene_entrez_id},
            "species_common_name" : {self.species_common_name},
            "structure_inactive_c" : {self.st_inact_c },
            "structure_intermediate_c" : {self.st_interm_c },
            "structure_active_c" : {self.st_act_c },
            "best_resolution_inactive" : {self.best_res_inact },
            "best_resolution_intermediate" : {self.best_res_interm },
            "best_resolution_active" : {self.best_res_act },
            "transducer_gio_count" : {self.transd_gio_c },
            "transducer_gq11_count" : {self.transd_gq11_c },
            "transducer_g1213_count" : {self.transd_g1213_c },
            "transducer_arrest_count" : {self.transd_arrest_c },
            "transducer_grks_count" : {self.transd_grks_c },
            "ligand_types_smallmolecule_count" : {self.ligand_types_smmol_c },
            "ligand_types_peptide_count" : {self.ligand_types_pep_c },
            "limodgrp_stimulatory_count" : {self.limodgrp_stim_c },
            "limodgrp_inhibitory_count" : {self.limodgrp_inhib_c },
            "limodgrp_unknown_count" : {self.limodgrp_unk_c },
            "limod_apo_count" : {self.limod_apo_c },
            "limod_agonist_count" : {self.limod_ag_c },
            "limod_agonist_partial_count" : {self.limod_ag_partial_c },
            "limod_allosteric_antagonist_count" : {self.limod_allo_antag_c },
            "limod_antagonist_count" : {self.limod_antag_c },
            "limod_inverse_agonist_count" : {self.limod_inv_ag_c },
            "limod_negative_allosteric_modulators_count" : {self.limod_nam_c },
            "limod_positive_allosteric_modulators_count" : {self.limod_pam_c },
            "limod_unk_count" : {self.limod_unk_c },
            '''
        return "{" + str + "}"


    class Meta:
        """Index definitions for the GpcrStructureStatisticsTable model to optimize query performance during table filtering"""
        
        db_table = "table_gpcr_structure_statistics"
        indexes = [
                models.Index(fields=["uniprot_entry_name"], name="gsst_uniprot_entry_name_idx"),
                models.Index(fields=["protein_name"], name="gsst_protein_name_idx"),
                models.Index(fields=["receptor_family"], name="gsst_receptor_family_idx"),
                models.Index(fields=["receptor_class"], name="gsst_receptor_class_idx"),
                models.Index(fields=["gene_name"], name="gsst_gene_idx"),
                models.Index(fields=["species_common_name"], name="gsst_species_idx"),
                models.Index(fields=["st_inact_c"], name="gsst_st_inact_c_idx"),
                models.Index(fields=["st_interm_c"], name="gsst_st_interm_c_idx"),
                models.Index(fields=["st_act_c"], name="gsst_st_act_c_idx"),
                models.Index(fields=["best_res_inact"], name="gsst_best_res_inact_idx"),
                models.Index(fields=["best_res_interm"], name="gsst_best_res_interm_idx"),
                models.Index(fields=["best_res_act"], name="gsst_best_res_act_idx"),
                models.Index(fields=["transd_gio_c"], name="gsst_transd_gio_c_idx"),
                models.Index(fields=["transd_gq11_c"], name="gsst_transd_gq11_c_idx"),
                models.Index(fields=["transd_g1213_c"], name="gsst_transd_g1213_c_idx"),
                models.Index(fields=["transd_arrest_c"], name="gsst_transd_arrest_c_idx"),
                models.Index(fields=["transd_grks_c"], name="gsst_transd_grks_c_idx"),
                models.Index(fields=["ligand_types_smmol_c"], name="gsst_ligand_types_smmol_c_idx"),
                models.Index(fields=["ligand_types_pep_c"], name="gsst_ligand_types_pep_c_idx"),
                models.Index(fields=["limodgrp_stim_c"], name="gsst_limodgrp_stim_c_idx"),
                models.Index(fields=["limodgrp_inhib_c"], name="gsst_limodgrp_inhib_c_idx"),
                models.Index(fields=["limodgrp_unk_c"], name="gsst_limodgrp_unk_c_idx"),
                models.Index(fields=["limod_apo_c"], name="gsst_limod_apo_c_idx"),
                models.Index(fields=["limod_ag_c"], name="gsst_limod_ag_c_idx"),
                models.Index(fields=["limod_ag_partial_c"], name="gsst_limod_ag_partial_c_idx"),
                models.Index(fields=["limod_allo_antag_c"], name="gsst_limod_allo_antag_c_idx"),
                models.Index(fields=["limod_antag_c"], name="gsst_limod_antag_c_idx"),
                models.Index(fields=["limod_inv_ag_c"], name="gsst_limod_inv_ag_c_idx"),
                models.Index(fields=["limod_nam_c"], name="gsst_limod_nam_c_idx"),
                models.Index(fields=["limod_pam_c"], name="gsst_limod_pam_c_idx"),
                models.Index(fields=["limod_unk_c"], name="gsst_limod_unk_c_idx"),
        ]


class GpcrStructureBrowserTable(models.Model):
    """Flat, per-structure table backing the structure browser (and, from round 2,
    the ligand-interaction PDB selection page). Built by
    `manage.py build_structure_browser_table` from structure.tables.structure_browser_table_query.

    Only fields that need real per-row computation (regex classification, picking
    the arrestin/G-alpha row, list-building, subqueries) are stored here — plain
    FK-chain scalars (family, class, species, state, pdb, resolution, ...) are read
    live via select_related off `structure` at request time instead, since a join
    over this table's row count (a few thousand) is effectively free.
    """

    structure = models.OneToOneField('structure.Structure', on_delete=models.CASCADE,
                                      primary_key=True, related_name='browser_row')
    arrestin_extra_protein = models.ForeignKey('structure.StructureExtraProteins', null=True,
                                                on_delete=models.SET_NULL, related_name='+')

    coverage = models.IntegerField(null=True)
    sodium_site = models.BooleanField(default=False)

    fusions = models.TextField(null=True)
    antibodies = models.TextField(null=True)

    auxiliary_molecules = models.TextField(null=True)
    auxiliary_molecule_type = models.TextField(null=True)
    auxiliary_molecule_function = models.TextField(null=True)

    ligands = JSONField(default=list)
    ligand_type = models.TextField(null=True)
    ligand_role = models.TextField(null=True)

    endo_ligands = JSONField(default=list)
    endo_type = models.TextField(null=True)

    gene_name = models.CharField(max_length=100, null=True)
    gene_entrez_url = models.TextField(null=True)
    iuphar_link = models.TextField(null=True)
    iuphar_index = models.CharField(max_length=50, null=True)

    # forward-looking, for the round-2 ligand-interaction selection page
    has_ligand_interactions = models.BooleanField(default=False)
    num_annotated_ligands = models.IntegerField(default=0)
    num_interactions = models.IntegerField(default=0)

    def __str__(self):
        return "{}".format(self.structure_id)

    class Meta:
        db_table = "table_gpcr_structure_browser"
        indexes = [
            models.Index(fields=["has_ligand_interactions"], name="gsbt_has_li_idx"),
        ]

