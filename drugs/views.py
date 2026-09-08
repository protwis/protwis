from django.shortcuts import render
from django.http import JsonResponse
from django.db.models import Count, Max, Case, When, IntegerField, Subquery, OuterRef
from django.core.cache import cache
from django.views.decorators.cache import cache_page
from django.views.generic import TemplateView
from structure.models import Structure
from common.models import WebResource, WebLink
from drugs.models import Drugs, Indication, ATCCodes, IndicationAssociation
from protein.views import get_sankey_data
from protein.models import Protein, ProteinFamily, TissueExpression, Gene
from mapper.views import DataMapperHome
from ligand.models import AssayExperiment, LigandID
from ligand.functions import standardize_smiles

import json
from collections import OrderedDict, defaultdict
from copy import deepcopy
import pandas as pd
import logging

def Venn(request, origin="both"):
    name_of_cache = 'venn_' + origin
    context = cache.get(name_of_cache)
    # NOTE cache disabled for development only!
    # context = None
    if context == None:
        context = OrderedDict()
        phases_dict = {}
        key_mapping = {
            1: "Phase I",
            2: "Phase II",
            3: "Phase III",
            4: "Phase IV"
        }

        if origin == "drugs":
            # Call to get drugs in each maximum phase
            drug_phases = Drugs.objects.all().values_list('indication_max_phase', 'ligand_id__name').distinct()
            for item in drug_phases:
                if item[0] not in phases_dict.keys():
                    phases_dict[item[0]] = []
                phases_dict[item[0]].append(item[1])
            phases_dict = {key_mapping[k]: phases_dict[k] for k in key_mapping if k in phases_dict}
            for key in list(phases_dict.keys()):
                count = len(phases_dict[key])
                new_key = f"{key} {count}"
                phases_dict[new_key] = '\n'.join(phases_dict[key])
                del phases_dict[key]
        else:
            # Call to get receptors in each maximum phase
            receptor_phases = Drugs.objects.all().values_list('indication_max_phase', 'target_id__entry_name').distinct()
            for item in receptor_phases:
                if item[0] not in phases_dict.keys():
                    phases_dict[item[0]] = []
                phases_dict[item[0]].append(item[1])
            phases_dict = {key_mapping[k]: phases_dict[k] for k in key_mapping if k in phases_dict}
            for key in list(phases_dict.keys()):
                count = len(phases_dict[key])
                new_key = f"{key} ({count})"
                phases_dict[new_key] = '\n'.join(phases_dict[key])
                del phases_dict[key]

        context["phases_dict"] = phases_dict
        context["phases_dict_keys"] = list(phases_dict.keys())

        # Fetch table data with all related information
        table_data = Drugs.objects.select_related(
            'target__family__parent__parent__parent',
            'ligand__ligand_type',
            'indication',
            'moa',
            'disease_association'
        ).annotate(
            #Fetch single gene name and entrez_id for each target using subqueries, prioritizing lowest entrez_id
            gene_name=Subquery(
                Gene.objects.filter(proteins=OuterRef('target__pk')).order_by('entrez_id').values('name')[:1]
            ),
            gene_entrez_id=Subquery(
                Gene.objects.filter(proteins=OuterRef('target__pk')).order_by('entrez_id').values('entrez_id')[:1]
            )
        ).values(
            'target',
            'target__entry_name',
            'gene_name',
            'gene_entrez_id',
            'target__name',
            'target__family__parent__name',
            'target__family__parent__parent__name',
            'target__family__parent__parent__parent__name',
            'ligand__id',
            'ligand__name',
            'ligand__ligand_type__name',
            'ligand__smiles',
            'ligand__mw',
            'ligand__sequence',
            'moa__name',
            'indication__title',
            'indication__slug',
            'indication__code',
            'indication_max_phase',
            'disease_association__association_score',
            'drug_status'
        )

        # Convert the table_data queryset to a list of dictionaries
        table_data_list = list(table_data)

        # Convert the list of dictionaries to a pandas DataFrame
        df = pd.DataFrame(table_data_list)

        # Rename the columns to your desired format
        df.rename(columns={
            'target': 'TargetID',
            'target__entry_name': "Uniprot entry",
            'gene_name': "Gene name",
            'gene_entrez_id': "Gene_entrez_id",
            'target__name': "Protein name",
            'target__family__parent__name': 'Receptor family',
            'target__family__parent__parent__name': 'Ligand type',
            'target__family__parent__parent__parent__name': 'Class',
            'ligand__id': 'LigandID',
            'ligand__name': 'Ligand name',
            'ligand__ligand_type__name': 'Drug type',
            'moa__name': 'Modality',
            'ligand__smiles': 'raw_smiles',
            'ligand__mw': 'mw',
            'ligand__sequence': 'sequence',
            'indication__slug': 'Indication Slug',
            'indication__title': 'Indication name',
            'indication__code': 'ICD11',
            'indication_max_phase': 'Phase',
            'disease_association__association_score': 'Association score',
            'drug_status': 'Approved'
        }, inplace=True)

        #Convert Gene_entrez_ids to web links
        entrez_websource = WebResource.objects.get(slug="entrez_gene")
        df['Gene_entrez_weblink'] = df['Gene_entrez_id'].apply(lambda id: str(WebLink(index=id, web_resource=entrez_websource)) if pd.notna(id) and str(id) != "" else "")
        df.drop(columns=['Gene_entrez_id'], inplace=True)

        # Preprocess SMILES data
        extra_df = df.apply(DrugSectionSelection.process_smiles, axis=1)

        # Merge the extra DataFrame with the main DataFrame
        df = pd.concat([df, extra_df], axis=1)

        # Convert 'Approved' from integer to 'Yes'/'No'
        df['Approved'] = df['Approved'].apply(lambda x: 'Yes' if x == "Approved" else 'No')

        # Drop rows with missing LigandID
        df.dropna(subset=['LigandID'], inplace=True)

        # Fill missing values for other columns
        df.fillna({
            'Ligand name': 'Unknown',
            'Phase': 'N/A',
            'Approved': 'No',
        }, inplace=True)

        # Fetch master indications
        master_indications = Indication.objects.filter(level=0).values('slug', 'title')
        master_mapping = {entry['slug']: entry['title'] for entry in master_indications}

        # Add master indication to DataFrame
        df['Master Indication'] = df['Indication Slug'].apply(
            lambda slug: master_mapping.get(slug.split("_")[0], "") if isinstance(slug, str) else "")

        # Fetch ATC data
        atc_data = ATCCodes.objects.values(
            'ligand', 'code__index', 'name_0', 'name_1'
        )
        atc_df = pd.DataFrame(list(atc_data))
        atc_df.rename(columns={
            'ligand': 'LigandID',
            'code__index': 'ATC Code',
            'name_0': 'ATC Parent Name',
            'name_1': 'ATC Name'
        }, inplace=True)

        # Merge ATC data with the main DataFrame
        df = df.merge(atc_df, on='LigandID', how='left')
        df.fillna("", inplace=True)
        df.drop_duplicates(inplace=True)

        # Add In Trial column
        df['In trial'] = df['Phase'].apply(lambda phase: 'Yes' if phase < 4 else 'No')

        # Convert DataFrame to JSON
        json_records = df.to_json(orient='records')

        # Pass the JSON data to the template context
        context['Full_data'] = json_records

    context["layout"] = origin

    return render(request, 'venn_diagrams.html', context)


##########################################
#### Drugs, Indications, and Targets #####
##########################################

class DrugSectionSelection(TemplateView):
    # variable setup #
    template_name = 'common/selection_drugs.html'

    # Initialize the page, title and description paramter
    page = 'Drugs'
    title = "Drug search"
    description = 'Search by drug name'

    def process_smiles(row):
        raw, mw = row['raw_smiles'], row['mw']
        _, smiles_for_image, picture_flag = standardize_smiles(raw, mw)
        return pd.Series({
            'smiles_for_image': smiles_for_image,
            'picture': picture_flag
        })

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)

        # Ensure that both title and description are initialized
        title = self.title
        description = self.description  # Set default description
        page = self.page

         # Modify the description and fetch data based on the page type
        if page == 'Drugs':
            description = 'Search by agent / drug name'

            # Fetch distinct ligands and their names
            table_data = Drugs.objects.select_related('ligand').values(
                'ligand',  # Agent/Drug id
                'ligand__name',  # Agent/Drug name
            ).distinct()

            # Convert `table_data` to a DataFrame
            table_df = pd.DataFrame(table_data)

            # Ensure we have a unique list of ligand IDs
            unique_ligand_ids = table_df['ligand'].unique().tolist()

            # Fetch related resources: `web_resource__name` and `index` for the unique ligand IDs
            resource_data = LigandID.objects.filter(
                ligand__in=unique_ligand_ids
            ).values(
                'ligand',  # Ligand ID
                'web_resource__name',  # Web resource name (e.g., "PubChem", "ChEMBL")
                'index',  # Identifier ID
            )

            # Convert `resource_data` queryset to a list of dictionaries
            resource_data_list = list(resource_data)

            def clean_index(val):
                val_str = str(val)
                if val_str.endswith('.0'):
                    return val_str[:-2]  # Removes the trailing ".0"
                return val_str

            # Transform `resource_data_list` into a structure for merging
            resource_dict = defaultdict(dict)
            for row in resource_data_list:
                ligand_id = row['ligand']
                resource_name = str(row['web_resource__name'])
                resource_index = clean_index(row['index'])
                resource_dict[ligand_id][resource_name] = resource_index

            # Add resource data to the table DataFrame by mapping from `resource_dict`
            table_df['Resources'] = table_df['ligand'].map(resource_dict)

            # Flatten the `Resources` dictionary into separate columns for each resource
            resource_df = table_df['Resources'].apply(pd.Series)

            # Merge `table_df` and `resource_df`
            final_df = pd.concat([table_df, resource_df], axis=1).drop(columns=['Resources'])

            # Rename columns for better readability
            final_df.rename(columns={
                'ligand': 'LigandID',
                'ligand__name': 'Ligand Name',
                'ChEMBL_compound_ids': 'ChEMBL',  # Rename for consistency
                'Guide To Pharmacology': 'GTP',  # Rename for consistency
                'Drug Central': 'DrugCentral',   # Rename for consistency
            }, inplace=True)

            # Convert the final DataFrame to JSON
            json_records = final_df.to_json(orient='records')

            # Pass the JSON data to the template context
            context['Drugs_identifier_table'] = json_records

        elif page == 'Targets':
            description = 'Search by target name'
            # Fetch distinct targets and create a dictionary of {target.name: target.id}
            search_data = Drugs.objects.all().prefetch_related('target').distinct('target')
            search_dict = {
                f"{drug.target.name.replace(' receptor', '')} ({drug.target.entry_name.replace('_human', '').upper()})": drug.target.id
                for drug in search_data
            }

            # Fetch ATC codes for the table
            ATC_data = ATCCodes.objects.select_related(
                'code'
            ).values(
                'ligand', # Agent/Drug id
                'code__index' # ATC codes
            )

            # Convert ATC_data queryset to a list of dictionaries
            ATC_data_list = list(ATC_data)

            # Create a DataFrame from the ATC data
            atc_df = pd.DataFrame(ATC_data_list)

            # Group ATC codes by 'Ligand ID' and concatenate them
            atc_df_grouped = atc_df.groupby('ligand')['code__index'].agg(', '.join).reset_index()

            # Rename columns for the ATC DataFrame
            atc_df_grouped.rename(columns={'ligand': 'LigandID', 'code__index': 'ATC'}, inplace=True)

            # Fetch table data with all related information
            table_data = Drugs.objects.select_related(
                'target',
                'ligand__ligand_type',
                'indication',
                'moa',
                'disease_association'
            ).values(
                'target',  # Target ID
                'target__entry_name', # Gene Name
                'target__name',  # Target name
                'ligand__name',  # Agent/Drug
                'ligand__ligand_type__name',  # "type'
                'ligand__smiles', # SMILES
                'ligand__mw',
                'ligand__sequence',
                'moa__name',  # Modality
                'indication__title',  # Disease name
                'indication__code',  # Disease ICD11 code
                'indication_max_phase',  # Max phase
                'ligand', # ligand ID
                'disease_association__association_score', # Disease association score
                'drug_status',  # Approval
            )

            # Convert the table_data queryset to a list of dictionaries
            table_data_list = list(table_data)

            # Convert the list of dictionaries to a pandas DataFrame
            df = pd.DataFrame(table_data_list)

            # Rename the columns to your desired format
            df.rename(columns={
                'target': 'Target ID',
                'target__entry_name': 'Gene name',
                'target__name': 'Target name',
                'ligand__name': 'Ligand name',
                'ligand__ligand_type__name': 'Type',
                'ligand__smiles': 'raw_smiles',
                'ligand__mw': 'mw',
                'ligand__sequence': 'sequence',
                'moa__name': 'Modality',
                'indication__title': 'Indication name',
                'indication__code': 'ICD11',
                'indication_max_phase': 'Phase',
                'ligand': 'LigandID',
                'disease_association__association_score' : 'Association score',
                'drug_status': 'Status'
            }, inplace=True)

            # Preprocess SMILES data
            extra_df = df.apply(DrugSectionSelection.process_smiles, axis=1)

            # Merge the extra DataFrame with the main DataFrame
            df = pd.concat([df, extra_df], axis=1)

            # Merge the ATC data into the main DataFrame (df) on 'Ligand ID'
            df = df.merge(atc_df_grouped, on='LigandID', how='left')
            # Fill NaN values in the 'ATC' column with None
            group_cols = ['Target ID', 'Target name','Gene name', 'LigandID', 'Ligand name', 'Indication name', 'Type', 'Modality', 'ICD11', 'ATC', 'Association score']
            df[group_cols] = df[group_cols].fillna({
                'Target name': 'Unknown',
                'Gene name': 'Unknown',
                'Ligand name': 'Unknown',
                'Indication name': 'Unknown',
                'Type': 'Unknown',
                'Modality': 'Unknown',
                'ICD11': 'Unspecified',
                'ATC': '',
                'Association score': 0
            })

            # Precompute flags for phase and approval status
            df['Is_Approved'] = (df['Status'] == 'Approved').astype(int)

            # Group by necessary columns and perform aggregation in one go
            grouped = df.groupby(
                ['Target ID', 'Target name','Gene name', 'LigandID', 'Ligand name', 'Indication name', 'Type', 'Modality', 'ICD11', 'ATC', 'Association score']
            )

            # Perform aggregation
            Modified_df = grouped.agg(
                Highest_phase=('Phase', 'max'),  # Get the highest phase for each group
                Approved=('Is_Approved', 'max'),  # Check if any row has 'Approved' status (max of binary flag)
                sequence=('sequence', 'first'),
                smiles_for_image=('smiles_for_image', 'first'),
                picture=('picture', 'first')
            ).reset_index()

            # Convert 'Approved' from integer to 'Yes'/'No'
            Modified_df['Approved'] = Modified_df['Approved'].apply(lambda x: 'Yes' if x == 1 else 'No')

            # Add In Trial column
            Modified_df['In trial'] = df['Phase'].apply(lambda phase: 'Yes' if phase < 4 else 'No')

            # Convert DataFrame to JSON
            json_records = Modified_df.to_json(orient='records')

            # Pass the JSON data to the template context
            context['Full_data'] = json_records
        elif page == 'Indications':
            description = 'Search by indication name'
            # Fetch distinct indications and create a dictionary of {indication.name: indication.id}
            search_data = Drugs.objects.all().prefetch_related('indication').distinct('indication')
            search_dict = {drug.indication.title: drug.indication.id for drug in search_data}

            # ###########################
            # Single Data Query
            # ###########################
            # Fetch all data in a single query
            table_data = Drugs.objects.select_related(
                'indication__code',
                'ligand__ligand_type',
                'target__family__parent__parent__parent',  # All target info
                'moa',
                'disease_association'
            ).annotate(
            #Fetch single gene name and entrez_id for each target using subqueries, prioritizing lowest entrez_id
            gene_name=Subquery(
                Gene.objects.filter(proteins=OuterRef('target_id')).order_by('entrez_id').values('name')[:1]
                ),
            gene_entrez_id=Subquery(
                Gene.objects.filter(proteins=OuterRef('target_id')).order_by('entrez_id').values('entrez_id')[:1]
                )
            ).values(
                'indication',  # Indication ID
                'indication__title',  # Indication name
                'indication__code',  # Disease ICD11 code
                'ligand__id',  # Ligand ID
                'ligand__name',  # Ligand name
                'ligand__smiles', # SMILES
                'ligand__mw',
                'ligand__sequence',
                'indication_max_phase',  # Max phase
                'drug_status',  # Approval
                'ligand__ligand_type__name',  # Molecule type
                'moa__name',  # Mode of action
                'target__entry_name',  # Uniprot name
                'gene_name', # Gene name
                'gene_entrez_id', # Gene id
                'target__name',  # Protein name
                'target__family__parent__name',  # Receptor family
                'target__family__parent__parent__name',  # Ligand type
                'target__family__parent__parent__parent__name',  # Class
                'disease_association__association_score', # Disease association score
                "disease_association__ot_genetics_portal", # OT Genetics
                "disease_association__gene_burden", # Gene Burden
                "disease_association__eva", # ClinVar
                "disease_association__genomics_england", # GEL PanelApp
                "disease_association__gene2phenotype", # Gene2phenotype
                "disease_association__uniprot_literature", # UniProt literature
                "disease_association__uniprot_variants", # UniProt curated variants
                "disease_association__orphanet", # Orphanet
                "disease_association__clingen", # ClinGen
                "disease_association__cancer_gene_census", # Cancer Gene Census
                "disease_association__intogen", # IntOGen
                "disease_association__eva_somatic", # Clinvar (somatic)
                "disease_association__cancer_biomarkers", # Cancer Biomarkers
                "disease_association__chembl", # ChEMBL
                "disease_association__crispr_screen", # CRISPR Screens
                "disease_association__crispr", # Project Score
                "disease_association__slapenrich", # SLAPenrich
                "disease_association__reactome", # Reactome
                "disease_association__sysbio", # Gene signatures
                "disease_association__europepmc", # Europe PMC
                "disease_association__expression_atlas", # Expression Atlas
                "disease_association__impc" # IMPC
            )

            # Convert the table_data queryset to a list of dictionaries
            table_data_list = list(table_data)

            # Convert the list of dictionaries to a pandas DataFrame
            df = pd.DataFrame(table_data_list)

            # Rename the columns to your desired format
            df.rename(columns={
                'indication': 'Indication ID',
                'indication__title': 'Indication name',
                'indication__code': 'ICD11',
                'ligand__id': 'LigandID',
                'ligand__name': 'Drug name',
                'ligand__smiles': 'raw_smiles',
                'ligand__mw': 'mw',
                'ligand__sequence': 'sequence',
                'indication_max_phase': 'Phase',
                'drug_status': 'Status',
                'ligand__ligand_type__name': 'Molecule_type',
                'moa__name': 'Mode of action',
                'target__entry_name': 'Uniprot entry',
                'gene_name' : 'Gene name',
                'gene_entrez_id' : 'Gene_entrez_id',
                'target__name': 'Protein name',
                'target__family__parent__name': 'Receptor family',
                'target__family__parent__parent__name': 'Ligand type',
                'target__family__parent__parent__parent__name': 'Class',
                'disease_association__association_score' : 'Association score',
                "disease_association__ot_genetics_portal": "OT Genetics",
                "disease_association__gene_burden": "Gene Burden",
                "disease_association__eva": "ClinVar",
                "disease_association__genomics_england": "GEL PanelApp",
                "disease_association__gene2phenotype": "Gene2phenotype",
                "disease_association__uniprot_literature": "UniProt literature",
                "disease_association__uniprot_variants": "UniProt curated variants",
                "disease_association__orphanet": "Orphanet",
                "disease_association__clingen": "ClinGen",
                "disease_association__cancer_gene_census": "Cancer Gene Census",
                "disease_association__intogen": "IntOGen",
                "disease_association__eva_somatic": "Clinvar (somatic)",
                "disease_association__cancer_biomarkers": "Cancer Biomarkers",
                "disease_association__chembl":"ChEMBL",
                "disease_association__crispr_screen": "CRISPR Screens",
                "disease_association__crispr": "Project Score",
                "disease_association__slapenrich": "SLAPenrich",
                "disease_association__reactome": "Reactome",
                "disease_association__sysbio": "Gene signatures",
                "disease_association__europepmc": "Europe PMC",
                "disease_association__expression_atlas": "Expression Atlas",
                "disease_association__impc": "IMPC"
            }, inplace=True)

            # 1) Apply standardization
            extra_df = df.apply(DrugSectionSelection.process_smiles, axis=1)

            # 2) Merge them back
            df = pd.concat([df, extra_df], axis=1)

            # Define MOA categories
            stim_moa = ['Partial agonist', 'Agonist', 'PAM']
            inhib_moa = ['Antagonist', 'Inverse agonist', 'NAM']

            #Convert Gene_entrez_ids to web links
            entrez_websource = WebResource.objects.get(slug="entrez_gene")
            df['Gene_entrez_weblink'] = df['Gene_entrez_id'].apply(lambda id: str(WebLink(index=id, web_resource=entrez_websource)) if pd.notna(id) and str(id) != "" else "")
            df.drop(columns=['Gene_entrez_id'], inplace=True)

            # Split the DataFrame into two: one for targets and one for drugs
            df_targets = df.copy()
            df_drugs = df.copy()

            # ###########################
            # Data Aggregation for Targets
            # ###########################
            # Update group_cols to include 'Indication name'
            group_cols = ['Indication ID', 'ICD11', 'Uniprot entry', 'Gene name', 'Gene_entrez_weblink', 'Indication name', 'Protein name', 'Receptor family', 'Ligand type', 'Class']

            # Define disease association columns to be added to the grouping
            disease_cols = [
                'Association score', 'OT Genetics', 'Gene Burden', 'ClinVar', 'GEL PanelApp',
                'Gene2phenotype', 'UniProt literature', 'UniProt curated variants', 'Orphanet',
                'ClinGen', 'Cancer Gene Census', 'IntOGen', 'Clinvar (somatic)', 'Cancer Biomarkers',
                'ChEMBL', 'CRISPR Screens', 'Project Score', 'SLAPenrich', 'Reactome', 'Gene signatures',
                'Europe PMC', 'Expression Atlas', 'IMPC'
            ]

            # Add these to group_cols
            group_cols += disease_cols

            df_targets[disease_cols] = df[disease_cols].fillna("")

            # Precompute classifications
            df_targets['Classification'] = df_targets['Status'].apply(lambda x: 'Drug' if x == 'Approved' else 'Agent')

            # Pre-filter the DataFrame for stimulatory and inhibitory MOAs
            stim_moa_df = df_targets[df_targets['Mode of action'].isin(stim_moa)]
            inhib_moa_df = df_targets[df_targets['Mode of action'].isin(inhib_moa)]


            # Create a function to get the max phase safely
            def get_max_phase(df, group_cols):
                if df.empty:
                    return pd.Series(dtype='float64')  # Return an empty Series if DataFrame is empty
                return df.groupby(group_cols)['Phase'].max()

            # Group the targets DataFrame by the updated columns
            grouped_targets = df_targets.groupby(group_cols)

            # Aggregate data using efficient vectorized operations
            agg_data_targets = grouped_targets.agg(
                All_max_phase=('Phase', 'max'),
                All_Drugs=('Classification', lambda x: (x == 'Drug').sum()),
                All_Agents=('Classification', lambda x: (x == 'Agent').sum()),
                sequence=('sequence', 'first')
            ).reset_index()

            # Compute the stimulatory and inhibitory max phase and counts efficiently
            stim_max_phase = get_max_phase(stim_moa_df, group_cols)
            inhib_max_phase = get_max_phase(inhib_moa_df, group_cols)

            # Add the precomputed data into the DataFrame
            agg_data_targets = agg_data_targets.merge(stim_max_phase.rename('Stimulatory_max_phase'), on=group_cols, how='left')
            agg_data_targets = agg_data_targets.merge(inhib_max_phase.rename('Inhibitory_max_phase'), on=group_cols, how='left')

           # Compute Stimulatory and Inhibitory Drugs and Agents counts
            # Reindex to ensure both 'Drug' and 'Agent' columns are present, even if they don't exist in the data
            stim_counts = stim_moa_df.groupby(group_cols)['Classification'].value_counts().unstack(fill_value=0).reindex(
                columns=['Drug', 'Agent'], fill_value=0).reset_index()

            # Rename columns to make sure they have meaningful names
            stim_counts.rename(columns={'Drug': 'Stimulatory_Drugs', 'Agent': 'Stimulatory_Agents'}, inplace=True)

            # Reindex to ensure both 'Drug' and 'Agent' columns are present for inhibitory counts
            inhib_counts = inhib_moa_df.groupby(group_cols)['Classification'].value_counts().unstack(fill_value=0).reindex(
                columns=['Drug', 'Agent'], fill_value=0).reset_index()

            # Rename columns for inhibitory MOAs
            inhib_counts.rename(columns={'Drug': 'Inhibitory_Drugs', 'Agent': 'Inhibitory_Agents'}, inplace=True)

            # Merge the counts back into the aggregated data
            agg_data_targets = agg_data_targets.merge(stim_counts[group_cols + ['Stimulatory_Drugs', 'Stimulatory_Agents']], on=group_cols, how='left').fillna(0)
            agg_data_targets = agg_data_targets.merge(inhib_counts[group_cols + ['Inhibitory_Drugs', 'Inhibitory_Agents']], on=group_cols, how='left').fillna(0)

            # Convert DataFrame to JSON
            json_records_targets = agg_data_targets.to_json(orient='records')
            context['Full_data_targets'] = json_records_targets

            # ###########################
            # Data Aggregation for Drugs
            # ###########################
            group_cols_drugs = ['Indication ID', 'ICD11', 'Uniprot entry', 'Gene name', 'Gene_entrez_weblink', 'Drug name', 'LigandID', 'Indication name',
                                'Protein name', 'Receptor family', 'Ligand type', 'Class', 'Molecule_type','Mode of action','Phase']

            # Precompute relevant columns
            df_drugs['Is_Approved'] = df_drugs['Status'].apply(lambda x: 1 if x == 'Approved' else 0)
            df_drugs['Is_Phase_I'] = (df_drugs['Phase'] == 1).astype(int)
            df_drugs['Is_Phase_II'] = (df_drugs['Phase'] == 2).astype(int)
            df_drugs['Is_Phase_III'] = (df_drugs['Phase'] == 3).astype(int)

            # Group the DataFrame only once
            grouped_drugs = df_drugs.groupby(group_cols_drugs)

            # Aggregate data using vectorized operations
            agg_data_drugs = grouped_drugs.agg(
                Highest_phase=('Phase', 'max'),
                Phase_I_trials=('Is_Phase_I', 'sum'),
                Phase_II_trials=('Is_Phase_II', 'sum'),
                Phase_III_trials=('Is_Phase_III', 'sum'),
                Approved=('Is_Approved', 'max'),
                sequence=('sequence', 'first'),
                smiles_for_image=('smiles_for_image', 'first'),
                picture=('picture', 'first')
            ).reset_index()

            # Convert 'Approved' from integer to 'Yes'/'No'
            agg_data_drugs['Approved'] = agg_data_drugs['Approved'].apply(lambda x: 'Yes' if x == 1 else 'No')

            # Add In Trial column
            agg_data_drugs['In trial'] = agg_data_drugs['Phase'].apply(lambda phase: 'Yes' if phase < 4 else 'No')

            # Convert DataFrame to JSON
            json_records_drugs = agg_data_drugs.to_json(orient='records')
            context['Full_data_drugs'] = json_records_drugs

        # Convert to JSON string and pass to context
        if page == 'Drugs':
            pass
        else:
            context['search_data'] = json.dumps(search_dict)
        context['page'] = page
        context['title'] = title
        context['description'] = description
        # context['table_data'] = table_data

        return context

    def get_entry_name_by_target_id(self, target_id):
        """
        Retrieves the protein.entry_name based on the given drug.target.id.

        Args:
            target_id (int): The ID of the target (Protein).

        Returns:
            str or None: The entry_name of the associated Protein, or None if not found.
        """
        try:
            protein = Protein.objects.get(id=target_id)
            return protein.entry_name
        except Protein.DoesNotExist:
            return None

# Set up logging
logger = logging.getLogger(__name__)

def fetch_sankey_data_view(request):
    logger.info("fetch_sankey_data_view called")  # Log the function call
    entry_name = request.GET.get('entry_name', None)  # Extract `entry_name` from the request

    if not entry_name:
        logger.error("No entry_name provided in the request")
        return JsonResponse({'error': 'Missing entry_name'}, status=400)

    logger.info(f"Received entry_name: {entry_name}")  # Log the entry name

    try:
        sankey_data = get_sankey_data(entry_name)  # Call your helper function
        if sankey_data:
            logger.info(f"Sankey data successfully retrieved for entry_name: {entry_name}")
            return JsonResponse({'sankey_data': sankey_data}, status=200)
        else:
            logger.warning(f"No data found for entry_name: {entry_name}")
            return JsonResponse({'error': 'No data found for this entry'}, status=404)
    except Exception:
        logger.exception(f"An error occurred while fetching sankey data for entry_name: {entry_name}")
        return JsonResponse({'error': 'An internal server error occurred'}, status=500)

class DruggedGPCRome(TemplateView):
    """
    Per class statistics of known ligands.
    """

    template_name = 'drugged_gpcrome.html'

    def get_context_data(self, **kwargs):

        context = super().get_context_data(**kwargs)

        drug_count_receptor_dict = {}

        # Use iterator to process the queryset in chunks
        drug_data = Drugs.objects.all().values_list('target__entry_name', 'indication_max_phase').distinct()

        all_proteins = Protein.objects.filter(species_id=1, parent_id__isnull=True, accession__isnull=False, family_id__slug__startswith='0').exclude(
                                            family_id__slug__startswith='007'
                                        ).exclude(
                                            family_id__slug__startswith='008'
                                        )
        for prot in all_proteins:
            drug_count_receptor_dict[prot.entry_name] = 0

        for pair in drug_data:
            if pair[0] in drug_count_receptor_dict.keys():
                if int(pair[1]) > int(drug_count_receptor_dict[pair[0]]):
                    drug_count_receptor_dict[pair[0]] = int(pair[1])

        # Mapping for Value2
        status_color_map = {
            1: '#f5bcbf',
            2: '#f17270',
            3: '#dd2628',
            4: '#2c87c8'
        }

        # Creating the new dictionary
        results_updated_dict = {
            key: {
                'Value1': value,
                'Value2': status_color_map.get(value, '#FFFFFF')  # Default to black if key not found
            }
            for key, value in drug_count_receptor_dict.items()
                }

        gpcr_data = DataMapperHome.GenerateGPCRomeDataStructure(data_type="Classic")
        updated_data = DataMapperHome.update_nested_GPCRome_data(gpcr_data["Data"], results_updated_dict)

        context['GPCRomeData'] = json.dumps(updated_data)

        #TREE SECTION
        drug_data = Drugs.objects.all().values_list('target__entry_name', 'ligand__name','indication_max_phase')

        drug_dict = {}

        # Populate the dictionary
        for drug in drug_data:
            target, ligand, phase = drug
            if target not in drug_dict:
                # Initialize with empty lists
                drug_dict[target] = {'Outer1': [], 'Outer2': [], 'Outer3': [], 'Outer4': [], 'Inner': []}

            if phase in [1, 2, 3, 4]:
                outer_key = f"Outer{phase}"
                drug_dict[target][outer_key].append(ligand)  # Add ligand to the corresponding Outer list
            else:
                drug_dict[target]['Inner'].append(ligand)  # Add ligand to the Inner list

        # Convert lists to unique counts
        for target in drug_dict:
            for key in drug_dict[target]:
                unique_entries = set(drug_dict[target][key])  # Get unique ligands
                # if target == 'drd2_human' and key == 'Outer4':
                #     print(len(unique_entries),sorted(unique_entries))
                drug_dict[target][key] = len(unique_entries)  # Replace the list with the count

        tree, tree_options, circles, receptors, genes = DataMapperHome.generate_tree_plot(drug_dict)
        #Remove 0 circles
        for key, outer_dict in circles.items():
            circles[key] = {k: v for k, v in outer_dict.items() if v != 0}

        context['tree'] = json.dumps(tree)
        context['tree_options'] = tree_options
        context['circles'] = json.dumps(circles)
        context['Tree_Receptor_dict'] = json.dumps(receptors)
        context['Tree_Entrez_dict'] = json.dumps(genes)

        #REPURPOSED TREE SECTION
        # Red: IndicationMaxPhase < 4 and MaxPhase (compound) < 4 [New agents in trial (lack approval)]
        # Purple: IndicationMaxPhase < 4 and MaxPhase = 4         [Drugs being repurposed in trials (have approval)]
        # Blue: IndicationMaxPhase and MaxPhase = 4               [All drugs (both those being repurposed in trials and not)]

        drug_data = Drugs.objects.all().values_list('ligand__name', 'target__entry_name', 'indication_max_phase').distinct()
        phase4 = {}
        for drug in drug_data:
            if drug[0] not in phase4.keys():
                phase4[drug[0]] = []
            if drug[2] == 4:
                phase4[drug[0]].append(drug[1])
        phase4 = {k: v for k, v in phase4.items() if v != []}
        # Red: IndicationMaxPhase < 4 and MaxPhase (compound) < 4 [New agents in trial (lack approval)]
        # Purple: IndicationMaxPhase < 4 and MaxPhase = 4         [Drugs being repurposed in trials (have approval)]
        # Blue: IndicationMaxPhase and MaxPhase = 4               [All drugs (both those being repurposed in trials and not)]
        drug_dict = {}
        for drug in drug_data:
            if drug[1] not in drug_dict.keys():
                drug_dict[drug[1]] = {'Outer1': 0, 'Outer2': 0, 'Outer3': 0, 'Outer4': 0, 'Inner': 0}
            #Approved Drug
            if drug[2] == 4:
                drug_dict[drug[1]]['Outer3'] += 1
            else:
                #Repurposed Drug
                if drug[0] in phase4.keys():
                    drug_dict[drug[1]]['Outer2'] += 1
                # New Agent
                else:
                    drug_dict[drug[1]]['Outer1'] += 1

        repurposed_tree, repurposed_tree_options, repurposed_circles, repurposed_receptors, repurposed_genes = DataMapperHome.generate_tree_plot(drug_dict)
        #Remove 0 circles
        for key, outer_dict in repurposed_circles.items():
            repurposed_circles[key] = {k: v for k, v in outer_dict.items() if v != 0}

        context['rep_tree'] = json.dumps(repurposed_tree)
        context['rep_tree_options'] = repurposed_tree_options
        context['rep_circles'] = json.dumps(repurposed_circles)
        context['rep_Receptor_dict'] = json.dumps(repurposed_receptors)
        context['rep_Entrez_dict'] = json.dumps(repurposed_genes)

        return context

class DiseaseOverview(TemplateView):
    """
    Per class statistics of known ligands.
    """

    template_name = 'disease_overview.html'

    def get_context_data(self, **kwargs):

        context = super().get_context_data(**kwargs)

        def find_max_values(node, node_name):
            """
            Given a node (a dictionary) and its name (a string), return:
            (max_red, red_node, max_purple, purple_node, max_blue, blue_node)
            for the entire subtree starting at `node`.

            If a color key is missing at a node, treat its value as 0.
            """
            # Get current node's color values (default to 0 if missing)
            local_red = node.get('red', 0)
            local_purple = node.get('purple', 0)
            local_blue = node.get('blue', 0)

            # Initialize maxima with the current node's values
            max_red = local_red
            red_node = node_name
            max_purple = local_purple
            purple_node = node_name
            max_blue = local_blue
            blue_node = node_name

            # Iterate over children. Children are identified as dictionary values.
            for key, value in node.items():
                if isinstance(value, dict):
                    # Recurse into child
                    r_val, r_k, p_val, p_k, b_val, b_k = find_max_values(value, key)

                    # Compare and update max_red
                    if r_val > max_red:
                        max_red = r_val
                        red_node = r_k

                    # Compare and update max_purple
                    if p_val > max_purple:
                        max_purple = p_val
                        purple_node = p_k

                    # Compare and update max_blue
                    if b_val > max_blue:
                        max_blue = b_val
                        blue_node = b_k

            return max_red, red_node, max_purple, purple_node, max_blue, blue_node

        def find_node_by_name(nested_dict, target_name):
            """
            Recursively search nested_dict for a node keyed by target_name.
            Returns the dictionary for that node if found, else None.
            """
            for key, value in nested_dict.items():
                if key == target_name:
                    # Found the node
                    return value
                if isinstance(value, dict):
                    # Recurse into children
                    result = find_node_by_name(value, target_name)
                    if result is not None:
                        return result
            return None

        def increment_color_value(nested_dict, target_leaf, color):
            """
            Search through nested_dict for target_leaf. Once found, increment the
            value associated with 'color' by 1.
            Returns True if successful, False if not found.
            """
            for key, value in nested_dict.items():
                if key == target_leaf:
                    # Found the leaf node
                    if color in value:
                        # Handle case if value is None or not initialized
                        if value[color] is None:
                            value[color] = 0
                        value[color] += 1
                    else:
                        # If the color key doesn't exist, initialize it at 1
                        value[color] = 1
                    return True

                # If the value is a nested dict, recurse into it
                if isinstance(value, dict):
                    if increment_color_value(value, target_leaf, color):
                        return True

            # If not found in this branch, return False
            return False

        def gradient_color(value_fraction, start_color=(255,255,255), end_color=(0,0,0)):

            r = int(start_color[0] + (end_color[0] - start_color[0]) * value_fraction)
            g = int(start_color[1] + (end_color[1] - start_color[1]) * value_fraction)
            b = int(start_color[2] + (end_color[2] - start_color[2]) * value_fraction)

            return f"#{r:02x}{g:02x}{b:02x}"

        def find_global_maxima(node, global_max_values):
            """
            Recursively traverse 'node' and update global_max_values.
            Any node can have 'red', 'purple', 'blue' keys.
            Missing keys default to 0.
            """
            red_val = node.get('red', 0)
            purple_val = node.get('purple', 0)
            blue_val = node.get('blue', 0)

            if red_val > global_max_values['red']:
                global_max_values['red'] = red_val
            if purple_val > global_max_values['purple']:
                global_max_values['purple'] = purple_val
            if blue_val > global_max_values['blue']:
                global_max_values['blue'] = blue_val

            # Recurse into child dictionaries
            for v in node.values():
                if isinstance(v, dict):
                    find_global_maxima(v, global_max_values)

        def apply_colors(node, global_max_values, color_targets):
            """
            Recursively apply gradient colors to each node based on its values and global maxima.
            Replace 'red', 'purple', 'blue' values with their corresponding gradient colors.
            """
            # Extract node's own values (default 0 if missing)
            red_val = node.get('red', 0)
            purple_val = node.get('purple', 0)
            blue_val = node.get('blue', 0)

            red_frac = red_val / global_max_values['red'] if global_max_values['red'] > 0 else 0
            purple_frac = purple_val / global_max_values['purple'] if global_max_values['purple'] > 0 else 0
            blue_frac = blue_val / global_max_values['blue'] if global_max_values['blue'] > 0 else 0

            # Compute gradient colors
            # The function gradient_color(frac, start_color, end_color) must be defined beforehand.
            new_node = {}
            new_node['red'] = gradient_color(red_frac, (255,255,255), color_targets['red'])
            new_node['purple'] = gradient_color(purple_frac, (255,255,255), color_targets['purple'])
            new_node['blue'] = gradient_color(blue_frac, (255,255,255), color_targets['blue'])

            # Recurse into children
            for k, v in node.items():
                if isinstance(v, dict):
                    new_node[k] = apply_colors(v, global_max_values, color_targets)
                else:
                    # For non-dict, non-color keys, if you need to preserve them, reassign here
                    # Check that you're not overwriting the colors keys we just set
                    if k not in ['red', 'purple', 'blue']:
                        new_node[k] = v

            return new_node

        def convert_dict_to_colors(data_dict):
            # Define the target colors for each metric
            color_targets = {
                "red": (255,0,0),
                "purple": (128,0,128),
                "blue": (0,0,255)
            }

            # 1. Find global maxima for each color across the entire nested structure
            global_max_values = {"red": 0, "purple": 0, "blue": 0}
            find_global_maxima(data_dict, global_max_values)

            # 2. Apply gradient colors throughout the nested dict
            result = apply_colors(data_dict, global_max_values, color_targets)
            return result

        def color_style(color):
            # If the color is white, add a border style. Otherwise, no border.
            if color == '#ffffff':
                return f"background-color: {color}; border: 1px solid #ccc;"
            else:
                return f"background-color: {color};"

        def propagate_max_colors(node):
            """
            Given a nested dictionary node:
            - Keys that are colors (red, purple, blue) have integer values.
            - Other keys may lead to child dictionaries.
            This function updates each node's (red, purple, blue) to be the maximum
            found in that node or any of its descendants.

            Returns a tuple (max_red, max_purple, max_blue) for the subtree.
            """

            # Initialize current node's colors. If missing, default to 0.
            current_red = node.get('red', 0)
            current_purple = node.get('purple', 0)
            current_blue = node.get('blue', 0)

            # Check children
            for value in node.values():
                if isinstance(value, dict):
                    # Recursively propagate max colors down the subtree
                    r_val, p_val, b_val = propagate_max_colors(value)

                    # Update current node's max based on child
                    if r_val > current_red:
                        current_red = r_val
                    if p_val > current_purple:
                        current_purple = p_val
                    if b_val > current_blue:
                        current_blue = b_val

            # Update the node's own values after checking children
            node['red'] = current_red
            node['purple'] = current_purple
            node['blue'] = current_blue

            return current_red, current_purple, current_blue

        def render_nested_structure(data):
            html = []
            html.append('<ul class="no-bullet-indent">')
            for key, value in data.items():
                if key in ('red', 'purple', 'blue'):
                    continue

                # Safely extract colors
                red = value.get('red', '#ffffff') if isinstance(value, dict) else '#ffffff'
                purple = value.get('purple', '#ffffff') if isinstance(value, dict) else '#ffffff'
                blue = value.get('blue', '#ffffff') if isinstance(value, dict) else '#ffffff'

                # Check if node has children
                has_children = False
                if isinstance(value, dict):
                    for v in value.values():
                        if isinstance(v, dict):
                            has_children = True
                            break

                if isinstance(value, dict) and has_children:
                    # Node with children
                    html.append(f"""
                    <li>
                        <details>
                            <summary>
                                <span class="dot" style="{color_style(red)}"></span>
                                <span class="dot" style="{color_style(purple)}"></span>
                                <span class="dot" style="{color_style(blue)}"></span>
                                <b class="hover-pointer">{key}</b>
                            </summary>
                            {render_nested_structure(value)}
                        </details>
                    </li>
                    """)
                else:
                    # Leaf node
                    if isinstance(value, dict):
                        # Leaf node with colors
                        html.append(f"""
                        <li>
                            <span class="dot" style="{color_style(red)}"></span>
                            <span class="dot" style="{color_style(purple)}"></span>
                            <span class="dot" style="{color_style(blue)}"></span>
                            {key}
                        </li>
                        """)
                    else:
                        # Leaf node without colors (just a string or value)
                        # No colors defined, but if you want, you can handle differently
                        html.append(f"<li>{key}: {value}</li>")
            html.append('</ul>')
            return ''.join(html)

        lvl_0_1 = Indication.objects.filter(level__in=[0,1]).order_by('slug')
        all_indis = Indication.objects.all().order_by('slug')
        indication_data = {}
        listdata = {}
        listplot_data_variables = {}
        Label_Conversion = {}
        data_types_list = {'Col1': "Continuouos", 'Col2': "Continuouos", 'Col3': "Continuouos", 'Col4': "Continuouos",}
        #Gather data and setup master dict

        for indi in lvl_0_1:
           if indi.level == 0 and indi.title not in indication_data.keys():
               indication_data[indi.title] = {}
               listdata[indi.title] = []
           if indi.level == 1:
               lvl0 = indi.get_level_0().title
               if lvl0 not in indication_data.keys():
                   indication_data[lvl0] = {}
                   listdata[lvl0] = []
               if indi.title not in indication_data[lvl0].keys():
                   indication_data[lvl0][indi.title] = {'red': 0, 'purple': 0, 'blue': 0}
                   listdata[lvl0].append(indi.title)
                   listplot_data_variables[indi.title] = {'Value1': 'Circle', 'Value2': 0, 'Value3': 'Circle', 'Value4': 0, 'Value5': 'Circle','Value6': 0}
                   Label_Conversion[indi.title] = indi.title
        #with master dict set, now we add data from drugs
        drug_data = Drugs.objects.all().prefetch_related('ligand', 'indication', 'indication__parent',
                                                         'indication__parent__parent', 'indication__parent__parent__parent',
                                                         'indication__parent__parent__parent__parent',
                                                         'indication__parent__parent__parent__parent__parent')

        drug_dict = {}
        for item in drug_data:
            if item.ligand.name not in drug_dict.keys():
                drug_dict[item.ligand.name] = {'Max Phase': 0}
            if item.indication.get_level_1().title not in drug_dict[item.ligand.name].keys():
                drug_dict[item.ligand.name][item.indication.get_level_1().title] = [item.indication_max_phase, item.indication.get_level_0().title]
            if item.indication_max_phase > drug_dict[item.ligand.name]['Max Phase']:
                drug_dict[item.ligand.name]['Max Phase'] = item.indication_max_phase

        # Red: IndicationMaxPhase < 4 and MaxPhase (compound) < 4 [New agents in trial (lack approval)]
        # Purple: IndicationMaxPhase < 4 and MaxPhase = 4         [Drugs being repurposed in trials (have approval)]
        # Blue: IndicationMaxPhase and MaxPhase = 4               [All drugs (both those being repurposed in trials and not)]

        for drug in drug_dict.keys():
            max_phase = drug_dict[drug]['Max Phase']
            indications = list(drug_dict[drug].keys())[1:]
            for ind in indications:
                if max_phase < 4 and drug_dict[drug][ind][0] < 4:
                    listplot_data_variables[ind]['Value2'] += 1
                elif max_phase == 4 and drug_dict[drug][ind][0] < 4:
                    listplot_data_variables[ind]['Value4'] += 1
                elif max_phase == 4 and drug_dict[drug][ind][0] == 4:
                    listplot_data_variables[ind]['Value6'] += 1

        counter = 1
        labeled_listdata = {}
        for key in listdata.keys():
            new_key = f"{counter:02d} {key}"
            labeled_listdata[new_key] = listdata[key]
            counter += 1

        #Manually purging irrelevant classes
        to_purge = ['23 External causes of morbidity or mortality',
                    # 'Injury, poisoning or certain other consequences of external causes',
                    '26 Supplementary Chapter Traditional Medicine Conditions',
                    '27 Supplementary section for functioning assessment',
                    '24 Factors influencing health status or contact with health services',
                    '28 Extension Codes']

        for item in to_purge:
            del labeled_listdata[item]

        List_Data_result = {"category_array": [], "final_array": []}
        for key in labeled_listdata.keys():
            List_Data_result['category_array'].append('ReceptorFamily')
            List_Data_result['final_array'].append(key)
            for value in labeled_listdata[key]:
                List_Data_result['category_array'].append('Receptor')
                List_Data_result['final_array'].append(value)

        ####### TESTING ########

        records = []
        for indi in all_indis:
            records.append({'title':indi.title, 'level':indi.level})
        # Setup the master nested dict
        hierarchy = {}
        # Start with a "virtual root" at level -1
        stack = [(-1, hierarchy)]
        for record in records:
            level = record['level']
            title = record['title']
            # Pop until the top of the stack is strictly less than the current record's level
            while stack and stack[-1][0] >= level:
                stack.pop()
            # Now the stack top should be the parent level
            parent_dict = stack[-1][1]
            # Create a new entry for the current record
            parent_dict[title] = {}
            # Push current node
            stack.append((level, parent_dict[title]))

        counter = 1
        labeled_hierarchy = {}
        for key in hierarchy.keys():
            new_key = f"{counter:02d} {key}"
            labeled_hierarchy[new_key] = hierarchy[key]
            counter += 1

        for item in to_purge:
            del labeled_hierarchy[item]

        drug_dict_leaf = {}
        for item in drug_data:
            if item.ligand.name not in drug_dict_leaf.keys():
                drug_dict_leaf[item.ligand.name] = {'Max Phase': 0}
            if item.indication.title not in drug_dict_leaf[item.ligand.name].keys():
                drug_dict_leaf[item.ligand.name][item.indication.title] = [item.indication_max_phase, item.indication.title]
            if item.indication_max_phase > drug_dict_leaf[item.ligand.name]['Max Phase']:
                drug_dict_leaf[item.ligand.name]['Max Phase'] = item.indication_max_phase

        #now need to populated them
        for drug in drug_dict_leaf.keys():
            max_phase = drug_dict_leaf[drug]['Max Phase']
            indications = list(drug_dict_leaf[drug].keys())[1:]
            for ind in indications:
                if max_phase < 4 and drug_dict_leaf[drug][ind][0] < 4:
                    increment_color_value(labeled_hierarchy, ind, 'red')
                elif max_phase == 4 and drug_dict_leaf[drug][ind][0] < 4:
                    increment_color_value(labeled_hierarchy, ind, 'purple')
                elif max_phase == 4 and drug_dict_leaf[drug][ind][0] == 4:
                    increment_color_value(labeled_hierarchy, ind, 'blue')

        max_labels = {}
        for primary_key in labeled_hierarchy:
            r_val, r_k, p_val, p_k, b_val, b_k = find_max_values(labeled_hierarchy[primary_key], primary_key)
            max_labels[primary_key] = {
                'red': r_k if r_val != 0 else '#ffffff',
                'purple': p_k if p_val != 0 else '#ffffff',
                'blue': b_k if b_val != 0 else '#ffffff'
            }

        propagate_max_colors(labeled_hierarchy)
        colored = convert_dict_to_colors(labeled_hierarchy)

        # Now we update max_labels with actual colors from colored_result.
        for major_key, color_dict in max_labels.items():
            for color_key, node_name in color_dict.items():
                if node_name == '#ffffff':
                    continue

                # Find the node in colored_result that matches node_name
                node = find_node_by_name(colored, node_name)
                if node is not None and color_key in node:
                    # Replace the string in max_labels with the actual color tuple
                    max_labels[major_key][color_key] = node[color_key]

        del colored['red']
        del colored['purple']
        del colored['blue']
        # Combine data and labels into a single structure:
        combined = []
        for main_key, keys in colored.items():
            main_label = max_labels.get(main_key, {"red": "#ffffff", "purple": "#ffffff", "blue": "#ffffff"})
            rendered_html = render_nested_structure(keys)
            combined.append({
                "main_key": main_key,
                "main_label": main_label,
                "rendered_items": rendered_html
            })

        context['data'] = json.dumps(indication_data)
        context['listplot_data'] = json.dumps(labeled_listdata)
        context['listplot_data_variables'] = json.dumps(listplot_data_variables)
        context['listplot_datatypes'] = json.dumps(data_types_list)
        context['Label_Conversion'] = json.dumps(Label_Conversion)
        context['List_Data_result'] = json.dumps(List_Data_result)
        context['combined'] = combined
        context['data'] = colored
        context['labels'] = max_labels

        return context

class TargetSelectionTool(TemplateView):
    # Template using this class #
    template_name = 'TargetSelectionTool.html'
    # Get context for hmtl usage #
    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)

        # Fetch all proteins
        all_proteins = Protein.objects.filter(
            species_id=1,
            parent_id__isnull=True,
            accession__isnull=False,
            family_id__slug__startswith='0'
        ).exclude(
            family_id__slug__startswith='007'
        ).exclude(
            family_id__slug__startswith='008'
        ).annotate(
            #Fetch single gene name and entrez_id for each target using subqueries, prioritizing lowest entrez_id
            gene_name=Subquery(
                Gene.objects.filter(proteins=OuterRef('pk')).order_by('entrez_id').values('name')[:1]
            ),
            gene_entrez_id=Subquery(
                Gene.objects.filter(proteins=OuterRef('pk')).order_by('entrez_id').values('entrez_id')[:1]
            )
        ).values(
            'id', 'entry_name', 'gene_name', 'gene_entrez_id', 'name', 'family__parent__name',
            'family__parent__parent__name', 'family__parent__parent__parent__name'
        )

        # Convert to DataFrame
        proteins_df = pd.DataFrame(list(all_proteins))
        proteins_df.rename(columns={
            'id': 'Target ID',
            'entry_name': 'Uniprot entry',
            'gene_name': 'Gene name',
            'gene_entrez_id': 'Gene_entrez_id',
            'name': 'Protein name',
            'family__parent__name': 'Receptor family',
            'family__parent__parent__name': 'Ligand type',
            'family__parent__parent__parent__name': 'Class'
        }, inplace=True)

        #Remove duplicates due to multiple gene synonyms on some target entries
        proteins_df = proteins_df.groupby(by=['Target ID'], as_index = False).first()

        #Convert Gene_entrez_ids to web links
        entrez_websource = WebResource.objects.get(slug="entrez_gene")
        proteins_df['Gene_entrez_weblink'] = proteins_df['Gene_entrez_id'].apply(lambda id: str(WebLink(index=id, web_resource=entrez_websource)) if id != "" else "")
        proteins_df.drop(columns=['Gene_entrez_id'], inplace=True)

        # Fetch all data in a single query
        table_data = Drugs.objects.select_related(
            'target__family__parent__parent__parent',
            'moa',
            'ligand__ligand_type'
        ).values(
            'target',
            'target__entry_name',
            'target__name',
            'target__family__parent__name',
            'target__family__parent__parent__name',
            'target__family__parent__parent__parent__name',
            'ligand',
            'ligand__ligand_type',
            'novelty_score',
            'drug_status',
            'moa__name',
            'indication_max_phase',
            'publication_count',
            'target_level'
        )

        # Convert the table_data queryset to a list of dictionaries
        table_data_list = list(table_data)

        # Convert the list of dictionaries to a pandas DataFrame
        df = pd.DataFrame(table_data_list)

        # Drop duplicates based on all columns
        # df.drop_duplicates(inplace=True)

        # Rename the columns to your desired format
        df.rename(columns={
            'target': 'Target ID', # Target ID
            'target__entry_name': "Uniprot entry", # Uniprot entry name
            'target__name': 'Protein name', # Protein name
            'target__family__parent__name': 'Receptor family', # Receptor family
            'target__family__parent__parent__name': 'Ligand type', # Ligand type
            'target__family__parent__parent__parent__name': 'Class', # Class
            'ligand': "Ligand ID", # Ligand ID
            'ligand__ligand_type': "LigandType", # Ligand type
            'novelty_score': 'Novelty (Pharos)', #Only novelty we have
            'drug_status': 'Status', #Approval
            'moa__name': 'Mode of action', #Modality
            'indication_max_phase': 'Phase', # phase
            'publication_count': 'Literature', # Pub count
            'target_level': 'IDG'
            }, inplace=True)

        # Merge the proteins DataFrame with the drugs DataFrame on 'Target ID'
        df = proteins_df.merge(df, on='Target ID', how='left', suffixes=('', '_dup'))

        # Drop duplicate columns introduced by the merge
        df = df.loc[:, ~df.columns.str.endswith('_dup')]

        # Fill missing values in the merged DataFrame with an empty string
        df.fillna("", inplace=True)

        # Assuming `target_ids` is a list of IDs from the Drugs data
        target_ids = df['Target ID'].unique()

        # Fetch only relevant `Structure` data for matching `target_ids`
        structure_data = Structure.objects.filter(
            protein_conformation__protein__parent__id__in=target_ids
        ).exclude(structure_type__slug__startswith='af-').values(
            'protein_conformation__protein__parent__id'
        ).annotate(
            Inactive=Count(Case(When(state_id=1, then=1), output_field=IntegerField())),
            Active=Count(Case(When(state_id=2, then=1), output_field=IntegerField())),
            Intermediate=Count(Case(When(state_id=3, then=1), output_field=IntegerField())),
            Total=Count('id')
        ).order_by('protein_conformation__protein__parent__id')

        # Convert structure data to a DataFrame
        structure_df = pd.DataFrame(list(structure_data))

        # Rename columns for better clarity
        structure_df.rename(columns={
            'protein_conformation__protein__parent__id': 'Target ID'
        }, inplace=True)
        # Merge structure counts into the drugs DataFrame on `Target ID`
        df = df.merge(structure_df, on='Target ID', how='left')

        # Fill any missing values for states with 0 (in case a target ID has no structure data)
        df[['Active', 'Inactive', 'Intermediate', 'Total']] = df[['Active', 'Inactive', 'Intermediate', 'Total']].fillna(0)

        # Define stimulatory and inhibitory MOA categories
        stim_moa = ['Partial agonist', 'Agonist', 'PAM']
        inhib_moa = ['Antagonist', 'Inverse agonist', 'NAM']

        # Precompute classifications: Label as 'Drug' if status is 'Approved', otherwise as 'Agent'
        df['Classification'] = df['Status'].apply(lambda x: 'Drug' if x == 'Approved' else 'Agent')

        # Filter data for stimulatory and inhibitory based on Mode of Action
        stim_moa = ['Partial agonist', 'Agonist', 'PAM']
        inhib_moa = ['Antagonist', 'Inverse agonist', 'NAM']

        stim_moa_df = df[df['Mode of action'].isin(stim_moa)]
        inhib_moa_df = df[df['Mode of action'].isin(inhib_moa)]

        # Define a helper function to compute the max phase safely
        def get_max_phase(df, group_cols):
            if df.empty:
                return pd.Series(dtype='float64')  # Return an empty Series if DataFrame is empty
            return df.groupby(group_cols)['Phase'].max()

        # Define a helper function to compute unique counts for drugs/agents
        def compute_unique_counts(df, group_cols, value_col, classification_col):
            """
            Aggregates the unique count of `value_col` grouped by `group_cols`
            and splits them by their classification.
            """
            grouped = df.groupby(group_cols)[[value_col, classification_col]].apply(
                lambda x: pd.Series({
                    'Drug': len(set(x.loc[x[classification_col] == 'Drug', value_col])),
                    'Agent': len(set(x.loc[x[classification_col] == 'Agent', value_col]))
                })
            )
            return grouped.reset_index()

        # Group columns for aggregation
        group_cols = ['Target ID']

        # Compute All_Max_Phase for each Target ID as the maximum phase across all drugs and agents
        all_max_phase = df.groupby(group_cols)['Phase'].max().reset_index().rename(columns={'Phase': 'All_Max_Phase'})

        # Compute stimulatory max phases
        stim_max_phase = get_max_phase(stim_moa_df, group_cols).rename('Stimulatory_max_phase')

        # Compute inhibitory max phases
        inhib_max_phase = get_max_phase(inhib_moa_df, group_cols).rename('Inhibitory_max_phase')

        # Compute unique stimulatory counts
        stim_counts_unique = compute_unique_counts(
            stim_moa_df, group_cols, 'Ligand ID', 'Classification'
        ).rename(columns={'Drug': 'Stimulatory_Drugs', 'Agent': 'Stimulatory_Agents'})

        # Compute unique inhibitory counts
        inhib_counts_unique = compute_unique_counts(
            inhib_moa_df, group_cols, 'Ligand ID', 'Classification'
        ).rename(columns={'Drug': 'Inhibitory_Drugs', 'Agent': 'Inhibitory_Agents'})

        # add in the ligand type aggragation

        # Filter only relevant ligand types
        ligand_type_df = df[df['LigandType'].isin([1, 2, 3, 8])]

        # Drop duplicates to ensure we're only counting each ligand once per target
        unique_ligands = ligand_type_df[['Target ID', 'Ligand ID', 'LigandType']].drop_duplicates()

        # Create boolean masks
        # small_mol = unique_ligands['LigandType'] == 1
        # biologics = unique_ligands['LigandType'].isin([2, 3, 8])

        # Count unique ligand IDs by type per target
        # ligand_type_counts = unique_ligands.groupby('Target ID').apply(
        #     lambda g: pd.Series({
        #         'Type_small_molecule': g[small_mol & (g['Target ID'] == g.name)]['Ligand ID'].nunique(),
        #         'Type_biologics': g[biologics & (g['Target ID'] == g.name)]['Ligand ID'].nunique()
        #     })
        # ).reset_index()

        # Group and count unique ligands of each type per Target ID
        ligand_type_counts = unique_ligands.groupby('Target ID').apply(lambda group: pd.Series({
            'Type_small_molecule': group[group['LigandType'] == 1]['Ligand ID'].nunique(),
            'Type_biologics': group[group['LigandType'].isin([2, 3, 8])]['Ligand ID'].nunique(),
        })).reset_index()

        # Merge these into the aggregated data
        agg_data_targets = all_max_phase.merge(stim_max_phase, on=group_cols, how='left')
        agg_data_targets = agg_data_targets.merge(inhib_max_phase, on=group_cols, how='left')
        agg_data_targets = agg_data_targets.merge(stim_counts_unique, on=group_cols, how='left').fillna(0)
        agg_data_targets = agg_data_targets.merge(inhib_counts_unique, on=group_cols, how='left').fillna(0)
        agg_data_targets = agg_data_targets.merge(ligand_type_counts, on=group_cols, how='left').fillna(0)

        # Compute total unique drugs and agents
        agg_data_targets['All_Drugs'] = agg_data_targets['Stimulatory_Drugs'] + agg_data_targets['Inhibitory_Drugs']
        agg_data_targets['All_Agents'] = agg_data_targets['Stimulatory_Agents'] + agg_data_targets['Inhibitory_Agents']

        # Ensure NaNs are replaced with 0 for all aggregated fields
        agg_data_targets.fillna({
            'Stimulatory_max_phase': 0, 'Inhibitory_max_phase': 0,
            'Stimulatory_Drugs': 0, 'Stimulatory_Agents': 0,
            'Inhibitory_Drugs': 0, 'Inhibitory_Agents': 0,
            'All_Drugs': 0, 'All_Agents': 0, 'Type_small_molecule': 0, 'Type_biologics': 0
        }, inplace=True)

        # Merge `agg_data_targets` back into the original DataFrame on `Target ID`
        df_first = df.merge(agg_data_targets, on=group_cols, how='left')

        # Keep only the specified columns in df_first
        keep_col_names = [
            'Target ID', 'Uniprot entry', 'Gene name', 'Gene_entrez_weblink',
            'Protein name', 'Receptor family', 'Ligand type',
            'Class', 'Literature', 'Novelty (Pharos)', 'IDG',
            'Total', 'Active', 'Inactive', 'All_Max_Phase', 'All_Drugs', 'All_Agents',
            'Stimulatory_max_phase', 'Stimulatory_Drugs', 'Stimulatory_Agents',
            'Inhibitory_max_phase', 'Inhibitory_Drugs', 'Inhibitory_Agents',
            'Type_small_molecule', 'Type_biologics']

        # Keep only the specified columns in df_first
        df_first = df_first[keep_col_names]

        # Drop duplicates based on the reduced set of columns
        df_first.drop_duplicates(inplace=True)

        # Fetch data from AssayExperiment with select_related for performance
        assay_data = (
            AssayExperiment.objects.select_related("protein", "ligand")
            .values("protein", "ligand", "value_type"))

        # Convert to DataFrame
        assay_df = pd.DataFrame(list(assay_data))

        # Step 1: Count occurrences of "pEC50" and "pIC50" for each target-ligand pair
        pair_counts = assay_df.groupby(["protein", "ligand", "value_type"]).size().unstack(fill_value=0)

        # Step 2: Identify and exclude pairs that have both "pEC50" and "pIC50"
        pairs_with_both = pair_counts[(pair_counts.get("pEC50", 0) > 0) & (pair_counts.get("pIC50", 0) > 0)].reset_index()
        remaining_pairs = assay_df[~assay_df.set_index(["protein", "ligand"]).index.isin(
            pairs_with_both.set_index(["protein", "ligand"]).index)]

        # Step 3: Count remaining unique target-ligand pairs and their "pEC50" and "pIC50"
        remaining_counts = (
            remaining_pairs.groupby(["protein", "ligand", "value_type"])
            .size()
            .unstack(fill_value=0)
            .reset_index()
        )

        # Step 4: Aggregate at the target level
        agg_results = remaining_counts.groupby("protein").agg(
            Total_Ligands=("ligand", "nunique"),
            pEC50_Count=("pEC50", "sum"),
            pIC50_Count=("pIC50", "sum"),
        ).reset_index()

        # Rename "protein" to "Target ID" to match df_first
        agg_results.rename(columns={"protein": "Target ID"}, inplace=True)

        # Merge the new aggregated data with df_first
        df_first = pd.merge(df_first, agg_results, on="Target ID", how="left")

        # Fill NaN values in the newly added columns with 0
        df_first[["Total_Ligands", "pEC50_Count", "pIC50_Count"]] = df_first[
            ["Total_Ligands", "pEC50_Count", "pIC50_Count"]
        ].fillna("")

        # Replace 0s with empty strings in relevant count columns
        df_first[["Total_Ligands", "pEC50_Count", "pIC50_Count",
                "Total", "Active", "Inactive", "All_Drugs", "All_Agents",
                "Stimulatory_max_phase", "Stimulatory_Drugs", "Stimulatory_Agents",
                "Inhibitory_max_phase", "Inhibitory_Drugs", "Inhibitory_Agents", "Type_small_molecule", "Type_biologics"]] = \
            df_first[["Total_Ligands", "pEC50_Count", "pIC50_Count",
                    "Total", "Active", "Inactive", "All_Drugs", "All_Agents",
                    "Stimulatory_max_phase", "Stimulatory_Drugs", "Stimulatory_Agents",
                    "Inhibitory_max_phase", "Inhibitory_Drugs", "Inhibitory_Agents", "Type_small_molecule", "Type_biologics"]].replace(0, "")

        # Step 1: Query the TissueExpression model
        tissue_datatable = TissueExpression.objects.select_related('tissue').values(
            'protein',        # Target ID
            'tissue__name',  # Tissue name
            'value'          # Expression value
        )

        # Step 2: Convert the QuerySet to a DataFrame
        df_tissue_table = pd.DataFrame(list(tissue_datatable))

        # Step 3: Rename columns for clarity
        df_tissue_table.rename(columns={
            'protein': 'Target ID',
            'tissue__name': 'Tissue Name',
            'value': 'Expression Value'
        }, inplace=True)

        # Step 4: Clean tissue names (e.g., replace underscores, capitalize)
        def clean_tissue_name(name):
            # Replace underscores with spaces, remove trailing "_X", and capitalize words
            name = name.replace("_", " ").rstrip(" 1234567890").title()
            return name

        df_tissue_table['Tissue Name'] = df_tissue_table['Tissue Name'].apply(clean_tissue_name)

        # Step 5: Pivot the DataFrame to make Tissue Names columns
        df_tissue_pivot = df_tissue_table.pivot_table(
            index='Target ID',         # Use Target ID as the index
            columns='Tissue Name',     # Use Tissue Name as columns
            values='Expression Value', # Use Expression Value as values
            aggfunc='first'            # If there are duplicates, take the first value
        )

        # Step 6: Reset index to turn Target ID back into a column
        df_tissue_pivot.reset_index(inplace=True)

        # Step 7: Fill NaN values with 0 (optional, based on your needs)
        df_tissue_pivot.fillna("", inplace=True)

        # Step 8: Generate a list of tissue column names (excluding Target ID)
        tissue_column_list = sorted([col for col in df_tissue_pivot.columns if col != 'Target ID'])

        # Merge df_first with df_tissue_pivot on 'Target ID'
        merged_df = pd.merge(df_first, df_tissue_pivot, on="Target ID", how="left")

        # Fill NaN values with empty strings for all columns
        merged_df.fillna("", inplace=True)

        # Add cancer to table 1

        # Add cancer data
        cancer_data = Protein.objects.select_related('Protein').values(
            'cancer__protein',
            'cancer__cancer__name',
            'cancer__expression__max_expression',
        )

        cancer_df = pd.DataFrame(list(cancer_data))
        cancer_df.rename(columns={
            'cancer__protein': 'Target ID',
            'cancer__cancer__name': 'Cancer',
            'cancer__expression__max_expression': 'Expression'
        }, inplace=True)

        # Filter only shared target IDs
        cancer_df = cancer_df[cancer_df['Target ID'].isin(target_ids)]

        # Ensure Target ID is an integer in cancer_df
        cancer_df['Target ID'] = cancer_df['Target ID'].astype(int)

       # Step 1: Calculate `NumberOfCancers` (sum of Medium and High for each Target ID)
        cancer_df['IsRelevantExpression'] = cancer_df['Expression'].isin(['Medium', 'High']).astype(int)
        number_of_cancers = (
            cancer_df.groupby('Target ID')['IsRelevantExpression']
            .sum()
            .reset_index(name='NumberOfCancers')
        )

        # Step 2: Pivot the cancer data
        cancer_pivot = cancer_df.pivot_table(
            index='Target ID',  # Rows: Target ID
            columns='Cancer',   # Columns: Cancer types
            values='Expression',  # Values: Expression levels
            aggfunc='first'      # Use the first value (only one row per Target ID + Cancer)
        )

        # Capitalize Cancer column names
        cancer_pivot.columns = [col.capitalize() for col in cancer_pivot.columns]

        # Reset index to make Target ID a column
        cancer_pivot.reset_index(inplace=True)

        # Step 3: Merge `NumberOfCancers` into the pivot table
        final_cancer_table = pd.merge(
            cancer_pivot,
            number_of_cancers,
            on='Target ID',
            how='left'
        )

        # Step 4: Merge cancer data with the main table
        Table1 = pd.merge(merged_df, final_cancer_table, on='Target ID', how='left')

        # Fill NaN values with empty strings
        Table1.fillna("", inplace=True)

        # Exclude 'Target ID' and sort the remaining columns
        cancer_column_list = sorted([col for col in cancer_pivot.columns if col != 'Target ID'])

        ########################################
        # Disease indications and associations #
        ########################################

        # Fetch the second data table with indication and master level
        table_data_2 = Drugs.objects.select_related(
            'target__family__parent__parent__parent',
            'moa',
            'indication'
        ).annotate(
            #Fetch single gene name and entrez_id for each target using subqueries, prioritizing lowest entrez_id
            gene_name=Subquery(
                Gene.objects.filter(proteins=OuterRef('target__pk')).order_by('entrez_id').values('name')[:1]
            ),
            gene_entrez_id=Subquery(
                Gene.objects.filter(proteins=OuterRef('target__pk')).order_by('entrez_id').values('entrez_id')[:1]
            )
        ).values(
            'target',  # Target ID
            'target__entry_name', # Uniprot entry name
            'target__genes__name', # Gene name/symbol
            'target__genes__entrez_id', # Gene NCBI/Entrez ID
            'target__name',  # Protein name
            'target__family__parent__name',  # Receptor family
            'target__family__parent__parent__name',  # Ligand type
            'target__family__parent__parent__parent__name',  # Class
            'ligand',  # Ligand ID
            'indication__title',  # Indication title
            'indication__slug',  # Indication slug
            'indication__code',
            'disease_association__association_score',  # Association score
            'drug_status'  # Status
        )

        # Convert to DataFrame
        df_table_2 = pd.DataFrame(list(table_data_2))
        df_table_2.rename(columns={
            'target': 'Target ID',
            'target__entry_name': "Uniprot entry",
            'target__genes__name': "Gene name",
            'target__genes__entrez_id': "Gene_entrez_id",
            'target__name': 'Protein name',
            'target__family__parent__name': 'Receptor family',
            'target__family__parent__parent__name': 'Ligand type',
            'target__family__parent__parent__parent__name': 'Class',
            'ligand': "Ligand ID",
            'indication__title': 'Indication',
            'indication__slug': 'Indication Slug',
            'indication__code': "ICD11",
            'disease_association__association_score': 'Association Score',
            'drug_status': 'Status'
        }, inplace=True)

        #Convert Gene_entrez_ids to web links
        df_table_2['Gene_entrez_weblink'] = df_table_2['Gene_entrez_id'].apply(lambda id: str(WebLink(index=id, web_resource=entrez_websource)) if id != "" else "")
        df_table_2.drop(columns=['Gene_entrez_id'], inplace=True)

        # Classify based on status
        df_table_2['Classification'] = df_table_2['Status'].apply(
            lambda x: 'Drug' if x == 'Approved' else 'Agent'
        )

        # Fetch top-level indications (master indications)
        master_indications = Indication.objects.filter(level=0).values('slug', 'title')

        # Create a mapping dictionary for master indications
        master_mapping = {entry['slug']: entry['title'] for entry in master_indications}

        # Map master titles to indications based on the top-level slug
        df_table_2['Master Indication'] = df_table_2['Indication Slug'].apply(
            lambda slug: master_mapping.get(slug.split("_")[0], "")
        )

        # Fetch ATC data
        atc_data = ATCCodes.objects.values(
            'ligand', 'code__index', 'name_0', 'name_1'
        )

        # Convert to DataFrame
        atc_df = pd.DataFrame(list(atc_data))
        atc_df.rename(columns={
            'ligand': 'Ligand ID',
            'code__index': 'ATC Code',
            'name_0': 'ATC Parent Name',
            'name_1': 'ATC Name'
        }, inplace=True)

        # Merge ATC data with the second table
        df_table_2 = df_table_2.merge(atc_df, on='Ligand ID', how='left')

        # Fill missing values with empty strings
        df_table_2.fillna("", inplace=True)

        # Aggregate the table efficiently
        def aggregate_table_fast(df):
            # Select unique rows for non-aggregated columns
            unique_df = df.drop_duplicates(subset=[
                'Target ID', 'Uniprot entry', 'Gene name', 'Protein name', 'Receptor family',
                'Ligand type', 'Class', 'Master Indication', 'Indication',
                'ATC Code', 'ATC Parent Name', 'ATC Name', 'Association Score','ICD11'
            ])

            # Compute counts for 'Drug' and 'Agent' using crosstab
            classification_counts = pd.crosstab(
                index=[df['Target ID'], df['Indication']],  # Grouping keys
                columns=df['Classification']
            ).reset_index()

            # Rename columns for clarity
            classification_counts.rename(columns={
                'Drug': 'Drug Count',
                'Agent': 'Agent Count'
            }, inplace=True)

            # Merge classification counts back into the unique rows
            result = pd.merge(
                unique_df,
                classification_counts,
                on=['Target ID', 'Indication'],
                how='left'
            )

            # Fill NaN counts with 0
            result[['Drug Count', 'Agent Count']] = result[['Drug Count', 'Agent Count']].fillna(0).astype(int)

            return result

        # Aggregate the table to collapse unique rows and compute counts quickly
        df_table_2 = aggregate_table_fast(df_table_2)

        def remove_duplicates_and_fill_gaps(df):
            # Replace empty strings with NaN for better handling of missing values
            df.replace("", pd.NA, inplace=True)

            # Sort by completeness (count of non-missing values) and drop exact duplicates
            df = df.sort_values(by=df.columns.to_list(), key=lambda col: col.notna(), ascending=False)
            df = df.drop_duplicates(subset=df.columns.difference(['ATC Code', 'ATC Name', 'ATC Parent Name']), keep='first')

            # Fill gaps by replacing NaN back with empty strings for display purposes
            df.fillna("", inplace=True)

            return df
        # Remove duplicates and prioritize completeness
        df_table_2 = remove_duplicates_and_fill_gaps(df_table_2)

        # Retain only the specified columns
        keep_cols = [
            'Target ID', 'Uniprot entry', 'Gene name', 'Gene_entrez_weblink', 'Protein name', 'Receptor family', 'Ligand type',
            'Class', 'Master Indication', 'Indication', 'ATC Code',
            'ATC Parent Name', 'ATC Name', 'Association Score', 'Drug Count', 'Agent Count','ICD11'
        ]
        df_table_2 = df_table_2[keep_cols]

        # Replace 0s with empty strings for display purposes
        df_table_2.replace(0, "", inplace=True)

        ####################
        # Untapped targets #
        ####################

        # Extract unique Target IDs from df_table_2
        existing_target_ids = df_table_2['Target ID'].unique()

        # Filter table_data_3 to exclude existing Target IDs
        filtered_table_data_3 = IndicationAssociation.objects.select_related(
            'target__family__parent__parent__parent',
            'indication'
        ).exclude(
            target__in=existing_target_ids  # Exclude targets that already exist in df_table_2
        ).annotate(
            #Fetch single gene name and entrez_id for each target using subqueries, prioritizing lowest entrez_id
            gene_name=Subquery(
                Gene.objects.filter(proteins=OuterRef('target__pk')).order_by('entrez_id').values('name')[:1]
            ),
            gene_entrez_id=Subquery(
                Gene.objects.filter(proteins=OuterRef('target__pk')).order_by('entrez_id').values('entrez_id')[:1]
            )
        ).values(
            'target',  # Target ID
            'target__entry_name',  # Uniprot entry
            'gene_name', # Gene symbol/name
            'gene_entrez_id', # Gene NCBI/entrez ID
            'target__name',  # Protein name
            'target__family__parent__name',  # Receptor family
            'target__family__parent__parent__name',  # Ligand type
            'target__family__parent__parent__parent__name',  # Class
            'indication__title',  # Indication title
            'indication__slug', # Indication slug
            'association_score',  # Association score
            'indication__code' #ICD11
        )

        # Convert the filtered QuerySet to a DataFrame
        df_table_3 = pd.DataFrame(list(filtered_table_data_3))

        df_table_3.rename(columns={
            'target': 'Target ID',
            'target__entry_name': 'Uniprot entry',
            'gene_name' : "Gene name",
            'gene_entrez_id' : "Gene_entrez_id", # Gene NCBI/entrez ID
            'target__name': 'Protein name',
            'target__family__parent__name': 'Receptor family',
            'target__family__parent__parent__name': 'Ligand type',
            'target__family__parent__parent__parent__name': 'Class',
            'indication__title': 'Indication',
            'indication__slug': 'Indication Slug',
            'association_score': 'Association Score',
            'indication__code': "ICD11"
        }, inplace=True)

        # Fetch master indications and create mapping
        master_indications = Indication.objects.filter(level=0).values('slug', 'title')
        master_mapping = {entry['slug']: entry['title'] for entry in master_indications}

        # Add Master Indication to df_table_3
        df_table_3['Master Indication'] = df_table_3['Indication Slug'].apply(
            lambda slug: master_mapping.get(slug.split("_")[0], "")
        )

        # Convert the final DataFrame to JSON
        json_records_targets = Table1.to_json(orient='records')
        json_records_table_2 = df_table_2.to_json(orient='records')
        json_records_table_3 = df_table_3.to_json(orient='records')
        context['Full_data'] = json_records_targets
        context['Disease_table'] = json_records_table_2
        context['Untapped_table'] = json_records_table_3
        context['Tissue_cols'] = tissue_column_list
        context['cancer_cols'] = cancer_column_list

        return context

@cache_page(60 * 60 * 24 * 28)
def drugmapping(request):
    context = dict()

    families = ProteinFamily.objects.all()
    lookup = {}
    for f in families:
        lookup[f.slug] = f.name.replace("receptors","").replace(" receptor","").replace(" hormone","").replace("/neuropeptide","/").replace(" (G protein-coupled)","").replace(" factor","").replace(" (LPA)","").replace(" (S1P)","").replace("GPR18, GPR55 and GPR119","GPR18/55/119").replace("-releasing","").replace(" peptide","").replace(" and oxytocin","/Oxytocin").replace("Adhesion class orphans","Adhesion orphans").replace("muscarinic","musc.").replace("-concentrating","-conc.")

    class_proteins = Protein.objects.filter(family__slug__startswith="0",source__name='SWISSPROT', species_id=1).prefetch_related('family').order_by('family__slug')

    temp = OrderedDict([
                    ('name',''),
                    ('trials', 0),
                    ('maxphase', 0),
                    ('approved', 0),
                    ('family_sum_approved', 0),
                    ('family_sum_trials' , 0),
                    ('establishment', 2),
                    ('children', OrderedDict())
                    ])

    coverage = OrderedDict()

    # Make the scaffold
    for p in class_proteins:
        #print(p,p.family.slug)
        fid = p.family.slug.split("_")
        if fid[0] not in coverage:
            coverage[fid[0]] = deepcopy(temp)
            coverage[fid[0]]['name'] = lookup[fid[0]]
        if fid[1] not in coverage[fid[0]]['children']:
            coverage[fid[0]]['children'][fid[1]] = deepcopy(temp)
            coverage[fid[0]]['children'][fid[1]]['name'] = lookup[fid[0]+"_"+fid[1]]
        if fid[2] not in coverage[fid[0]]['children'][fid[1]]['children']:
            coverage[fid[0]]['children'][fid[1]]['children'][fid[2]] = deepcopy(temp)
            coverage[fid[0]]['children'][fid[1]]['children'][fid[2]]['name'] = lookup[fid[0]+"_"+fid[1]+"_"+fid[2]][:28]
        if fid[3] not in coverage[fid[0]]['children'][fid[1]]['children'][fid[2]]['children']:
            coverage[fid[0]]['children'][fid[1]]['children'][fid[2]]['children'][fid[3]] = deepcopy(temp)
            coverage[fid[0]]['children'][fid[1]]['children'][fid[2]]['children'][fid[3]]['name'] = p.entry_name.split("_")[0] #[:10]

    # # POULATE WITH DATA
    total_approved = 0
    drugtargets_approved_class = Protein.objects.filter(drugs__status='approved').values('family_id__parent__parent__parent__slug').annotate(value=Count('drugs__name', distinct = True))
    for i in drugtargets_approved_class:
        fid = i['family_id__parent__parent__parent__slug'].split("_")
        coverage[fid[0]]['family_sum_approved'] += i['value']
        total_approved += i['value']

    drugtargets_approved_type = Protein.objects.filter(drugs__status='approved').values('family_id__parent__parent__slug').annotate(value=Count('drugs__name', distinct = True))
    for i in drugtargets_approved_type:
        fid = i['family_id__parent__parent__slug'].split("_")
        coverage[fid[0]]['children'][fid[1]]['family_sum_approved'] += i['value']

    drugtargets_approved_family = Protein.objects.filter(drugs__status='approved').values('family_id__parent__slug').annotate(value=Count('drugs__name', distinct = True))
    for i in drugtargets_approved_family:
        fid = i['family_id__parent__slug'].split("_")
        coverage[fid[0]]['children'][fid[1]]['children'][fid[2]]['family_sum_approved'] += i['value']
    drugtargets_approved_target = Protein.objects.filter(drugs__status='approved').values('family_id__slug').annotate(value=Count('drugs__name', distinct = True)).annotate(maxphase=Max('drugs__phase'))
    for i in drugtargets_approved_target:
        fid = i['family_id__slug'].split("_")
        coverage[fid[0]]['children'][fid[1]]['children'][fid[2]]['children'][fid[3]]['approved'] += i['value']
        if int(i['maxphase']) > coverage[fid[0]]['children'][fid[1]]['children'][fid[2]]['children'][fid[3]]['maxphase']:
            coverage[fid[0]]['children'][fid[1]]['children'][fid[2]]['children'][fid[3]]['maxphase'] = int(i['maxphase'])
        if i['value'] > 0:
            coverage[fid[0]]['children'][fid[1]]['children'][fid[2]]['children'][fid[3]]['establishment'] = 4

    total_trials = 0
    drugtargets_trials_class = Protein.objects.filter(drugs__status__in=['in trial']).values('family_id__parent__parent__parent__slug').annotate(value=Count('drugs__name', distinct = True))
    for i in drugtargets_trials_class:
        fid = i['family_id__parent__parent__parent__slug'].split("_")
        coverage[fid[0]]['family_sum_trials'] += i['value']
        total_trials += i['value']

    drugtargets_trials_type = Protein.objects.filter(drugs__status__in=['in trial']).values('family_id__parent__parent__slug').annotate(value=Count('drugs__name', distinct = True))
    for i in drugtargets_trials_type:
        fid = i['family_id__parent__parent__slug'].split("_")
        coverage[fid[0]]['children'][fid[1]]['family_sum_trials'] += i['value']

    drugtargets_trials_family = Protein.objects.filter(drugs__status__in=['in trial']).values('family_id__parent__slug').annotate(value=Count('drugs__name', distinct = True))
    for i in drugtargets_trials_family:
        fid = i['family_id__parent__slug'].split("_")
        coverage[fid[0]]['children'][fid[1]]['children'][fid[2]]['family_sum_trials'] += i['value']

    drugtargets_trials_target = Protein.objects.filter(drugs__status__in=['in trial']).values('family_id__slug').annotate(value=Count('drugs__name', distinct = True)).annotate(maxphase=Max('drugs__phase'))
    # add highest reached trial here
    for i in drugtargets_trials_target:
        fid = i['family_id__slug'].split("_")
        coverage[fid[0]]['children'][fid[1]]['children'][fid[2]]['children'][fid[3]]['trials'] += i['value']
        if int(i['maxphase']) > coverage[fid[0]]['children'][fid[1]]['children'][fid[2]]['children'][fid[3]]['maxphase']:
            coverage[fid[0]]['children'][fid[1]]['children'][fid[2]]['children'][fid[3]]['maxphase'] = int(i['maxphase'])
        if i['value'] > 0 and coverage[fid[0]]['children'][fid[1]]['children'][fid[2]]['children'][fid[3]]['establishment'] == 2:
            coverage[fid[0]]['children'][fid[1]]['children'][fid[2]]['children'][fid[3]]['establishment'] = 7

    # MAKE THE TREE
    tree = OrderedDict({'name':'GPCRome', 'family_sum_approved': total_approved, 'family_sum_trials': total_trials,'children':[]})
    i = 0
    n = 0
    for c,c_v in coverage.items():
        c_v['name'] = c_v['name'].split("(")[0]
        if c_v['name'].strip() == 'Other GPCRs':
            # i += 1
            continue
            # pass
        children = []
        for lt,lt_v in c_v['children'].items():
            children_rf = []
            for rf,rf_v in lt_v['children'].items():
                rf_v['name'] = rf_v['name'].split("<")[0]
                # if rf_v['name'].strip() == 'Taste 2':
                    # continue
                children_r = []
                for r,r_v in rf_v['children'].items():
                    r_v['sort'] = n
                    children_r.append(r_v)
                    n += 1
                rf_v['children'] = children_r
                rf_v['sort'] = n
                children_rf.append(rf_v)
            lt_v['children'] = children_rf
            lt_v['sort'] = n
            children.append(lt_v)
        c_v['children'] = children
        c_v['sort'] = n
        tree['children'].append(c_v)
        i += 1

    jsontree = json.dumps(tree)

    context["drugdata"] = jsontree

    return render(request, 'drugmapping.html', {'drugdata':context})

# @cache_page(60 * 60 * 24 * 28)
def indication_detail(request, code):

    code = code.upper()
    context = dict()
    #code = '4A8Z'
    indication_data = Drugs.objects.filter(indication__code=code).prefetch_related('ligand',
                                                                                        'target',
                                                                                        'indication',
                                                                                        'indication__uri')

    indication_name = Indication.objects.filter(code=code).values_list('title', flat=True).distinct()[0]

    sankey = {"nodes": [],
              "links": []}
    caches = {'indication':[],
              'level_0': [],
              'ligands': [],
              'targets': [],
              'entries': []}

    node_counter = 0
    # Initialize the path source matrix
    path_matrix = []
    link_id = 0  # Unique identifier for each link
    row_id = 0
    for record in indication_data:
        #assess the values for indication/ligand/protein
        indication_code = record.indication.title.capitalize()
        indication_uri = record.indication.uri.index if record.indication.uri else ''
        ligand_name = record.ligand.name.capitalize()
        indication_0 = record.indication.get_level_0().title
        uri = record.indication.uri.index if record.indication.uri else ''
        ligand_id = record.ligand.id
        protein_name = record.target.name
        target_name = record.target.entry_name
        #check for each value if it exists and retrieve the source node value
        if indication_code not in caches['indication']:
            sankey['nodes'].append({"node": node_counter, "name": indication_code, "url":'https://icd.who.int/browse/2024-01/mms/en#'+indication_uri,"column":"x4"})
            node_counter += 1
            caches['indication'].append(indication_code)
        indi_node = next((item['node'] for item in sankey['nodes'] if item['name'] == indication_code), None)

        if indication_0 not in caches['level_0']:
            sankey['nodes'].append({"node": node_counter, "name": indication_0, "url":'https://icd.who.int/browse/2024-01/mms/en#'+uri,"column":"x3"})
            node_counter += 1
            caches['level_0'].append(indication_0)
        level_0_node = next((item['node'] for item in sankey['nodes'] if item['name'] == indication_0), None)

        if [ligand_name, ligand_id] not in caches['ligands']:
            sankey['nodes'].append({"node": node_counter, "name": ligand_name, "url":'/ligand/'+str(ligand_id)+'/group',"column":"x2"})
            node_counter += 1
            caches['ligands'].append([ligand_name, ligand_id])
        lig_node = next((item['node'] for item in sankey['nodes'] if item['name'] == ligand_name), None)

        if protein_name not in caches['targets']:
            sankey['nodes'].append({"node": node_counter, "name": protein_name, "url":'/protein/'+str(target_name),"column":"x1"})
            node_counter += 1
            caches['targets'].append(protein_name)
            caches['entries'].append(target_name)
        prot_node = next((item['node'] for item in sankey['nodes'] if item['name'] == protein_name), None)

        # create matrix row
        row = {
            'row_id': row_id,
            'x1': prot_node,
            'x2': lig_node,
            'x3': level_0_node,
            'x4': indi_node,
            'x1_name': protein_name,
            'x2_name': ligand_name,
            'x3_name': indication_0,
            'x4_name': indication_name
        }

        sankey['links'].append({"source": prot_node, "target": lig_node, "value": 1, "ligtrace": protein_name, "prottrace": indication_name, "linkage_key": "primary","link_identifier": link_id})  # x1 -> x2 (Protein → Ligand)
        link_id += 1
        sankey['links'].append({"source": lig_node, "target": level_0_node, "value": 1, "ligtrace": protein_name, "prottrace": indication_0, "linkage_key": "primary","link_identifier": link_id})  # x2 -> x3 (Ligand → Level 0)
        link_id += 1
        sankey['links'].append({"source": level_0_node, "target": indi_node, "value": 1, "ligtrace": protein_name, "prottrace": indication_name, "linkage_key": "primary","link_identifier": link_id})  # x3 -> x4 (Level 0 → Indication)
        link_id += 1
        path_matrix.append(row)
        row_id += 1

    #Fixing redundancy in sankey['links']
    unique_combinations = {}

    for d in sankey['links']:
        # Create a key based on source and target for identifying unique combinations
        key = (d['source'], d['target'])

        if key in unique_combinations:
            # If the combination exists, add the value to the existing entry
            unique_combinations[key]['value'] += d['value']
        else:
            # If it's a new combination, add it to the dictionary
            unique_combinations[key] = d

    # Convert the unique_combinations back to a list of dictionaries
    sankey['links'] = list(unique_combinations.values())

    # Find all nodes that are actually used in the links
    used_nodes = set()
    for link in sankey['links']:
        used_nodes.add(link['source'])
        used_nodes.add(link['target'])

    # Build new nodes list and a mapping from old to new node indices
    old_to_new = {}
    new_nodes = []
    new_index = 0
    for node_dict in sankey['nodes']:
        old_idx = node_dict["node"]
        if old_idx in used_nodes:
            # Only keep nodes that are actually used
            new_node_dict = {
                "node": new_index,
                "name": node_dict["name"],
                "url": node_dict["url"],
                "column": node_dict['column']
            }
            new_nodes.append(new_node_dict)
            old_to_new[old_idx] = new_index
            new_index += 1

    # Update the source and target indices in the links
    for link in sankey['links']:
        link["source"] = old_to_new[link["source"]]
        link["target"] = old_to_new[link["target"]]

    # Replace the nodes list with the new one
    sankey['nodes'] = new_nodes

    # Needed code for sankey in future updates (matrix solution):
    # def update_max_paths(sankey_nodes, path_matrix):

    #     # Function to count unique paths leading **to** the node
    #     def count_unique_paths_to(previous_col, current_col, node_name):
    #         unique_paths = set()
    #         # Filter rows to include only those leading to the current node
    #         filtered_rows = [row for row in path_matrix if row[current_col + '_name'] == node_name]
    #         for row in filtered_rows:
    #             from_node = row[previous_col]
    #             to_node = row[current_col]
    #             unique_paths.add((from_node, to_node))
    #         return len(unique_paths)

    #     # Function to count unique paths leading **from** the node
    #     def count_unique_paths_from(current_col, next_col, node_name):
    #         unique_paths = set()
    #         # Filter rows to include only those originating from the current node
    #         filtered_rows = [row for row in path_matrix if row[current_col + '_name'] == node_name]
    #         for row in filtered_rows:
    #             from_node = row[current_col]
    #             to_node = row[next_col]
    #             unique_paths.add((from_node, to_node))
    #         return len(unique_paths)

    #     # Iterate over all nodes and calculate max_paths
    #     for node in sankey_nodes:
    #         node_name = node['name']
    #         node_column = node['column']

    #         # Determine which transition to calculate based on column
    #         if node_column == 'x1':
    #             # Only count outgoing paths
    #             node['max_paths'] = count_unique_paths_from('x1', 'x2', node_name)
    #         elif node_column == 'x2':
    #             # Count both incoming and outgoing paths
    #             node['max_paths'] = max(
    #                 count_unique_paths_to('x1', 'x2', node_name),  # x1 -> x2
    #                 count_unique_paths_from('x2', 'x3', node_name) # x2 -> x3
    #             )
    #         elif node_column == 'x3':
    #             # Count both incoming and outgoing paths
    #             node['max_paths'] = max(
    #                 count_unique_paths_to('x2', 'x3', node_name),  # x2 -> x3
    #                 count_unique_paths_from('x3', 'x4', node_name) # x3 -> x4
    #             )
    #         elif node_column == 'x4':
    #             # Only count incoming paths
    #             node['max_paths'] = count_unique_paths_to('x3', 'x4', node_name)

    # # Call the function to update max_paths
    # update_max_paths(sankey['nodes'], path_matrix)

    total_points = len(caches['targets']) + len(caches['targets']) + 1
    if len(caches['ligands']) > len(caches['targets']):
        context['nodes_nr'] = len(caches['ligands'])
    else:
        context['nodes_nr'] = len(caches['targets'])
    context['indication_code'] = code
    context['indication_uri'] = indication_uri
    context['indication'] = indication_name.capitalize()
    context['sankey'] = json.dumps(sankey)
    context['path_matrix'] = json.dumps(path_matrix)
    context['points'] = total_points
    context['targets'] = list(caches['entries'])
    context['ligands'] = list(caches['ligands'])
    return render(request, 'indication_detail.html', context)

def fetch_sankey_indi_data_view(request):
    logger.info("fetch_sankey_data_view called")  # Log the function call
    code = request.GET.get('code', None)  # Extract `entry_name` from the request

    if not code:
        logger.error("No ICD11 provided in the request")
        return JsonResponse({'error': 'Missing ICD11'}, status=400)

    logger.info(f"Received ICD11: {code}")  # Log the entry name

    try:
        sankey_data = get_sankey_indi_data(code)  # Call your helper function
        if sankey_data:
            logger.info(f"Sankey data successfully retrieved for ICD11: {code}")
            return JsonResponse({'sankey_data': sankey_data}, status=200)
        else:
            logger.warning(f"No data found for ICD11: {code}")
            return JsonResponse({'error': 'No data found for this ICD11'}, status=404)
    except Exception:
        logger.exception(f"An error occurred while fetching sankey data for ICD11: {code}")
        return JsonResponse({'error': 'An internal server error occurred'}, status=500)

def get_sankey_indi_data(code):

    code = code.upper()

    indication_data = Drugs.objects.filter(indication__code=code).prefetch_related('ligand',
                                                                                        'target',
                                                                                        'indication',
                                                                                        'indication__uri')

    indication_name = Indication.objects.filter(code=code).values_list('title', flat=True).distinct()[0]

    sankey = {"nodes": [],
              "links": []}
    caches = {'indication':[],
              'level_0': [],
              'ligands': [],
              'targets': [],
              'entries': []}

    node_counter = 0
    # Initialize the path source matrix
    path_matrix = []
    link_id = 0  # Unique identifier for each link
    row_id = 0
    for record in indication_data:
        #assess the values for indication/ligand/protein
        indication_code = record.indication.title.capitalize()
        indication_uri = record.indication.uri.index if record.indication.uri else ''
        ligand_name = record.ligand.name.capitalize()
        indication_0 = record.indication.get_level_0().title
        uri = record.indication.uri.index if record.indication.uri else ''
        ligand_id = record.ligand.id
        protein_name = record.target.name
        target_name = record.target.entry_name
        #check for each value if it exists and retrieve the source node value
        if indication_code not in caches['indication']:
            sankey['nodes'].append({"node": node_counter, "name": indication_code, "url":'https://icd.who.int/browse/2024-01/mms/en#'+indication_uri,"column":"x4"})
            node_counter += 1
            caches['indication'].append(indication_code)
        indi_node = next((item['node'] for item in sankey['nodes'] if item['name'] == indication_code), None)

        if indication_0 not in caches['level_0']:
            sankey['nodes'].append({"node": node_counter, "name": indication_0, "url":'https://icd.who.int/browse/2024-01/mms/en#'+uri,"column":"x3"})
            node_counter += 1
            caches['level_0'].append(indication_0)
        level_0_node = next((item['node'] for item in sankey['nodes'] if item['name'] == indication_0), None)

        if [ligand_name, ligand_id] not in caches['ligands']:
            sankey['nodes'].append({"node": node_counter, "name": ligand_name, "url":'/ligand/'+str(ligand_id)+'/group',"column":"x2"})
            node_counter += 1
            caches['ligands'].append([ligand_name, ligand_id])
        lig_node = next((item['node'] for item in sankey['nodes'] if item['name'] == ligand_name), None)

        if protein_name not in caches['targets']:
            sankey['nodes'].append({"node": node_counter, "name": protein_name, "url":'/protein/'+str(target_name),"column":"x1"})
            node_counter += 1
            caches['targets'].append(protein_name)
            caches['entries'].append(target_name)
        prot_node = next((item['node'] for item in sankey['nodes'] if item['name'] == protein_name), None)

        # create matrix row
        row = {
            'row_id': row_id,
            'x1': prot_node,
            'x2': lig_node,
            'x3': level_0_node,
            'x4': indi_node,
            'x1_name': protein_name,
            'x2_name': ligand_name,
            'x3_name': indication_0,
            'x4_name': indication_name
        }

        sankey['links'].append({"source": prot_node, "target": lig_node, "value": 1, "ligtrace": protein_name, "prottrace": indication_name, "linkage_key": "primary","link_identifier": link_id})  # x1 -> x2 (Protein → Ligand)
        link_id += 1
        sankey['links'].append({"source": lig_node, "target": level_0_node, "value": 1, "ligtrace": protein_name, "prottrace": indication_0, "linkage_key": "primary","link_identifier": link_id})  # x2 -> x3 (Ligand → Level 0)
        link_id += 1
        sankey['links'].append({"source": level_0_node, "target": indi_node, "value": 1, "ligtrace": protein_name, "prottrace": indication_name, "linkage_key": "primary","link_identifier": link_id})  # x3 -> x4 (Level 0 → Indication)
        link_id += 1
        path_matrix.append(row)
        row_id += 1

    #Fixing redundancy in sankey['links']
    unique_combinations = {}

    for d in sankey['links']:
        # Create a key based on source and target for identifying unique combinations
        key = (d['source'], d['target'])

        if key in unique_combinations:
            # If the combination exists, add the value to the existing entry
            unique_combinations[key]['value'] += d['value']
        else:
            # If it's a new combination, add it to the dictionary
            unique_combinations[key] = d

    # Convert the unique_combinations back to a list of dictionaries
    sankey['links'] = list(unique_combinations.values())

    # Find all nodes that are actually used in the links
    used_nodes = set()
    for link in sankey['links']:
        used_nodes.add(link['source'])
        used_nodes.add(link['target'])

    # Build new nodes list and a mapping from old to new node indices
    old_to_new = {}
    new_nodes = []
    new_index = 0
    for node_dict in sankey['nodes']:
        old_idx = node_dict["node"]
        if old_idx in used_nodes:
            # Only keep nodes that are actually used
            new_node_dict = {
                "node": new_index,
                "name": node_dict["name"],
                "url": node_dict["url"],
                "column": node_dict['column']
            }
            new_nodes.append(new_node_dict)
            old_to_new[old_idx] = new_index
            new_index += 1

    # Update the source and target indices in the links
    for link in sankey['links']:
        link["source"] = old_to_new[link["source"]]
        link["target"] = old_to_new[link["target"]]

    # Replace the nodes list with the new one
    sankey['nodes'] = new_nodes

    total_points = len(caches['targets']) + len(caches['targets']) + 1
    if len(caches['ligands']) > len(caches['targets']):
        nodes_nr = len(caches['ligands'])
    else:
        nodes_nr = len(caches['targets'])

    sankey_data = {
        'sankey': sankey,
        'total_points': total_points,
        'nodes_nr': nodes_nr,
        'path_matrix': path_matrix
    }

    return sankey_data