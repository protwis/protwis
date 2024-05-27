from django.shortcuts import get_object_or_404, render, redirect
from django.views import generic
from django.http import JsonResponse, HttpResponse
from django.db.models import Q, F, Func, Value, Prefetch
from django.core.cache import cache
from django.views.decorators.cache import cache_page
from django.urls import reverse
from django.shortcuts import render
from django.conf import settings
from django.views.decorators.cache import cache_page
from django.views.generic import TemplateView
from django.http import HttpResponseNotAllowed
from django import forms


from protein.models import Protein, ProteinConformation, ProteinAlias, ProteinFamily, Gene, ProteinSegment
from residue.models import Residue
from structure.models import Structure, StructureModel, StructureExtraProteins
from interaction.models import ResidueFragmentInteraction,StructureLigandInteraction
from common.selection import Selection
from common.views import AbsBrowseSelection
from ligand.models import Ligand, LigandID
from common.phylogenetic_tree import PhylogeneticTreeGenerator

import json
from copy import deepcopy
from collections import OrderedDict
import umap
import numpy as np
import pandas as pd
import random
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler

import openpyxl
import os


class LandingPage(TemplateView):
    template_name = 'mapper/data_mapper_landing.html'

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        # input_dict = LandingPage.parse_data_from_xls()
        method = 'tsne' #will be defined by user input
        list_plot = LandingPage.generate_list_plot()
        tree_plot, tree_options = LandingPage.generate_tree_plot()
        cluster_plot = LandingPage.generate_cluster(method)

        context['list_data'] = json.dumps(list_plot)
        context['tree_dict'] = json.dumps(tree_plot)
        context['tree_options'] = tree_options
        context['cluster_data'] = cluster_plot
        context['plot_type'] = method
        # context['HeatMapData'] = json.dumps(HeatMapData)
        return context

    @staticmethod
    def keep_by_names(data, names_to_keep):
        if isinstance(data, list):
            # Process each item in the list
            kept_items = [keep_by_names(item, names_to_keep) for item in data]
            # Return only non-None items
            return [item for item in kept_items if item is not None]
        elif isinstance(data, OrderedDict):
            if data.get('name') not in names_to_keep:
                if 'children' in data:
                    # Recursively process children
                    data['children'] = LandingPage.keep_by_names(data['children'], names_to_keep)
                    # Remove the 'children' key if it's empty after processing
                    if not data['children']:
                        return None
                else:
                    return None
            else:
                # If the name is in the keep list, process children if present
                if 'children' in data:
                    data['children'] = LandingPage.keep_by_names(data['children'], names_to_keep)
                    if not data['children']:
                        del data['children']
            return data
        return data

    @staticmethod
    def convert_keys(datatree, conversion):
        new_tree = {}
        for key, value in datatree.items():
            # Convert the key using the conversion dictionary
            new_key = conversion.get(key, key)  # Fallback to the original key if no conversion is found
            if isinstance(value, dict):
                # Recursively convert keys of nested dictionaries
                new_tree[new_key] = LandingPage.convert_keys(value, conversion)
            else:
                # If the value is a list (end of the branch), just assign it
                new_tree[new_key] = value
        return new_tree

    @staticmethod
    def generate_list_plot(): #ADD AN INPUT FILTER DICTIONARY
        # Generate the master dict of protein families
        families = ProteinFamily.objects.all()
        datatree = {}
        conversion = {}
        for item in families:
            if len(item.slug) == 3 and item.slug not in datatree.keys():
                datatree[item.slug] = {}
                conversion[item.slug] = item.name
            if len(item.slug) == 7 and item.slug not in datatree[item.slug[:3]].keys():
                datatree[item.slug[:3]][item.slug[:7]] = {}
                conversion[item.slug] = item.name
            if len(item.slug) == 11 and item.slug not in datatree[item.slug[:3]][item.slug[:7]].keys():
                datatree[item.slug[:3]][item.slug[:7]][item.slug[:11]] = []
                conversion[item.slug] = item.name
            if len(item.slug) == 15 and item.slug not in datatree[item.slug[:3]][item.slug[:7]][item.slug[:11]]:
                datatree[item.slug[:3]][item.slug[:7]][item.slug[:11]].append(item.name)

        datatree2 = LandingPage.convert_keys(datatree, conversion)
        datatree2.pop('Parent family', None)
        return datatree2

    @staticmethod
    def generate_tree_plot(): #ADD AN INPUT FILTER DICTIONARY
        ### TREE SECTION
        tree = PhylogeneticTreeGenerator()
        class_a_data = tree.get_tree_data(ProteinFamily.objects.get(name='Class A (Rhodopsin)'))
        class_a_options = deepcopy(tree.d3_options)
        class_b1_data = tree.get_tree_data(ProteinFamily.objects.get(name__startswith='Class B1 (Secretin)'))
        class_b1_options = deepcopy(tree.d3_options)
        class_b2_data = tree.get_tree_data(ProteinFamily.objects.get(name__startswith='Class B2 (Adhesion)'))
        class_b2_options = deepcopy(tree.d3_options)
        class_c_data = tree.get_tree_data(ProteinFamily.objects.get(name__startswith='Class C (Glutamate)'))
        class_c_options = deepcopy(tree.d3_options)
        class_f_data = tree.get_tree_data(ProteinFamily.objects.get(name__startswith='Class F (Frizzled)'))
        class_f_options = deepcopy(tree.d3_options)
        class_t2_data = tree.get_tree_data(ProteinFamily.objects.get(name__startswith='Class T (Taste 2)'))
        class_t2_options = deepcopy(tree.d3_options)
        ### GETTING NODES
        data_a = class_a_data.get_nodes_dict(None)
        data_b1 = class_b1_data.get_nodes_dict(None)
        data_b2 = class_b2_data.get_nodes_dict(None)
        data_c = class_c_data.get_nodes_dict(None)
        data_f = class_f_data.get_nodes_dict(None)
        data_t2 = class_t2_data.get_nodes_dict(None)
        #Collating everything into a single tree
        general_options = {'depth': 4,
                           'branch_length': {1: 'Class A (Rhodopsin)',
                                             2: 'Alicarboxylic acid',
                                             3: 'Gonadotrophin-releasing hormone',
                                             4: ''},
                           'branch_trunc': 0,
                           'leaf_offset': 30,
                           'anchor': "tree_plot",
                           'label_free': []}
        master_dict = OrderedDict([('name', ''),
                                   ('value', 3000),
                                   ('color', ''),
                                   ('children',[])])
        class_a_dict = OrderedDict([('name', 'Class A (Rhodopsin)'),
                                   ('value', 0),
                                   ('color', 'Red'),
                                   ('children',data_a['children'])])
        class_b1_dict = OrderedDict([('name', 'Class B1 (Secretin)'),
                                   ('value', 0),
                                   ('color', 'Green'),
                                   ('children',data_b1['children'])])
        class_b2_dict = OrderedDict([('name', 'Class B2 (Adhesion)'),
                                  ('value', 0),
                                  ('color', 'Blue'),
                                  ('children',data_b2['children'])])
        class_c_dict = OrderedDict([('name', 'Class C (Glutamate)'),
                                  ('value', 0),
                                  ('color', 'Purple'),
                                  ('children',data_c['children'])])
        class_f_dict = OrderedDict([('name', 'Class F (Frizzled)'),
                                  ('value', 0),
                                  ('color', 'Grey'),
                                  ('children',data_f['children'])])
        class_t2_dict = OrderedDict([('name', 'Class T (Taste 2)'),
                                  ('value', 0),
                                  ('color', 'Orange'),
                                  ('children',data_t2['children'])])
        ### APPENDING TO MASTER DICT
        master_dict['children'].append(class_a_dict)
        master_dict['children'].append(class_b1_dict)
        master_dict['children'].append(class_b2_dict)
        master_dict['children'].append(class_c_dict)
        master_dict['children'].append(class_f_dict)
        master_dict['children'].append(class_t2_dict)

        # master_dict = LandingPage.keep_by_names(filter_dict)

        return master_dict, general_options

    @staticmethod
    def generate_circles_data():
        gpcrs = Protein.objects.filter(species_id=1).values_list('entry_name', flat=True)
        circles = {}
        for gpcr in gpcrs:
            circles[gpcr] = {'Inner': 0,
                             'Outer1': 0,
                             'Outer2': 0,
                             'Outer3': 0,
                             'Outer4': 0,
                             'Outer5': 0}
        return circles

        # @staticmethod
        # def generate_heatmap():

    @staticmethod
    def generate_cluster(method):
        # Initialize the test dictionary (should represent the xls data)
        nested_dict = {}
        proteins = list(Protein.objects.filter(species_id=1, accession__isnull=False).values_list('entry_name', flat=True).distinct())
        # Generate the nested dictionary (should represent the xls data)
        for i in proteins:
            main_key = i.split('_human')[0]
            nested_dict[main_key] = {}
            for j in range(1, 81):
                nested_key = f'Variable{j}'
                nested_dict[main_key][nested_key] = round(random.uniform(0, 100), 2)

        # Convert the nested dictionary to a DataFrame
        data = pd.DataFrame(nested_dict).T

        def reduce_and_cluster(data, method='umap', n_components=2, n_clusters=10):
            if method == 'umap':
                reducer = umap.UMAP(n_components=n_components, random_state=42)
            elif method == 'tsne':
                reducer = TSNE(n_components=n_components, random_state=42)
            elif method == 'pca':
                reducer = PCA(n_components=n_components, random_state=42)
            else:
                raise ValueError("Method should be either 'umap' or 'tsne'")

            reduced_data = reducer.fit_transform(data)

            # Clustering the reduced data
            kmeans = KMeans(n_clusters=n_clusters, random_state=42)
            clusters = kmeans.fit_predict(reduced_data)

            # Prepare the data for D3.js
            df = pd.DataFrame(reduced_data, columns=['x', 'y'])
            df['cluster'] = clusters
            df['label'] = data.index

            return df

        # Example usage
        reduced_df = reduce_and_cluster(data, method=method)

        # Prepare the data for visualization
        data_json = reduced_df.to_json(orient='records')

        return data_json


    @staticmethod
    def parse_data_from_xls():
        data = {}
        #CODE
        #CODE
        #CODE
        #CODE
        return data

    def post(self, request, *args, **kwargs):

        #############################################################
        ### This method handles POST requests for form submission ###
        #############################################################

        if request.method == 'POST':

            #################################
            # Utilize ExcelUploadForm class #
            #################################

            form = ExcelUploadForm(request.POST,request.FILES)

            ####################
            # If form is valid #
            ####################

            if form.is_valid():

                ####################
                # Get cleaned data #
                ####################

                file = form.cleaned_data['file']

                ##########################
                # Check if file is .xlsx #
                ##########################

                if not file.name.endswith('.xlsx'):
                    return render(request, self.template_name, {'upload_status': 'Failed'})
                else:

                    ##################################################
                    # Load excel file (workbook) and get sheet names #
                    ##################################################

                    try:
                        workbook = openpyxl.load_workbook(filename=file,read_only=False)

                        sheet_names = workbook.sheetnames

                        ##################################
                        # For each sheet in the workbook #
                        ##################################

                        for sheet_name in sheet_names:

                            ########################
                            # Initialize worksheet #
                            ########################

                            worksheet = workbook.get_sheet_by_name(sheet_name)

                            ###############################
                            # check if headers is correct #
                            ###############################

                            header_list = [worksheet['A1'].value, worksheet['B1'].value, worksheet['C1'].value]

                            ###################################################
                            # If first sheet is receptor with correct headers #
                            ###################################################

                            if sheet_name == 'Phylogenetic Tree' and header_list == ['Receptor (Uniprot)', '1. Feature (Inner cicle)', '2. Order (Outer cicle 1)']:

                                print("success!")

                                #########################################
                                ## if everything is good in do a query ##
                                ##  Fetch: Protein, family and class   ##
                                #########################################

                                protein_data = Protein.objects.filter(entry_name__endswith='_human').prefetch_related('family__parent__parent', 'family__parent__parent__parent').distinct('entry_name')

                                ######################################
                                # Initialize sets for unique entries #
                                ######################################

                                unique_entry_names = set()
                                unique_protein_families = set()
                                unique_protein_classes = set()

                                ##################################
                                # for each entry add to sets of  #
                                # proteins, families and classes #
                                ##################################

                                for entry in protein_data:

                                    ########################
                                    # Initiate the entries #
                                    ########################

                                    protein = str(entry)
                                    protein_family = str(entry.family.parent.parent.name)
                                    protein_class = str(entry.family.parent.parent.parent.name)

                                    ##########################
                                    # Populate the set-lists #
                                    ##########################

                                    unique_entry_names.add(protein)
                                    unique_protein_families.add(protein_family)
                                    unique_protein_classes.add(protein_class)

                                #########################
                                # convert sets to lists #
                                #########################

                                list_unique_entry_names = list(unique_entry_names)

                                ###################################
                                ## Retrieve all cell values from ##
                                ## columns A, B, and C as lists  ##
                                ###################################

                                ############################
                                # Initiate lists and dicts #
                                ############################

                                headers = [cell.value for cell in worksheet[1]]
                                data_types = [cell.value for cell in worksheet[2]]

                                Correct_values = {}
                                Incorrect_values = {}
                                missing_data_columns = {}

                                ##########################
                                # Run through excel file #
                                ##########################

                                # Identify completely empty columns and skip them
                                empty_columns = set()
                                for col_idx in range(len(headers)):
                                    if all(cell is None for cell in worksheet.iter_cols(min_col=col_idx + 1, max_col=col_idx + 1, min_row=3, values_only=True)):
                                        empty_columns.add(col_idx)
                                print(empty_columns)
                                # Iterate through rows starting from the third row
                                for index, row in enumerate(worksheet.iter_rows(min_row=3, values_only=True), start=3):
                                    # Check the "Receptor (Uniprot)" column for correct values
                                    if row[0] in list_unique_entry_names:
                                        Correct_values[index] = 'Entry correct'
                                    else:
                                        Incorrect_values[index] = 'Wrong entry'

                                    # Check each column for data points, boolean values, and float values
                                    for col_idx, value in enumerate(row):
                                        if col_idx == 0 or col_idx in empty_columns:
                                            continue  # Skip the "Receptor (Uniprot)" column and completely empty columns

                                        if value is not None:
                                            if col_idx not in missing_data_columns:
                                                missing_data_columns[col_idx] = True  # Column has data points
                                            if data_types[col_idx] == 'Boolean' and str(value).lower() not in ['yes', 'no', '1', '0']:
                                                if col_idx not in Incorrect_values:
                                                    Incorrect_values[col_idx] = 'Non-boolean value found'
                                            elif data_types[col_idx] == 'Float':
                                                try:
                                                    float(value)
                                                except ValueError:
                                                    if col_idx not in Incorrect_values:
                                                        Incorrect_values[col_idx] = 'Non-float value found'
                                        else:
                                            if col_idx not in missing_data_columns:
                                                missing_data_columns[col_idx] = False  # Column has missing data points
                                            elif missing_data_columns[col_idx] is True:
                                                missing_data_columns[col_idx] = False  # Column has mixed data points

                                # Prepare the final output
                                for col_idx, has_data in missing_data_columns.items():
                                    if not has_data:
                                        if col_idx not in Incorrect_values:
                                            Incorrect_values[col_idx] = 'Missing data points found'

                                print("Correct Values: ", Correct_values)
                                print("Incorrect Values: ", Incorrect_values)

                                if Incorrect_values:
                                    print("Incorrect values found")
                                    return render(request, self.template_name,{'upload_status': 'Success', 'report_status': 'Failed', 'incorrect_values': Incorrect_values})
                                else:
                                    print("All values are correct")
                                    return render(request, self.template_name,{'upload_status': 'Success', 'report_status': 'Success'})
                    except:
                        return render(request, self.template_name, {'upload_status': 'Failed'})
            else:
                # Return a 405 Method Not Allowed response if not a POST request
                #return HttpResponseNotAllowed(['POST'])
                return render(request, self.template_name, {'upload_status': 'Failed'})

# def LandingPage(request):
#     return render(request, 'mapper/data_mapper_landing.html')

#######################
## Excel upload form ##
#######################

class ExcelUploadForm(forms.Form):
    file = forms.FileField()
