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
from django import template


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
import umap.umap_ as umap
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
        ################################
        # Will be implemented later !! #
        ################################
        # method = 'umap' #will be defined by user input
        # list_plot = LandingPage.generate_list_plot()
        # tree_plot, tree_options = LandingPage.generate_tree_plot()
        # cluster_plot = LandingPage.generate_cluster(method)

        # context['list_data'] = json.dumps(list_plot)
        # context['tree_dict'] = json.dumps(tree_plot)
        # context['tree_options'] = tree_options
        # context['cluster_data'] = cluster_plot
        # context['plot_type'] = method
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
    def generate_cluster(method,input):
        

        # # Initialize the test dictionary (should represent the xls data)
        # nested_dict = {}
        # # Generate the nested dictionary (should represent the xls data)
        # for i in range(1, 21):
        #     main_key = f'GPCR{i}'
        #     nested_dict[main_key] = {}
        #     for j in range(1, 41):
        #         nested_key = f'Variable{j}'
        #         nested_dict[main_key][nested_key] = round(random.uniform(0, 100), 2)

        # Convert the nested dictionary to a DataFrame
        data = pd.DataFrame(input).T

        def reduce_and_cluster(data, method='umap', n_components=2, n_clusters=5):
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
                    return render(request, self.template_name, {'upload_status': 'Failed','Error_message': "The uploaded file is not an .xlsx file."})
                else:
                    try:
                        workbook = openpyxl.load_workbook(filename=file,read_only=False)
                    except:
                        return render(request, self.template_name, {'upload_status': 'Failed','Error_message': "Unable to load excel file, might be corrupted or not inline with the template file."})

                    if workbook:
                        ######################
                        # Fetch protein data #
                        ######################

                        protein_data = list(Protein.objects.filter(species=1).values_list('entry_name', flat=True).distinct())

                        ##################################################
                        # Load excel file (workbook) and get sheet names #
                        ##################################################

                        sheet_names = workbook.sheetnames
                        
                        ######################
                        # Sheets and headers #
                        ######################

                        # This should probably be excluded depending on what we want from the setup #

                        Phylogenetic_Tree_Header = ['Receptor (Uniprot)', '1. Feature (Inner cicle)', '2. Order (Outer cicle 1)', '3. Order (Outer cicle 2)', '4. Order (Outer cicle 3)', '5. Order (Outer cicle 4)', '6. Order (Outer cicle 5)']
                        Cluster_Analysis_Header = ['Receptor (Uniprot)','Feature 1','Feature 2','Feature 3','Feature 4']
                        List_Plot_Header = ['Receptor (Uniprot)','Feature 1','Feature 2']
                        Heatmap_Header = ['Receptor (Uniprot)','Feature 1','Feature 2','Feature 3','Feature 4','Feature 5']
                        Sheet_Header_pass_check = [False,False,False,False,False]
                        
                        # Check all sheet names, headers and subheaders (needs to be implemented) #

                        for sheet_name in sheet_names:
                            worksheet = workbook[sheet_name]
                            header_list = [cell.value for cell in worksheet[1]]
                            if sheet_name == 'Info':
                                Sheet_Header_pass_check[0] = True
                            elif sheet_name == 'Phylogenetic Tree' and header_list == Phylogenetic_Tree_Header:
                                Sheet_Header_pass_check[1] = True
                            elif sheet_name == 'Cluster Analysis' and header_list[:5] == Cluster_Analysis_Header:
                                Sheet_Header_pass_check[2] = True
                            elif sheet_name == 'List Plot' and header_list == List_Plot_Header:
                                Sheet_Header_pass_check[3] = True
                            elif sheet_name == 'Heatmap' and header_list == Heatmap_Header:
                                Sheet_Header_pass_check[4] = True
                            else:
                                pass
                        
                        if not all(Sheet_Header_pass_check):
                            # Add addition for the different sheets.
                            return render(request, self.template_name, {'upload_status': 'Failed','Error_message': "The excel file is not structured as the template file. There are incorrect headers and subheaders."})
                        else:
                            
                            # Init incorrect values #
                            
                            plot_names = ['Phylogenetic Tree', 'Cluster Analysis', 'List Plot', 'Heatmap']
                            Data = {}
                            Incorrect_values = {}

                            for key in plot_names:
                                Data[key] = {}
                                Incorrect_values[key] = {}
                            
                            Plot_parser = ['Failed','Failed','Failed','Failed']

                            ##################################
                            # For each sheet in the workbook #
                            ##################################

                            for sheet_name in sheet_names:

                                ########################
                                # Initialize worksheet #
                                ########################

                                worksheet = workbook[sheet_name]
                                
                                try:
                                    header_list = [cell.value for cell in worksheet[1]]
                                except:
                                    return render(request, self.template_name, {'upload_status': 'Failed','Error_message': "Corrupted excel, headers not inline with the template file."})
                                
                                
                                ###################################################
                                # If first sheet is receptor with correct headers #
                                ###################################################

                                if sheet_name == 'Phylogenetic Tree':
                                    
                                    # Initialize dictionaries
                                    data_types = [cell.value for cell in worksheet[2]]
                                    for key in header_list:
                                        Incorrect_values[sheet_name][key] = {}

                                    #######################################
                                    # Run through Phylogenetic tree sheet #
                                    #######################################
                                    try:

                                        empty_sheet = True  # Initialize the flag

                                        # Iterate over rows starting from the second row (excluding the header row)
                                        for row in worksheet.iter_rows(min_row=3, values_only=True):
                                            # Check only the columns that have headers, skipping the first column
                                            if any(row[i] is not None for i, header in enumerate(header_list[1:], start=1) if header):
                                                empty_sheet = False
                                                break

                                        if empty_sheet:
                                            pass
                                        else:
                                            # Iterate through rows starting from the third row
                                            for index, row in enumerate(worksheet.iter_rows(min_row=3, values_only=True), start=3):
                                                # Check the "Receptor (Uniprot)" column for correct values
                                                if row[0] not in protein_data:
                                                    Incorrect_values[sheet_name][header_list[0]][index] = '"{}" is a invalid entry'.format(row[0])
                                                else:
                                                    if row[0] not in Data[sheet_name]:
                                                        Data[sheet_name][row[0]] = {}
                                                    
                                                    # Check each column for data points, boolean values, and float values
                                                    for col_idx, value in enumerate(row):
                                                        if col_idx == 0:
                                                            continue  # Skip the "Receptor (Uniprot)" column and completely empty columns
                                                        elif data_types[col_idx] not in ['Boolean','Number','Text']:
                                                            Incorrect_values[sheet_name][header_list[col_idx]] = 'Incorrect datatype'
                                                        else:
                                                            if value is not None:
                                                                if data_types[col_idx] == 'Boolean':
                                                                    if str(value).lower() not in ['yes', 'no', '1', '0']:
                                                                        Incorrect_values[sheet_name][header_list[col_idx]][index] = 'Non-Boolean Value'
                                                                    else:
                                                                        Data[sheet_name][row[0]]['Value{}'.format(col_idx)] = value
                                                                elif data_types[col_idx] == 'Number':
                                                                    try:
                                                                        float_value = float(value)
                                                                        Data[sheet_name][row[0]]['Value{}'.format(col_idx)] = float_value
                                                                    except ValueError:
                                                                        Incorrect_values[sheet_name][header_list[col_idx]][index] = 'Non-Number Value'
                                                                else:
                                                                    pass
                                                            else:
                                                                pass

                                        # Check if any values are incorrect #
                                        status = 'Success'

                                        if empty_sheet:
                                            status = 'Empty sheet'
                                        elif Data[sheet_name]:
                                            # print(Data[sheet_name])
                                            for col_idx in Incorrect_values[sheet_name]:
                                                # Check if there are any assigned index values for this col_idx
                                                if any(Incorrect_values[sheet_name][col_idx].values()):
                                                    # If any index is assigned, set status to 'Partially_success' and break out of the loop
                                                    status = 'Failed'
                                                    break
                                        else:
                                            status = 'Failed'
                                        
                                        ## Update Plot parser ##
                                        Plot_parser[0] = status
                                    except:
                                        print("phylo failed")
                                
                                ### Cluster Analysis ###
                                elif sheet_name == 'Cluster Analysis':
                                    
                                    # Initialize dictionaries
                                    data_types = [cell.value for cell in worksheet[2]]
                                    for key in header_list:
                                        Incorrect_values[sheet_name][key] = {}
                                    try:
                                        
                                        empty_sheet = True  # Initialize the flag

                                        # Iterate over rows starting from the second row (excluding the header row)
                                        for row in worksheet.iter_rows(min_row=3, values_only=True):
                                            # Check only the columns that have headers, skipping the first column
                                            if any(row[i] is not None for i, header in enumerate(header_list[1:], start=1) if header):
                                                empty_sheet = False
                                                break

                                        if empty_sheet:
                                            pass
                                        else:
                                        
                                            # Iterate through rows starting from the second row
                                            for index, row in enumerate(worksheet.iter_rows(min_row=3, values_only=True), start=3):
                                                # Check the "Receptor (Uniprot)" column for correct values
                                                if row[0] not in protein_data:
                                                    Incorrect_values[sheet_name][header_list[0]][index] = '"{}" is a invalid entry'.format(row[0])
                                                else:
                                                    if row[0] not in Data[sheet_name]:
                                                        Data[sheet_name][row[0]] = {}
                                                        
                                                    # Check each column for data points, boolean values, and float values #
                                                    for col_idx, value in enumerate(row):
                                                        if col_idx == 0:
                                                            continue  # Skip the "Receptor (Uniprot)" column and completely empty columns #
                                                        elif data_types[col_idx] not in ['Boolean','Number','Text']:
                                                            Incorrect_values[sheet_name][header_list[col_idx]] = 'Incorrect datatype'
                                                        else:
                                                            if value is not None:
                                                                # Handle the 3 different types of input for Cluster analysis (Boolean, Number, and Text) #
                                                                if data_types[col_idx] == 'Boolean':
                                                                    if str(value).lower() not in ['yes', 'no', '1', '0']:
                                                                        Incorrect_values[sheet_name][header_list[col_idx]][index] = 'Non-Boolean Value'
                                                                    else:
                                                                        Data[sheet_name][row[0]]['Value{}'.format(col_idx)] = value
                                                                elif data_types[col_idx] == 'Number':
                                                                    try:
                                                                        float_value = float(value)
                                                                        Data[sheet_name][row[0]]['Value{}'.format(col_idx)] = float_value
                                                                    except ValueError:
                                                                        Incorrect_values[sheet_name][header_list[col_idx]][index] = 'Non-Number Value'
                                                                elif data_types[col_idx] == 'Text':
                                                                    if isinstance(value, str):
                                                                        Data[sheet_name][row[0]]['Value{}'.format(col_idx)] = value
                                                                    else:
                                                                        Incorrect_values[sheet_name][header_list[col_idx]][index] = 'Non-Text Value'
                                                                else:
                                                                    pass
                                                            else:
                                                                pass
                                            
                                        # Check if any values are incorrect #
                                        status = 'Success'

                                        if empty_sheet:
                                            status = 'Empty sheet'
                                        elif Data[sheet_name]:
                                            for col_idx in Incorrect_values[sheet_name]:
                                                # Check if there are any assigned index values for this col_idx
                                                if any(Incorrect_values[sheet_name][col_idx].values()):
                                                    # If any index is assigned, set status to 'Partially_success' and break out of the loop
                                                    status = 'Failed'
                                                    break
                                        else:
                                            status = 'Failed'
                                        
                                        ## Update Plot_parser for Cluster Analysis
                                        Plot_parser[1] = status
                                    except:
                                        print("Cluster failed")
                                
                                ### List Plot ###
                                elif sheet_name == 'List Plot':
                                    
                                    # Initialize dictionaries
                                    data_types = [cell.value for cell in worksheet[2]]
                                    for key in header_list:
                                        Incorrect_values[sheet_name][key] = {}
                                    try:
                                        empty_sheet = True  # Initialize the flag

                                        # Iterate over rows starting from the second row (excluding the header row)
                                        for row in worksheet.iter_rows(min_row=3, values_only=True):
                                            # Check only the columns that have headers, skipping the first column
                                            if any(row[i] is not None for i, header in enumerate(header_list[1:], start=1) if header):
                                                empty_sheet = False
                                                break

                                        if empty_sheet:
                                            pass
                                        else:
                                            # Iterate through rows starting from the second row
                                            for index, row in enumerate(worksheet.iter_rows(min_row=3, values_only=True), start=3):
                                                # Check the "Receptor (Uniprot)" column for correct values
                                                if row[0] not in protein_data:
                                                    Incorrect_values[sheet_name][header_list[0]][index] = '"{}" is a invalid entry'.format(row[0])
                                                else:
                                                    if row[0] not in Data[sheet_name]:
                                                        Data[sheet_name][row[0]] = {}
                                                        
                                                    # Check each column for data points, boolean values, and float values #
                                                    for col_idx, value in enumerate(row):
                                                        if col_idx == 0:
                                                            continue  # Skip the "Receptor (Uniprot)" column and completely empty columns #
                                                        elif data_types[col_idx] not in ['Boolean','Number']:
                                                            Incorrect_values[sheet_name][header_list[col_idx]] = 'Incorrect datatype'
                                                        else:
                                                            if value is not None:
                                                                # Handle the 2 different types of input for Cluster analysis (Boolean or Number) #
                                                                if data_types[col_idx] == 'Boolean':
                                                                    if str(value).lower() not in ['yes', 'no', '1', '0']:
                                                                        Incorrect_values[sheet_name][header_list[col_idx]][index] = 'Non-Boolean Value'
                                                                    else:
                                                                        Data[sheet_name][row[0]]['Value{}'.format(col_idx)] = value
                                                                elif data_types[col_idx] == 'Number':
                                                                    try:
                                                                        float_value = float(value)
                                                                        Data[sheet_name][row[0]]['Value{}'.format(col_idx)] = float_value
                                                                    except ValueError:
                                                                        Incorrect_values[sheet_name][header_list[col_idx]][index] = 'Non-Number Value'
                                                                else:
                                                                    pass
                                                            else:
                                                                pass
                                            
                                        # Check if any values are incorrect #
                                        status = 'Success'

                                        if empty_sheet:
                                            status = 'Empty sheet'
                                        elif Data[sheet_name]:
                                            for col_idx in Incorrect_values[sheet_name]:
                                                # Check if there are any assigned index values for this col_idx
                                                if any(Incorrect_values[sheet_name][col_idx].values()):
                                                    # If any index is assigned, set status to 'Partially_success' and break out of the loop
                                                    status = 'Failed'
                                                    break
                                        else:
                                            status = 'Failed'
                                    
                                        ## Update Plot_parser for Cluster Analysis
                                        Plot_parser[2] = status
                                    except:
                                        print("List plot failed")
                                
                                ### Heatmap ###
                                elif sheet_name == 'Heatmap':

                                    # Initialize dictionaries
                                    data_types = [cell.value for cell in worksheet[2]]
                                    for key in header_list:
                                        Incorrect_values[sheet_name][key] = {}
                                    try:

                                        empty_sheet = True  # Initialize the flag

                                        # Iterate over rows starting from the second row (excluding the header row)
                                        for row in worksheet.iter_rows(min_row=3, values_only=True):
                                            # Check only the columns that have headers, skipping the first column
                                            if any(row[i] is not None for i, header in enumerate(header_list[1:], start=1) if header):
                                                empty_sheet = False
                                                break

                                        if empty_sheet:
                                            pass
                                        else:

                                            # Iterate through rows starting from the second row
                                            for index, row in enumerate(worksheet.iter_rows(min_row=3, values_only=True), start=3):
                                                # Check the "Receptor (Uniprot)" column for correct values
                                                if row[0] not in protein_data:
                                                    Incorrect_values[sheet_name][header_list[0]][index] = '"{}" is a invalid entry'.format(row[0])
                                                else:
                                                    if row[0] not in Data[sheet_name]:
                                                        Data[sheet_name][row[0]] = {}
                                                        
                                                    # Check each column for data points, boolean values, and float values #
                                                    for col_idx, value in enumerate(row):
                                                        if col_idx == 0:
                                                            continue  # Skip the "Receptor (Uniprot)" column and completely empty columns #
                                                        elif data_types[col_idx] not in ['Number']:
                                                            Incorrect_values[sheet_name][header_list[col_idx]] = 'Incorrect datatype'
                                                        else:
                                                            if value is not None:
                                                                # Handle the 1 different types of input for Heatmap (Number) #
                                                                if data_types[col_idx] == 'Number':
                                                                    try:
                                                                        float_value = float(value)
                                                                        Data[sheet_name][row[0]]['Value{}'.format(col_idx)] = float_value
                                                                    except ValueError:
                                                                        Incorrect_values[sheet_name][header_list[col_idx]][index] = 'Non-Number Value'
                                                                else:
                                                                    pass
                                                            else:
                                                                pass
                                            
                                        # Check if any values are incorrect #
                                        status = 'Success'

                                        if empty_sheet:
                                            status = 'Empty sheet'
                                        elif Data[sheet_name]:
                                            for col_idx in Incorrect_values[sheet_name]:
                                                # Check if there are any assigned index values for this col_idx
                                                if any(Incorrect_values[sheet_name][col_idx].values()):
                                                    # If any index is assigned, set status to 'Partially_success' and break out of the loop
                                                    status = 'Failed'
                                                    break
                                        else:
                                            status = 'Failed'
                                        
                                        ## Update Plot_parser for Cluster Analysis
                                        Plot_parser[3] = status
                                    except:
                                        print("Heatmap Failed")
                            ## Return all values for plotparser and correctly (or partially) succesful plots ##
                            
                            plot_names = ['Phylogenetic Tree', 'Cluster Analysis', 'List Plot', 'Heatmap']
                            plot_data = {}
                            plot_incorrect_data = {}

                            for plot_name, plot_status in zip(plot_names, Plot_parser):
                                if plot_status == 'Success':
                                    plot_data[plot_name] = Data[plot_name]
                                elif plot_status == 'Failed':
                                    plot_incorrect_data[plot_name] = Incorrect_values[plot_name]
                            
                            plot_data_json = json.dumps(plot_data, indent=4, sort_keys=True) if plot_data else None
                            plot_incorrect_data_json = json.dumps(plot_incorrect_data, indent=4, sort_keys=True) if plot_incorrect_data else None
                            
                            Plot_parser_json = json.dumps([status == 'Success' for status in Plot_parser])
                            
                            plots_status = [{'status': status, 'plot_name': plot_name} for status, plot_name in zip(Plot_parser, plot_names)]
                            
                            context = {'upload_status': 'Success', 
                                       'report_status': 'Failed',
                                       'Plot_parser':Plot_parser,
                                       'Plot_parser_json':Plot_parser_json,
                                       'plot_names':plot_names,
                                       'plots_status':plots_status}
                            print(Plot_parser)
                            if plot_data:
                                if 'Success' in Plot_parser and not all(status == 'Success' for status in Plot_parser):
                                    print("partially")
                                    context['report_status'] = 'Partially_success'
                                    context['Data'] = plot_data_json
                                elif all(status == 'Success' for status in Plot_parser):
                                    print("success")
                                    context['report_status'] = 'Success'
                                    context['Data'] = plot_data_json
                            else:
                                context['Data'] = "No Data"

                            if plot_incorrect_data:
                                context['Incorrect_data_json'] = plot_incorrect_data
                            else:
                                print("No incorrect data")
                                context['Incorrect_data_json'] = "No incorrect data"
                            return render(request, self.template_name, context)

                                    
                    else:
                        return render(request, self.template_name, {'upload_status': 'Failed','Error_message': "Unable to load excel file, might be corrupted or not inline with the template file."})
                    # Needs to send a response if everything is handled #
                    #return render(request, self.template_name, {'upload_status': 'Success'})
                                
            else:
                # Return a 405 Method Not Allowed response if not a POST request
                #return HttpResponseNotAllowed(['POST'])
                return render(request, self.template_name, {'upload_status': 'Failed','Error_message': "Not a valid excel file. Please try and use the template excel file."})

# def LandingPage(request):
#     return render(request, 'mapper/data_mapper_landing.html')

#######################
## Excel upload form ##
#######################

class ExcelUploadForm(forms.Form):
    file = forms.FileField()

class plotrender(TemplateView):
    template_name = 'mapper/data_mapper_plotrender.html'

    def post(self, request, *args, **kwargs):
        # Retrieve the sample data from the POST request
        Plot_evaluation_json = request.POST.get('Plot_evaluation')
        Data_json = request.POST.get('Data')
        # If Plot_evaluation_json is not None, parse it as JSON
        if Plot_evaluation_json and Data_json:
            try:
                Plot_evaluation = json.loads(Plot_evaluation_json)
                Data = json.loads(Data_json)
            except json.JSONDecodeError:
                # Handle the case when the JSON data is invalid
                return HttpResponse("Invalid JSON data")

            # Add the sample data to the context
            #output = LandingPage.generate_cluster('umap',Data)
            context = {'Plot_evaluation_json': Plot_evaluation,'Data': Data}

            # Return the context dictionary
            return self.render_to_response(context)
        else:
            # Handle the case when Plot_evaluation_json is None
            # This could happen if the form was submitted without the JSON data
            return HttpResponse("Missing sample data")