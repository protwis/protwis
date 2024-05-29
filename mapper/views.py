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
        # cluster_plot = LandingPage.generate_test_cluster(method)
        # context['list_data'] = json.dumps(list_plot)
        # context['tree_dict'] = json.dumps(tree_plot)
        # context['tree_options'] = tree_options
        # context['cluster_data'] = cluster_plot
        # context['plot_type'] = method
        # context['HeatMapData'] = json.dumps(HeatMapData)
        return context

    @staticmethod
    def keep_by_names(data, names_to_keep):
        data_copy = deepcopy(data)
        if isinstance(data_copy, list):
            # Process each item in the list
            kept_items = [LandingPage.keep_by_names(item, names_to_keep) for item in data_copy]
            # Return only non-None items
            return [item for item in kept_items if item is not None]
        elif isinstance(data_copy, OrderedDict):
            name = data_copy.get('name')
            if name not in names_to_keep.keys():
                if 'children' in data_copy:
                    # Recursively process children
                    data_copy['children'] = LandingPage.keep_by_names(data_copy['children'], names_to_keep)
                    # Remove the 'children' key if it's empty after processing
                    if not data_copy['children']:
                        return None
                else:
                    return None
            else:
                # If the name is in the keep list, update the 'value' field
                data_copy['value'] = names_to_keep[name]['Inner']
                # Process children if present
                if 'children' in data_copy:
                    data_copy['children'] = LandingPage.keep_by_names(data_copy['children'], names_to_keep)
                    if not data_copy['children']:
                        del data_copy['children']
            return data_copy
        return data_copy

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
    def generate_tree_plot(input_data): #ADD AN INPUT FILTER DICTIONARY
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

        updated_data = {key.replace('_human', ''): value for key, value in input_data.items()}
        circles = {key.replace('_human', '').upper(): {k: v for k, v in value.items() if k != 'Inner'} for key, value in input_data.items()}
        master_dict = LandingPage.keep_by_names(master_dict, updated_data)
        whole_receptors = Protein.objects.prefetch_related("family", "family__parent__parent__parent")
        whole_rec_dict = {}
        for rec in whole_receptors:
            rec_uniprot = rec.entry_short()
            rec_iuphar = rec.family.name.replace("receptor", '').replace("<i>", "").replace("</i>", "").strip()
            if (rec_iuphar[0].isupper()) or (rec_iuphar[0].isdigit()):
                whole_rec_dict[rec_uniprot] = [rec_iuphar]
            else:
                whole_rec_dict[rec_uniprot] = [rec_iuphar.capitalize()]

        return master_dict, general_options, circles, whole_rec_dict

    @staticmethod
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

    @staticmethod
    def generate_cluster(method, input):
        # Convert the nested dictionary to a DataFrame
        input = {key.replace('_human', '').upper(): value for key, value in input.items()}
        data = pd.DataFrame(input).T
        # Example usage
        reduced_df = LandingPage.reduce_and_cluster(data, method=method)
        # Prepare the data for visualization
        data_json = reduced_df.to_json(orient='records')

        return data_json

    # @staticmethod
    # def generate_test_cluster(method):
    #     # Initialize the test dictionary (should represent the xls data)
    #     nested_dict = {}
    #     # Generate the nested dictionary (should represent the xls data)
    #     proteins = list(Protein.objects.filter(entry_name__endswith='_human').values_list('entry_name', flat=True).distinct())
    #     for i in proteins:
    #         main_key = i.split('_')[0].upper()
    #         nested_dict[main_key] = {}
    #         for j in range(1, 81):
    #             nested_key = f'Variable{j}'
    #             nested_dict[main_key][nested_key] = round(random.uniform(0, 100), 2)
    #     # Convert the nested dictionary to a DataFrame
    #     data = pd.DataFrame(nested_dict).T
    #     # Example usage
    #     reduced_df = LandingPage.reduce_and_cluster(data, method=method)
    #     # Prepare the data for visualization
    #     data_json = reduced_df.to_json(orient='records')
    #     return data_json

    @staticmethod
    def map_to_quartile(value, quartiles):
        if value <= quartiles[0.25]:
            return 10
        elif value <= quartiles[0.5]:
            return 20
        elif value <= quartiles[0.75]:
            return 30
        else:
            return 40

    def post(self, request, *args, **kwargs):
        ### This method handles POST requests for form submission ###

        if request.method == 'POST':

            # Utilize ExcelUploadForm class #
            form = ExcelUploadForm(request.POST,request.FILES)

            # If form is valid #
            if form.is_valid():

                # Get cleaned data #
                file = form.cleaned_data['file']

                # Check if file is .xlsx #
                if not file.name.endswith('.xlsx'):
                    return render(request, self.template_name, {'upload_status': 'Failed','Error_message': "The uploaded file is not an .xlsx file."})
                else:
                    try:
                        workbook = openpyxl.load_workbook(filename=file,read_only=False)
                    except:
                        return render(request, self.template_name, {'upload_status': 'Failed','Error_message': "Unable to load excel file, might be corrupted or not inline with the template file."})

                    if workbook:

                        # Fetch protein data #
                        protein_data = Protein.objects.filter(entry_name__endswith='_human').prefetch_related('family__parent__parent', 'family__parent__parent__parent').distinct('entry_name')

                        # Initialize sets for unique entries #
                        unique_entry_names = set()
                        unique_protein_families = set()
                        unique_protein_classes = set()

                        # for each entry add to sets of  #
                        # proteins, families and classes #
                        for entry in protein_data:

                            # Initiate the entries #
                            protein = str(entry)
                            protein_family = str(entry.family.parent.parent.name)
                            protein_class = str(entry.family.parent.parent.parent.name)

                            # Populate the set-lists #
                            unique_entry_names.add(protein)
                            unique_protein_families.add(protein_family)
                            unique_protein_classes.add(protein_class)

                        # convert sets to lists #
                        list_unique_entry_names = list(unique_entry_names)

                        # Load excel file (workbook) and get sheet names #
                        sheet_names = workbook.sheetnames

                        # Sheets and headers #
                        Phylogenetic_Tree_Header = ['Receptor (Uniprot)', 'Inner', 'Outer 1', 'Outer 2', 'Outer 3', 'Outer 4', 'Outer 5']
                        Cluster_Analysis_Header = ['Receptor (Uniprot)','Feature 1','Feature 2','Feature 3','Feature 4']
                        List_Plot_Header = ['Receptor (Uniprot)','Feature 1','Feature 2','Feature 3']
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
                            elif sheet_name == 'Cluster Analysis' and header_list == Cluster_Analysis_Header:
                                Sheet_Header_pass_check[2] = True
                            elif sheet_name == 'List Plot' and header_list == List_Plot_Header:
                                Sheet_Header_pass_check[3] = True
                            elif sheet_name == 'Heatmap' and header_list == Heatmap_Header:
                                Sheet_Header_pass_check[4] = True
                            else:
                                pass

                        if not all(Sheet_Header_pass_check):
                            return render(request, self.template_name, {'upload_status': 'Failed','Error_message': "The excel file is not structured as the template file. There are incorrect headers and subheaders."})
                        else:

                            # Init incorrect values #
                            Incorrect_values = {}
                            Incorrect_values['Phylogenetic Tree'] = {}
                            Incorrect_values['Cluster Analysis'] = {}
                            Incorrect_values['List Plot'] = {}
                            Incorrect_values['Heatmap'] = {}

                            Plot_parser = ['Failed','Failed','Failed','Failed']

                            # For each sheet in the workbook #
                            for sheet_name in sheet_names:

                                # Initialize worksheet #
                                worksheet = workbook[sheet_name]

                                try:
                                    header_list = [cell.value for cell in worksheet[1]]
                                except:
                                    return render(request, self.template_name, {'upload_status': 'Failed','Error_message': "Corrupted excel, headers not inline with the template file."})

                                # If first sheet is receptor with correct headers #
                                if sheet_name == 'Phylogenetic Tree':

                                    # Locate the column index for the "Inner" header
                                    header = next(worksheet.iter_rows(min_row=1, max_row=1, values_only=True))
                                    inner_col_idx = header.index('Inner') + 1  # openpyxl uses 1-based indexing
                                    # Check the first value under the "Inner" header
                                    first_value = worksheet.cell(row=2, column=inner_col_idx).value
                                    if first_value != "Boolean":
                                        # Extract all values from the "Inner" column, skipping the header
                                        inner_values = []
                                        for row in worksheet.iter_rows(min_row=3, min_col=inner_col_idx, max_col=inner_col_idx, values_only=True):
                                            if row[0] is not None:
                                                inner_values.append(row[0])

                                        inner_series = pd.Series(inner_values)
                                        # Calculate the quartiles
                                        quartiles = inner_series.quantile([0.25, 0.5, 0.75])
                                        # Apply the function to the series, passing quartiles as an argument
                                        mapped_values = inner_series.apply(lambda x: self.map_to_quartile(x, quartiles))

                                        # If you need to update the worksheet with these values
                                        for i, value in enumerate(mapped_values, start=3):  # start=3 to skip the header row
                                            worksheet.cell(row=i, column=inner_col_idx).value = value

                                    data_types = [cell.value for cell in worksheet[2]]
                                    Phylogenetic_Tree_data = {}
                                    for col_idx in range(len(header_list)):
                                        Incorrect_values['Phylogenetic Tree'][col_idx] = {}

                                    # Run through Phylogenetic tree sheet #
                                    # Iterate through rows starting from the third row
                                    for index, row in enumerate(worksheet.iter_rows(min_row=3, values_only=True), start=3):
                                        # Check the "Receptor (Uniprot)" column for correct values
                                        if row[0] not in list_unique_entry_names:
                                            Incorrect_values['Phylogenetic Tree'][0][index] = 'Wrong entry'
                                        else:
                                            if row[0] not in Phylogenetic_Tree_data:
                                                Phylogenetic_Tree_data[row[0]] = {}

                                            # Check each column for data points, boolean values, and float values
                                            for col_idx, value in enumerate(row):
                                                if col_idx == 0:
                                                    continue  # Skip the "Receptor (Uniprot)" column and completely empty columns
                                                if value is not None:
                                                    if data_types[col_idx] == 'Boolean':
                                                        if str(value).lower() not in ['yes', 'no', '1', '0', 'Yes', 'No', 'X']:
                                                            Incorrect_values['Phylogenetic Tree'][col_idx][index] = 'Non-boolean value'
                                                        else:
                                                            #Formatting it to 2000/0 (red/white)
                                                            if value in ['yes', '1', 'Yes', 'X']:
                                                                value = 1
                                                            else:
                                                                value = 0
                                                            if col_idx == 1:
                                                                Phylogenetic_Tree_data[row[0]]['Inner'] = value
                                                            else:
                                                                Phylogenetic_Tree_data[row[0]]['Outer{}'.format(col_idx-1)] = value
                                                    elif data_types[col_idx] == 'Float':
                                                        try:
                                                            float_value = float(value)
                                                            if col_idx == 1:
                                                                Phylogenetic_Tree_data[row[0]]['Inner'] = value
                                                            else:
                                                                Phylogenetic_Tree_data[row[0]]['Outer{}'.format(col_idx-1)] = float_value
                                                        except ValueError:
                                                            Incorrect_values['Phylogenetic Tree'][col_idx][index] = 'Non-float value'
                                                else:
                                                    pass
                                    # Check if any values are incorrect #
                                    status = 'success'

                                    if Phylogenetic_Tree_data:
                                        Data_Phylogenetic_Tree = json.dumps(Phylogenetic_Tree_data, indent=4, sort_keys=True)
                                        for col_idx in Incorrect_values['Phylogenetic Tree']:
                                            # Check if there are any assigned index values for this col_idx
                                            if any(Incorrect_values['Phylogenetic Tree'][col_idx].values()):
                                                # If any index is assigned, set status to 'Partially_success' and break out of the loop
                                                status = 'Partially_success'
                                                break
                                    else:
                                        status = 'Failed'

                                    if status == 'success':
                                        # Needs to send a response if everything is handled #
                                        print('success')
                                        Plot_parser[0] = status
                                        sample_data_json = json.dumps([True,True,False,False])
                                        return render(request, self.template_name, {'upload_status': 'Success','report_status':'Success','Data_Phylogenetic_Tree': Data_Phylogenetic_Tree,'sample_data_json':sample_data_json})
                                    elif status == 'Partially_success':
                                        Plot_parser[0] = status
                                        print('Partially_success')
                                        Data_Incorrect_Phylogenetic_Tree = json.dumps(Incorrect_values['Phylogenetic Tree'],indent=4, sort_keys=True)
                                        return render(request, self.template_name, {'upload_status': 'Success','report_status':'Partially_success','Data_Phylogenetic_Tree': Data_Phylogenetic_Tree,'Data_Incorrect_Phylogenetic_Tree':Data_Incorrect_Phylogenetic_Tree})
                                    elif status == 'Failed':
                                        Plot_parser[0] = status
                                        return render(request, self.template_name, {'upload_status': 'Success','report_status':'Failed'})

                                elif sheet_name == 'Cluster Analysis':
                                    print("Cluster!")
                    # Needs to send a response if everything is handled #
                    return render(request, self.template_name, {'upload_status': 'Success'})

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
        sample_data_json = request.POST.get('sample_data')
        Phylogenetic_data_json = request.POST.get('Phylogenetic_data')
        # If sample_data_json is not None, parse it as JSON
        if sample_data_json and Phylogenetic_data_json:
            try:
                sample_data = json.loads(sample_data_json)
                Phylogenetic_data = json.loads(Phylogenetic_data_json)
            except json.JSONDecodeError:
                # Handle the case when the JSON data is invalid
                return HttpResponse("Invalid JSON data")
            # Add the sample data to the context
            tree, tree_options, circles, receptors = LandingPage.generate_tree_plot(Phylogenetic_data)
            output = LandingPage.generate_cluster('umap', Phylogenetic_data)
            context = {'sample_data_json': sample_data,'Phylogenetic_data': output}
            context['tree'] = json.dumps(tree)
            context['tree_options'] = tree_options
            context['circles'] = json.dumps(circles)
            context['whole_dict'] = json.dumps(receptors)

            # Return the context dictionary
            return self.render_to_response(context)
        else:
            # Handle the case when sample_data_json is None
            # This could happen if the form was submitted without the JSON data
            return HttpResponse("Missing sample data")
