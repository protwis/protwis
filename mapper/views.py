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
import random
import pandas as pd
import random
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
import umap
import numpy as np



class LandingPage(TemplateView):
    template_name = 'mapper/data_mapper_landing.html'

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        # input_dict = LandingPage.parse_data_from_xls()
        method = 'umap' #will be defined by user input
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
        # Generate the nested dictionary (should represent the xls data)
        for i in range(1, 21):
            main_key = f'key{i}'
            nested_dict[main_key] = {}
            for j in range(1, 41):
                nested_key = f'nestedkey{j}'
                nested_dict[main_key][nested_key] = round(random.uniform(0, 10), 2)

        # Convert the nested dictionary to a DataFrame
        data = pd.DataFrame(nested_dict).T

        # Dimensionality reduction functions
        def apply_pca(data, n_components=2):
            pca = PCA(n_components=n_components)
            reduced_data = pca.fit_transform(data)
            return reduced_data

        def apply_tsne(data, n_components=2):
            tsne = TSNE(n_components=n_components)
            reduced_data = tsne.fit_transform(data)
            return reduced_data

        def apply_umap(data, n_components=2):
            umap_model = umap.UMAP(n_components=n_components)
            reduced_data = umap_model.fit_transform(data)
            return reduced_data

        # Main function to select and apply dimensionality reduction
        def reduce_dimensions(data, method='pca', n_components=2):
            if method == 'pca':
                return apply_pca(data, n_components)
            elif method == 'tsne':
                return apply_tsne(data, n_components)
            elif method == 'umap':
                return apply_umap(data, n_components)
            else:
                raise ValueError("Method not recognized. Choose 'pca', 'tsne', or 'umap'.")

        # Example usage
        reduced_data = reduce_dimensions(data, method)

        # Prepare the data for visualization
        visualization_data = pd.DataFrame(reduced_data, columns=['x', 'y'])
        visualization_data['label'] = data.index
        json_data = visualization_data.to_json(orient='records')

        return json_data


    @staticmethod
    def parse_data_from_xls():
        data = {}
        #CODE
        #CODE
        #CODE
        #CODE
        return data
