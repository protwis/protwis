"""
A set of functions for generating statistics trees.
Annotates crystalized targets and number of ligands/target available in ChEMBL.
"""
from django.db.models import Count

from interaction.models import ResidueFragmentInteraction, StructureLigandInteraction
from ligand.models import AssayExperiment, BiasedData
from protein.models import Protein, ProteinFamily
from structure.models import Structure

import json
from collections import  OrderedDict
from copy import deepcopy


class PhylogeneticTreeNode(object):

    def __init__(self, name='', color=''):

        self.name = name
        self.color = color
        self.children = OrderedDict()
        self.exp_data = {
            'crystals': 0,
            'mutations': 0,
            'ligands': 0,
            'ligand_bias': 0,
            'ligand_bias_bal': 0,
            'pathway_pref': 0,
            'subtype' : 0,
            'subtype_bal' : 0,
            }

    def get_value(self, param):
        """
        Function returning a parameter for coloring tree leaves.
        @param: a parameter based on which a color value will be set.

        TODO: Implement a scheme for mutations.
        """
        try:
            return self.exp_data[param]
        except KeyError:
            return 0

    def increment_value(self, param, value=1):
        """
        Function returning a parameter for coloring tree leaves.
        @param: a parameter based on which a color value will be set.

        TODO: Implement a scheme for mutations.
        """
        self.exp_data[param] += value

    def update_exp_data(self, data):

        for key, value in data.items():
            if self.exp_data[key] > value: continue
            self.exp_data[key] = value


    def get_nodes_dict(self, param):
        if param == None:
            return OrderedDict([
                ('name', self.name),
                ('value', 3000),
                ('color', self.color),
                ('children', [
                    y.get_nodes_dict('crystals') for x,y in self.children.items() if self.children != OrderedDict()
                    ]),
                    ])
        else:
            return OrderedDict([
                ('name', self.name),
                ('value', self.get_value(param)),
                ('color', self.color),
                ('children', [
                    y.get_nodes_dict(param) for x,y in self.children.items() if self.children != OrderedDict()
                    ]),
                    ])


class PhylogeneticTree(object):

    def __init__(self, root_lvl, depth, family):

        self.tree = PhylogeneticTreeNode()


    def add_data(self, path, data):

        parent_path = path.split('_')[1:-1]
        tmp = self.tree.children
        tmp_path = [path.split('_')[0]]
        while parent_path != []:
            tmp_path.append(parent_path.pop(0))
            try:
                tmp = tmp['_'.join(tmp_path)].children
            except KeyError:
                tmp['_'.join(tmp_path)] = data
        try:
            tmp[path].update_exp_data(data.exp_data)
        except KeyError:
            tmp[path] = data


    def get_data(self, path):

        parent_path = path.split('_')[1:-1]
        tmp = self.tree.children
        tmp_path = [path.split('_')[0]]
        while parent_path != []:
            tmp_path.append(parent_path.pop(0))
            try:
                tmp = tmp['_'.join(tmp_path)].children
            except KeyError:
                print("You're screwed")
        try:
            print(path)
            print(tmp[path].name)
            print(tmp[path].exp_data)
        except:
            pass


    def get_nodes(self, level):

        ref = self.tree.children
        tmp = OrderedDict()
        while level:
            for child in ref.items():
                for grandchild in child[1].children.items():
                    tmp.update({grandchild[0]: grandchild[1]})
            ref = tmp
            tmp = OrderedDict()
            level = level -1
        return ref

    def get_nodes_dict(self, param):

        return self.tree.get_nodes_dict(param)


class PhylogeneticTreeGenerator(object):

    #TODO: This should go to settings as it is GPCR-specific.
    #Dict keys are the Class - Protein family pairs. '' means 'any'.
    #CSS_COLORS = {
    #   ("Class F (Frizzled)", '') : 'SteelBlue',
    #   ('Class A (Rhodopsin)', 'Protein') : 'SteelBlue',
    #   ('Class A (Rhodopsin)', 'Alicarboxylic acid') : 'Red',
    #   ('Class B2 (Adhesion)', '') : 'SteelBlue',
    #   ('Class A (Rhodopsin)', 'Peptide') : 'SkyBlue',
    #   ('Class B1 (Secretin)', '') : 'SkyBlue',
    #   ('Class A (Rhodopsin)', 'Lipid') : 'LightGreen',
    #   ('', 'Orphan') : 'Orange',
    #   ('Class A (Rhodopsin)', 'Sensory') : 'DarkGray',
    #   ('Class C (Glutamate)', 'Sensory') : 'DarkGray',
    #   ('Taste 2', '') : 'DarkGray',
    #   ('Class A (Rhodopsin)', 'Nucleotide') : 'Purple',
    #   }

    CSS_COLORS = {
        ("", "Adhesion receptors") : 'Crimson',
        ("", "Alicarboxylic acid receptors") : 'Red',
        ("", "Aminergic receptors") : 'OrangeRed',
        ("", "Amino acid receptors") : 'Orange',
        ("", "Ion receptors") : 'GoldenRod',
        ("", "Lipid receptors") : 'Gold',
        ("", "Melatonin receptors") : 'Yellow',
        ("", "Nucleotide receptors") : 'YellowGreen',
        ("", "Orphan receptors") : 'Gold',
        ("", "Other") : 'Green',
        ("", "Peptide receptors") : 'SkyBlue',
        ("", "Protein receptors") : 'SteelBlue',
        ("", "Sensory receptors") : 'Indigo',
        ("", "Steroid receptors") : 'Purple',
        ("Class B2 (Adhesion)", "") : 'LimeGreen',
        }
    #List of tree levels that should be sorted alphabetically
    SORTED_BRANCHES = [2,3]

#   o	Dark blue: class A Protein ligand type and whole of classes Adhesion and class F (they also have this ligand type)
#   o	Light blue: Peptide and whole of class B1 (it also has this ligand type)
#   o	Green: Lipid receptors
#   o	Orange: Orphan receptors
#   o	Dark grey: Sensory (class A opsins, class C Taste1 and whole of class Taste2)
#   o	Purple: Nucleotide (class A P2Y and adenosine)
#   o	Black: All other


#   http://www.d3noob.org/2014/01/tree-diagrams-in-d3js_11.html

    def __init__(self, root_lvl=1, depth=3):

        self.root_lvl = root_lvl
        self.tree_depth = depth

        self.lookup = { x: {} for x in range(self.tree_depth+1)}
        self.aux_data = {
                        'crystals': [],
                        'mutations': [],
                        'ligands': {},
                        'ligand_bias': {},
                        'ligand_bias_bal': {},
                        'pathway_pref': {},
                        'subtype': {},
                        'subtype_bal': {}
                        }

        self.get_aux_data()
        self.d3_options = {
            'depth': self.tree_depth,
            'branch_length': {},
            'branch_trunc': 0,
            'leaf_offset': 30
            }

        self.families = ProteinFamily.objects.all().prefetch_related('parent')
        for family in self.families:
            if family.slug == '000':
                self.lookup[0]['000'] = family
                continue
            tree_lvl = len(family.slug.split('_'))
            if tree_lvl > self.tree_depth:
                continue
            if family.slug == '005_001_002':
                continue
            self.lookup[tree_lvl][family.slug] = family

        self.color_mapping = {}
        self.map_family_colors()

        self.proteins = Protein.objects.filter(
            family__slug__startswith="00",
            source__name='SWISSPROT'
            ).prefetch_related(
                'family',
                'family__parent'
                ).order_by('family__slug', 'species_id') #should fix CXCR4

        self.proteins_index = {}
        for p in self.proteins:
            path = p.family.parent.slug
            if not path in self.proteins_index:
                self.proteins_index[path] = []
            self.proteins_index[path].append(p)


    def get_aux_data(self):

        self.aux_data['crystals'] = [x.protein_conformation.protein.parent.id for x in
                                                 Structure.objects.all().exclude(structure_type__slug__startswith='af-')
                                                 .distinct
                                                 ('protein_conformation__protein__parent').prefetch_related('protein_conformation__protein__parent')
                                                 ]

        ligand_data = AssayExperiment.objects.values(
            'protein',
            'protein__entry_name'
            ).annotate(num_ligands=Count('ligand', distinct=True))

        self.aux_data['ligands'] = {
            100 : [x['protein'] for x in ligand_data if x['num_ligands'] <= 100],
            500 : [x['protein'] for x in ligand_data if 100 < x['num_ligands'] <= 500],
            1000 : [x['protein'] for x in ligand_data if 500 < x['num_ligands'] <= 1000],
            2000 : [x['protein'] for x in ligand_data if x['num_ligands'] > 1000] #more than 1000
            }

        ligand_bias_data = BiasedData.objects.filter(
            physiology_biased__isnull=False).values(
            'receptor_id',
            'receptor_id__entry_name'
            ).annotate(num_ligands=Count('ligand_id', distinct=True))

        self.aux_data['ligand_bias'] = {
            10 : [x['receptor_id'] for x in ligand_bias_data if x['num_ligands'] <= 10],
            20 : [x['receptor_id'] for x in ligand_bias_data if 10 < x['num_ligands'] <= 20],
            30 : [x['receptor_id'] for x in ligand_bias_data if 20 < x['num_ligands'] <= 30],
            40 : [x['receptor_id'] for x in ligand_bias_data if x['num_ligands'] > 30] #more than 1000
            }

        ligand_balanced_data = BiasedData.objects.filter(
            pathway_biased__isnull=False).values(
            'receptor_id',
            'receptor_id__entry_name'
            ).annotate(num_ligands=Count('ligand_id', distinct=True))

        self.aux_data['ligand_bias_bal'] = {
            10 : [x['receptor_id'] for x in ligand_balanced_data if x['num_ligands'] <= 10],
            20 : [x['receptor_id'] for x in ligand_balanced_data if 10 < x['num_ligands'] <= 20],
            30 : [x['receptor_id'] for x in ligand_balanced_data if 20 < x['num_ligands'] <= 30],
            40 : [x['receptor_id'] for x in ligand_balanced_data if x['num_ligands'] > 30] #more than 1000
            }

        pathway_pref_data = BiasedData.objects.filter(
            pathway_preferred__isnull=False).values(
            'receptor_id',
            'receptor_id__entry_name'
            ).annotate(num_ligands=Count('ligand_id', distinct=True))

        self.aux_data['pathway_pref'] = {
            10 : [x['receptor_id'] for x in pathway_pref_data if x['num_ligands'] <= 10],
            20 : [x['receptor_id'] for x in pathway_pref_data if 10 < x['num_ligands'] <= 20],
            30 : [x['receptor_id'] for x in pathway_pref_data if 20 < x['num_ligands'] <= 30],
            40 : [x['receptor_id'] for x in pathway_pref_data if x['num_ligands'] > 30] #more than 1000
            }

        subtype_data = BiasedData.objects.filter(
            subtype_biased__isnull=False).values(
            'receptor_id',
            'receptor_id__entry_name'
            ).annotate(num_ligands=Count('ligand_id', distinct=True))

        self.aux_data['subtype'] = {
            10 : [x['receptor_id'] for x in subtype_data if x['num_ligands'] <= 10],
            20 : [x['receptor_id'] for x in subtype_data if 10 < x['num_ligands'] <= 20],
            30 : [x['receptor_id'] for x in subtype_data if 20 < x['num_ligands'] <= 30],
            40 : [x['receptor_id'] for x in subtype_data if x['num_ligands'] > 30] #more than 1000
            }

        subtype_balanced_data = BiasedData.objects.filter(
            pathway_subtype_biased__isnull=False).values(
            'receptor_id',
            'receptor_id__entry_name'
            ).annotate(num_ligands=Count('ligand_id', distinct=True))

        self.aux_data['subtype_bal'] = {
            10 : [x['receptor_id'] for x in subtype_balanced_data if x['num_ligands'] <= 10],
            20 : [x['receptor_id'] for x in subtype_balanced_data if 10 < x['num_ligands'] <= 20],
            30 : [x['receptor_id'] for x in subtype_balanced_data if 20 < x['num_ligands'] <= 30],
            40 : [x['receptor_id'] for x in subtype_balanced_data if x['num_ligands'] > 30] #more than 1000
            }

    def map_family_colors(self):

        for x,y in self.CSS_COLORS.items():
            lvl1_slug = [slug for slug, fam in self.lookup[1].items() if (x[0] == fam.name or x[0] == '')]
            lvl2_slug = []
            for slug, fam in self.lookup[2].items():
                if fam.name.startswith(x[1]) and slug[:3] in lvl1_slug:
                    self.color_mapping[slug] = y

    def get_color(self, slug):
        try:
            return self.color_mapping[slug[:7]]
        except KeyError:
            return 'Black'


    def get_tree_data(self, family):
        """
        Prepare data for coverage diagram. Iterative aproach.
        """
        self.d3_options['branch_length'] = {}
        coverage = PhylogeneticTree(self.root_lvl, self.tree_depth, family)

        for lvl in range(self.root_lvl, self.tree_depth+1):
            if lvl+1 not in self.d3_options['branch_length']:
                self.d3_options['branch_length'][lvl] = ''

            if lvl == self.tree_depth:
                for path, branch in coverage.get_nodes(lvl-2).items():
                    tmp_prots = self.proteins_index[path]
                    for protein in tmp_prots:
                        tmp_node = PhylogeneticTreeNode(
                            protein.entry_name.split("_")[0],
                            self.get_color(protein.family.slug)
                            )
                        if protein.id in self.aux_data['crystals']:
                            tmp_node.increment_value('crystals')
                        for key in self.aux_data['ligands']:
                            if protein.id in self.aux_data['ligands'][key]:
                                tmp_node.increment_value('ligands', key)
                        for key in self.aux_data['ligand_bias']:
                            if protein.id in self.aux_data['ligand_bias'][key]:
                                tmp_node.increment_value('ligand_bias', key)
                        for key in self.aux_data['ligand_bias_bal']:
                            if protein.id in self.aux_data['ligand_bias_bal'][key]:
                                tmp_node.increment_value('ligand_bias_bal', key)
                        for key in self.aux_data['pathway_pref']:
                            if protein.id in self.aux_data['pathway_pref'][key]:
                                tmp_node.increment_value('pathway_pref', key)
                        for key in self.aux_data['subtype']:
                            if protein.id in self.aux_data['subtype'][key]:
                                tmp_node.increment_value('subtype', key)
                        for key in self.aux_data['subtype_bal']:
                            if protein.id in self.aux_data['subtype_bal'][key]:
                                tmp_node.increment_value('subtype_bal', key)
                        coverage.add_data(protein.family.slug, tmp_node)
                return coverage
            children = OrderedDict()
            if lvl+1 in self.SORTED_BRANCHES:
                for slug, node in sorted(self.lookup[lvl+1].items(), key=lambda y: y[1].name.lower()):
                    if node.parent.slug.startswith(family.slug):
                        name = node.name.replace('receptors','').replace('<sub>',' ').replace('</sub>','').strip()
                        children[slug] = PhylogeneticTreeNode(name, self.get_color(node.slug))

                        if len(name) > len(self.d3_options['branch_length'][lvl]):
                            self.d3_options['branch_length'][lvl] = name
            else:
                for slug, node in self.lookup[lvl+1].items():
                    if node.parent.slug.startswith(family.slug):
                        name = node.name.replace('receptors','').replace('<sub>',' ').replace('</sub>','').strip()
                        children[slug] = PhylogeneticTreeNode(name, self.get_color(node.slug))

                        if len(name) > len(self.d3_options['branch_length'][lvl]):
                            self.d3_options['branch_length'][lvl] = name
            for path, data in children.items():
                coverage.add_data(path, data)
        return coverage


    def get_coverage_tree(self, family, coverage=PhylogeneticTreeNode()):
        """
        Prepare data for coverage diagram.
        """
        print('\n'.join([x[1].name for x in coverage.children.items()]))
        tmp_root = len(family.slug.split('_'))

        if tmp_root < self.root_lvl:
            return

        if tmp_root == self.tree_depth:
            tmp_prots = self.proteins.filter(family__parent=family)
            tmp_crystals = self.crystal_proteins.filter(family__parent=family)
            for protein in tmp_prots:
                if tmp_root - len(protein.family.slug.split('_')) == 1:
                    tmp_node = PhylogeneticTreeNode(protein.entry_name.split("_")[0])
                    if self.crystal_proteins.filter(id=protein.id):
                        tmp_node.increment_value("crystalized")

                    coverage.children[protein.family.parent.slug].children[protein.family] = deepcopy(tmp_node)
            return coverage

        if tmp_root+1 in self.SORTED_BRANCHES:
            coverage.children = OrderedDict((x[0], PhylogeneticTreeNode(x[1].name))
                                    for x in sorted(
                                        self.lookup[tmp_root+1].items(),
                                        key=lambda y: y[1].name.lower())
                                    if x[1].parent == family
                                    )
        else:
            coverage.children = OrderedDict({x: PhylogeneticTreeNode(self.lookup[tmp_root+1][x].name)
                                    for x in self.lookup[tmp_root+1]
                                    if self.lookup[tmp_root+1][x].parent == family
                                    })

        for slug, branch in coverage.children.items():
            branch.children = self.get_coverage_tree(self.families.get(slug=slug), deepcopy(coverage)).children

        return coverage
