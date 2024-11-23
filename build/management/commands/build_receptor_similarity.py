from build.management.commands.base_build import Command as BaseBuild

from django.db.models import F,Q


from protein.models import Protein, ProteinSegment, ProteinFamily, Species

from common.alignment import Alignment

from collections import OrderedDict
import os
import sys
import logging
import re

import csv

from datetime import datetime, date
import time



starttime = datetime.now()
logger = logging.getLogger('class_similarity')
hdlr = logging.FileHandler('./logs/class_similarity.log')
formatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s')
hdlr.setFormatter(formatter)
logger.addHandler(hdlr)
logger.setLevel(logging.INFO)

class_prefix_re = re.compile(r'^(Class)\s+', flags=re.I)

from protein.models import CLASSLESS_PARENT_GPCR_SLUGS

FIELDNAMES = ['class1','receptor1_name_short',
                'receptor1_name','receptor1_entry_name',
                'class2','receptor2_name_short',
                'receptor2_name','receptor2_entry_name',
                'identity', 'similarity']

matrix_header_max_length = 20

build_date = date.today()

import warnings
warnings.filterwarnings("ignore")

class Command(BaseBuild):
    help = 'Build cross-class receptor similarity and identity.'
    

    def add_arguments(self, parser):
        super(Command, self).add_arguments(parser=parser)
        parser.add_argument('--output',type=str, help="Output file path. Default 'gpcr_cross_class_similarity_data.csv'.", default='gpcr_cross_class_similarity_data.csv', action='store')
        parser.add_argument('--verbose', help='Prints progress in stdout.', default=False, action='store_true')
        parser.add_argument('--force', help="Overwrites output file. If --output is not used, overwrites the file 'gpcr_similarity_data.csv'.", default=False, action='store_true')
        parser.add_argument('--limit',type=int, help='Use only any indicated number of GPCRs per class.', default=False, action='store')

    def get_parent_gpcr_families(self,exclude_classless_artificial_class=True,include_classless_natural_classes=True):
        parent_family = ProteinFamily.objects.get(slug='000') 
        parent_gpcr_families = ProteinFamily.objects.filter(parent_id=parent_family.pk, slug__startswith='0').exclude(pk=parent_family.pk)
        if exclude_classless_artificial_class:
            for slug in CLASSLESS_PARENT_GPCR_SLUGS:
                parent_gpcr_families = parent_gpcr_families.exclude(slug__startswith=slug)
        if include_classless_natural_classes:
            classless_protein_families = self.get_classless_bottom_protein_families()
        parent_gpcr_families = list(parent_gpcr_families)+classless_protein_families
        return sorted(parent_gpcr_families,key=lambda f: (int(f.slug.split('_')[0])))
    
    def get_human_species(self):
        return Species.objects.get(common_name__iexact='Human')
    
    def get_yeast_species(self):
        return Protein.objects.filter(entry_name__iendswith='_yeast')[0].species
    
    def __slug_tree_branch(self, slug_parts,slug_tree_dict):
        if len(slug_parts) == 1:
            slug_tree_dict[slug_parts[0]] = None
            return slug_tree_dict
        else:
            if slug_parts[0] in slug_tree_dict:
                slug_subtree_dict = slug_tree_dict[slug_parts[0]]
                if slug_subtree_dict is None:
                    slug_subtree_dict = {}
            else:
                slug_subtree_dict = {}
            slug_tree_dict[slug_parts[0]] = self.__slug_tree_branch(slug_parts[1:],slug_subtree_dict)
        return slug_tree_dict

    def __sort_slug_tree_branch(self, slug_tree_dict):
        slug_tree_ordered_dict = OrderedDict()
        for slug in sorted(sorted(slug_tree_dict.keys(),key = lambda x: int(x[1:])),key = lambda x: x[0]):
            slug_subtree_dict = slug_tree_dict[slug]
            if slug_subtree_dict is not None:
                slug_tree_ordered_dict[slug] = self.__sort_slug_tree_branch(slug_subtree_dict)
            else:
                slug_tree_ordered_dict[slug] = None
        return slug_tree_ordered_dict
    
    def __parse_slug_tree_(self, slug_tree_dict,slug_list_list,slug_list):
        for slug,subtree in slug_tree_dict.items():
            slug_list.append(slug)
            if subtree is not None:
                self.__parse_slug_tree_(subtree,slug_list_list,slug_list)
            else:
                slug_list_list.append(slug_list.copy())
            slug_list.pop()
    
    def get_classless_bottom_protein_families(self):
        classless_parent_gpcrs_slugs_list = sorted(sorted(CLASSLESS_PARENT_GPCR_SLUGS,key = lambda x: int(x[1:])),key = lambda x: x[0])
        parent_gpcr_families = ProteinFamily.objects.filter(slug__startswith=classless_parent_gpcrs_slugs_list[0])
        
        for slug in classless_parent_gpcrs_slugs_list[1:]:
            parent_gpcr_families = parent_gpcr_families.filter(slug__startswith=slug)
        
        slug_2_family_dict = {}
        for f in parent_gpcr_families:
            slug_2_family_dict[f.slug] = f
        family_slug_tree_dict = {}
        for slug in slug_2_family_dict.keys():
            slug_parts = slug.split('_')
            self.__slug_tree_branch(slug_parts,family_slug_tree_dict)
        family_slug_tree_ordered_dict = self.__sort_slug_tree_branch(family_slug_tree_dict)
        slug_list_list = []
        slug_list = []
        self.__parse_slug_tree_(family_slug_tree_ordered_dict,slug_list_list,slug_list)
        classless_bottom_slugs_list = ['_'.join(slug_list) for slug_list in slug_list_list]
        return [slug_2_family_dict[slug] for slug in classless_bottom_slugs_list]
    
    def filter_out_non_species_parent_gpcr_families(self,parent_gpcr_families,species):
        """ Filters out parent GPCR families as a list of ProteinFamily objects that belong to a species.
            parent_gpcr_families: a list of protein.ProteinFamily objects
            species: protein.Species object
        """
        new_parent_gpcr_families_slugs_set = set()
        species2 = species
        try:
            species_iterator = iter(species)
        except TypeError as te:
             species2 = [species]
        parent_gpcr_families_slugs = []
        for family in parent_gpcr_families:
            parent_gpcr_families_slugs.append(family.slug)

        for slug in parent_gpcr_families_slugs:
            q = Protein.objects.annotate(family_slug=F('family__slug')).filter(family_slug__startswith=slug,species__in=species2)
            if q.exists():
                new_parent_gpcr_families_slugs_set.add(slug)
        return [f for f in parent_gpcr_families if f.slug in new_parent_gpcr_families_slugs_set]
        
    
    def filter_out_non_human_parent_gpcr_families(self,parent_gpcr_families):
        return self.filter_out_non_species_parent_gpcr_families(parent_gpcr_families,self.get_human_species())

    def filter_out_non_yeast_parent_gpcr_families(self,parent_gpcr_families):
        return self.filter_out_non_species_parent_gpcr_families(parent_gpcr_families,self.get_yeast_species())
    
    def handle(self, *args, **options):
        if os.path.lexists(options['output']) and not options['force']:
            print('Path "'+options['output']+' already exists. Aborting.')
            exit(1)
        with open(options['output'],'w') as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=FIELDNAMES)
            writer.writeheader()

            initial_step1 = 380 #If alignment fails, please, set this to a lower value
            initial_step2 = 380 #If alignment fails, please, set this to a lower value
            start_time = time.time()

            parent_families = self.get_parent_gpcr_families(exclude_classless_artificial_class=True,include_classless_natural_classes=True)
            human_parent_gpcr_families = self.filter_out_non_human_parent_gpcr_families(parent_families)
            human_species = self.get_human_species()



            yeast_parent_gpcr_families = self.filter_out_non_yeast_parent_gpcr_families(parent_families)
            yeast_non_human_parent_gpcr_families_set = set(yeast_parent_gpcr_families) - set(human_parent_gpcr_families)
            yeast_non_human_parent_gpcr_families = [gpcr_family for gpcr_family in yeast_parent_gpcr_families if gpcr_family in yeast_non_human_parent_gpcr_families_set]
            yeast_species = self.get_yeast_species()

            gpcr_segments = ProteinSegment.objects.filter(Q(proteinfamily='GPCR') & (Q(slug__regex='TM[1-7]') | Q(slug='H8')))


            human_parent_gpcr_families_protein = {}
            human_parent_gpcr_families_protein_num = {}
            for gpcr_class in human_parent_gpcr_families:
                gpcr_class_proteins = Protein.objects.all().annotate(family_slug=F('family__slug')).filter(species=human_species,family_slug__startswith=gpcr_class.slug)
                gpcr_class_proteins = gpcr_class_proteins.exclude(accession=None).order_by('family_slug','entry_name')
                human_parent_gpcr_families_protein[gpcr_class] = list(gpcr_class_proteins)
                human_parent_gpcr_families_protein_num[gpcr_class] = len(gpcr_class_proteins)

            yeast_non_human_parent_gpcr_families_protein = {}
            yeast_non_human_parent_gpcr_families_protein_num = {} 
            # for gpcr_class in yeast_non_human_parent_gpcr_families:
            #     gpcr_class_proteins = Protein.objects.all().annotate(family_slug=F('family__slug')).filter(species=yeast_species,family_slug__startswith=gpcr_class.slug)
            #     gpcr_class_proteins = gpcr_class_proteins.exclude(accession=None).order_by('family_slug','entry_name')
            #     yeast_non_human_parent_gpcr_families_protein[gpcr_class] = list(gpcr_class_proteins)
            #     yeast_non_human_parent_gpcr_families_protein_num[gpcr_class] = len(gpcr_class_proteins)

            # selected_parent_gpcr_families = human_parent_gpcr_families + yeast_non_human_parent_gpcr_families
            selected_parent_gpcr_families = human_parent_gpcr_families
            selected_parent_gpcr_families_protein = {}
            selected_parent_gpcr_families_protein_num = {}
            for gpcr_class in selected_parent_gpcr_families:
                if gpcr_class in human_parent_gpcr_families_protein and gpcr_class in yeast_non_human_parent_gpcr_families_protein:
                    selected_parent_gpcr_families_protein[gpcr_class] = human_parent_gpcr_families_protein[gpcr_class] + yeast_non_human_parent_gpcr_families_protein[gpcr_class]
                    selected_parent_gpcr_families_protein_num[gpcr_class] = human_parent_gpcr_families_protein_num[gpcr_class] + yeast_non_human_parent_gpcr_families_protein_num[gpcr_class]
                elif gpcr_class in human_parent_gpcr_families_protein:
                    selected_parent_gpcr_families_protein[gpcr_class] = human_parent_gpcr_families_protein[gpcr_class]
                    selected_parent_gpcr_families_protein_num[gpcr_class] = human_parent_gpcr_families_protein_num[gpcr_class]
                elif gpcr_class in yeast_non_human_parent_gpcr_families_protein:
                    selected_parent_gpcr_families_protein[gpcr_class] = yeast_non_human_parent_gpcr_families_protein[gpcr_class]
                    selected_parent_gpcr_families_protein_num[gpcr_class] = yeast_non_human_parent_gpcr_families_protein_num[gpcr_class]



            step1=int(initial_step1) #If alignment fails, please, set this to a lower value
            step2=int(initial_step2) #If alignment fails, please, set this to a lower value
            step_halved = False 

            unique_keys_set = set()
            for gpcr_class in selected_parent_gpcr_families:
                gpcr_class1_name = class_prefix_re.sub(r'',gpcr_class.name.replace('<i>','').replace('</i>',''))
                while_loop_continue = False
                while True:
                    if options['limit']: 
                        if selected_parent_gpcr_families_protein_num[gpcr_class] < options['limit']:
                            protein_num1 = selected_parent_gpcr_families_protein_num[gpcr_class]
                        else:
                            protein_num1 = options['limit']
                    else:
                        protein_num1 = selected_parent_gpcr_families_protein_num[gpcr_class]
                    
                    for clim in range(0,protein_num1,step1):
                        
                        gpcr_class_proteins = selected_parent_gpcr_families_protein[gpcr_class]
                        gpcr_class_proteins = gpcr_class_proteins[clim:clim+step1]
                        if options['limit']:
                            gpcr_class_proteins = gpcr_class_proteins[:options['limit']]
                        
                        for gpcr_class2 in selected_parent_gpcr_families:
                            gpcr_class2_name = class_prefix_re.sub(r'',gpcr_class2.name.replace('<i>','').replace('</i>',''))
                            unique_list = [gpcr_class1_name,gpcr_class2_name]
                            unique_list.sort()
                            unique_key = '@'.join(unique_list)
                            if unique_key in unique_keys_set:
                                continue
                            if options['limit']: 
                                if selected_parent_gpcr_families_protein_num[gpcr_class2] < options['limit']:
                                    protein_num2 = selected_parent_gpcr_families_protein_num[gpcr_class2]
                                else:
                                    protein_num2 = options['limit']
                            else:
                                protein_num2 = selected_parent_gpcr_families_protein_num[gpcr_class2]


                            for clim2 in range(0,protein_num2,step2):
                                if options['verbose']:
                                        print(gpcr_class,"from:"+str(clim+1),"to:"+str(clim+step1),"(of:{})".format(protein_num1),\
                                        'vs',gpcr_class2,"from:"+str(clim2+1),"to:"+str(clim2+step2),"(of:{})".format(protein_num2))
                                gpcr_class2_proteins = selected_parent_gpcr_families_protein[gpcr_class2]
                                gpcr_class2_proteins = gpcr_class2_proteins[clim2:clim2+step2]

                                if options['limit']:
                                    gpcr_class2_proteins = gpcr_class2_proteins[:options['limit']]

                                if gpcr_class1_name == gpcr_class2_name and clim2 == clim and clim+step1 == clim2+step2:
                                    proteins = gpcr_class_proteins
                                else:
                                    proteins = gpcr_class_proteins + gpcr_class2_proteins

                                cs_alignment = Alignment()
                                cs_alignment.load_proteins(proteins)
                                cs_alignment.load_segments(gpcr_segments)
                                build_alignment_return_value = cs_alignment.build_alignment()
                                if build_alignment_return_value == "Too large":
                                    print('Alignment too large. Retrying...', file=sys.stderr)
                                    while_loop_continue = True
                                    break
                                cs_alignment.remove_non_generic_numbers_from_alignment() 
                                cs_alignment.calculate_similarity_matrix()

                                entry_name_2_position = {}
                                pos = 0
                                for protein in cs_alignment.proteins:
                                    entry_name_2_position[protein.protein.entry_name] = pos
                                    pos += 1

                                for protein1 in gpcr_class_proteins:
                                    for protein2 in gpcr_class2_proteins:
                                        if protein1.entry_name == protein2.entry_name:
                                            continue
                                        pos_s = entry_name_2_position[protein1.entry_name]
                                        similarity_value = int(cs_alignment.similarity_matrix[protein2.entry_name]['values'][pos_s][0])
                                        pos_i = entry_name_2_position[protein2.entry_name]
                                        identity_value = int(cs_alignment.similarity_matrix[protein1.entry_name]['values'][pos_i][0])
                                        protein_tuple = (protein1,protein2)
                                        similarity_value = similarity_value
                                        identity_value = identity_value
                                        protein1_name_short = Protein(name=protein1.name).short().replace('<i>','').replace('</i>','')
                                        protein1_name = protein1.name.replace('<i>','').replace('</i>','')
                                        protein2_name_short = Protein(name=protein2.name).short().replace('<i>','').replace('</i>','')
                                        protein2_name = protein2.name.replace('<i>','').replace('</i>','')
                                        receptor1_entry_name = protein1.entry_name
                                        receptor2_entry_name = protein2.entry_name
                                        data = {'class1':gpcr_class1_name,'receptor1_name_short':protein1_name_short,
                                                'receptor1_name':protein1_name,'receptor1_entry_name':receptor1_entry_name,
                                                'class2':gpcr_class2_name,'receptor2_name_short':protein2_name_short,
                                                'receptor2_name': protein2_name,'receptor2_entry_name':receptor2_entry_name,
                                                'identity':identity_value, 'similarity': similarity_value}
                                        writer.writerow(data)
                                csvfile.flush()

                                        
                            unique_keys_set.add(unique_key)
                            if while_loop_continue:
                                break
                        try:
                            del proteins
                        except Exception:
                            pass
                        try:
                            del gpcr_class_proteins
                        except Exception:
                            pass
                        if while_loop_continue:
                            break
                        selected_parent_gpcr_families_protein[gpcr_class][clim:clim+step1] = [None for i1 in range(clim,clim+step1)]
                    if while_loop_continue:
                        if options['verbose']: print("Halving step1 and step2...") 
                        step1 = step1 // 2
                        step2 = step2 // 2
                        step_halved = True
                        while_loop_continue = False
                        if options['verbose']: print("Retrying last class similarity computation with the new steps...")
                        continue
                    elif step_halved:
                        print("Restoring initial step1 and step2...")
                        step1=initial_step1
                        step2=initial_step2
                        step_halved = False
                    del selected_parent_gpcr_families_protein[gpcr_class]
                    break    
                

        end_time = time.time()
        elapsed_time = end_time - start_time
        if options['verbose']:
            print('Execution time:', elapsed_time, 'seconds')
            print('Done.')

     





             



