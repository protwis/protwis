from django.core.management.base import BaseCommand, CommandError
from django.conf import settings
from django.utils.text import slugify
from django.utils.html import strip_tags

from protein.models import (Protein, ProteinConformation, ProteinState, ProteinSequenceType, ProteinSegment,
ProteinFusion, ProteinFusionProtein, ProteinSource)
from residue.models import Residue
from construct.models import *

from optparse import make_option
from datetime import datetime
import logging, os
import yaml


class Command(BaseCommand):
    help = 'Reads source data and creates protein records for constructs'
    
    def add_arguments(self, parser):
        parser.add_argument('--filename', action='append', dest='filename',
            help='Filename to import. Can be used multiple times')
        parser.add_argument('--purge', action='store_true', dest='purge', default=False,
            help='Purge existing construct records')

    logger = logging.getLogger(__name__)

    # source file directory
    construct_data_dir = os.sep.join([settings.DATA_DIR, 'structure_data', 'constructs'])

    def handle(self, *args, **options):
        # delete any existing construct data
        if options['purge']:
            try:
                self.purge_constructs()
            except Exception as msg:
                print(msg)
                self.logger.error(msg)
        
        # import the structure data
        try:
            self.create_constructs(options['filename'])
        except Exception as msg:
            print(msg)
            self.logger.error(msg)
    
    def purge_constructs(self):
        try:
            pst = ProteinSequenceType.objects.get(slug='mod')
            Protein.objects.filter(sequence_type=pst).delete()
        except ProteinSequenceType.DoesNotExist:
            self.logger.warning('ProteinSequenceType mod not found: nothing to delete.')

    def create_constructs(self, filenames):
        self.logger.info('CREATING CONSTRUCTS')
        
        # what files should be parsed?
        if not filenames:
            filenames = os.listdir(self.construct_data_dir)

        # parse files
        for source_file in filenames:
            source_file_path = os.sep.join([self.construct_data_dir, source_file])
            if os.path.isfile(source_file_path) and source_file[0] != '.':
                self.logger.info('Reading file {}'.format(source_file_path))
                # read the yaml file
                with open(source_file_path, 'r') as f:
                    sd = yaml.load(f, Loader=yaml.FullLoader)

                    # is a protein specified?
                    if 'protein' not in sd:
                        self.logger.error('Protein not specified for construct, skipping')
                        continue

                    # fetch the parent protein
                    try:
                        ppc = ProteinConformation.objects.select_related('protein__family', 'protein__species',
                            'protein__residue_numbering_scheme').get(protein__entry_name=sd['protein'],
                            state__slug=settings.DEFAULT_PROTEIN_STATE)
                    except ProteinConformation.DoesNotExist:
                        # abort if parent protein is not found
                        self.logger.error('Parent protein {} for construct {} not found, aborting!'.format(
                            sd['protein'], sd['name']))
                        continue

                    if not Protein.objects.filter(name=sd['name']).exists():
                        # create a protein record
                        p = Protein()
                        p.parent = ppc.protein
                        p.family = ppc.protein.family
                        p.species = ppc.protein.species
                        p.residue_numbering_scheme = ppc.protein.residue_numbering_scheme
                        p.sequence_type, created = ProteinSequenceType.objects.get_or_create(slug='mod',
                            defaults={'name': 'Modified'})
                        p.source, created = ProteinSource.objects.get_or_create(name='OTHER')
                        p.entry_name = slugify(strip_tags(sd['name']))
                        p.name = sd['name']
                        p.sequence = ppc.protein.sequence
                        # save protein (construct)
                        try:
                            p.save()
                            self.logger.info('Created construct {} with parent protein {}'.format(p.name,
                                ppc.protein.entry_name))
                        except Exception as e:
                            print(e)
                            self.logger.error('Failed creating construct {} with parent protein {}'.format(p.name,
                                ppc.protein.entry_name))
                            continue
                    else:
                        p = Protein.objects.get(name=sd['name'])

                    if not ProteinConformation.objects.filter(protein=p).exists():
                        # create protein conformation record
                        pc = ProteinConformation()
                        pc.protein = p
                        pc.state = ProteinState.objects.get(slug=settings.DEFAULT_PROTEIN_STATE)
                        try:
                            pc.save()
                            self.logger.info('Created conformation {} of protein {}'.format(pc.state.name, p.name))
                        except:
                            self.logger.error('Failed creating conformation {} of protein {}'.format(pc.state.name,
                                p.entry_name))

                    # # create residue records
                    # deletions = []
                    # deletions_list = []
                    # if 'deletions' in sd and sd['deletions']:
                    #     for t in sd['deletions']:
                    #         deletions += list(range(t[0],t[1]+1))
                    #         deletions_list.append(str(t[0])+'-'+str(t[1])) 
                    # s = ","
                    # deletion_string = s.join(deletions_list)
                         

                    # mutations = {}
                    # if 'mutations' in sd and sd['mutations']:
                    #     for m in sd['mutations']:
                    #         res_num = m[1:-1]
                    #         mutations[res_num] = {
                    #             'wt_res': m[0],
                    #             'mut_res': m[-1],
                    #             'full': m,
                    #         }
                    
                    # # Create construct record
                    # c = Construct()            
                    # c.protein_conformation = pc
                    # c.deletions =  deletion_string
                    # c.save()
                      

                    # Create Auxiliary proteins
#                    if 'auxiliary_proteins' in sd and sd['auxiliary_proteins']:
#                        ap = AuxProtein()
#                        ap.construct = c
#                        apct = AuxProteinType.objects.create()
                       # ap.protein_type = apct 
#                        apct.save()
#                        if 'remarks' in sd['auxiliary_proteins']:
#                            ap.remarks = sd['auxiliary_proteins']['remarks']
#                        ap.save()
 

#                        for step in sd['auxiliary_proteins']:
#                            if 'type' in step and 'name' in step and'sequence' in step:
#                                ap.protein_type = apct
                 #              ap.protein_type, created = AuxProteinType.objects.get_or_create()
#                                ap.name = sd['auxiliary_proteins']['name']
#                                ap.uniprot_id = sd['auxiliary_proteins']['uniprot_id']
#                                ap.sequence = sd['auxiliary_proteins']['sequence']
                                #mutations if any to be included from mutation model along with reason of mutation
#                                ap.position = sd['auxiliary_proteins']['position']
#                                ap.deletions = sd['auxiliary_proteins']['deletions']
                                
#                            else:
#                                self.logger.error('Auxiliary protein step incorrectly defined for {}'.format(p))



                   #   # create expression records
                   #  if 'expression_sys' in sd and sd['expression_sys']:
                   #      ce = ConstructExpression()           
                   #      ce.construct = c
                   #      ce.expression_system, created = ConstructExpressionSystem.objects.get_or_create(expression_method = sd['expression_sys']['expression_method'], host_cell_type = sd['expression_sys']['host_cell_type'], host_cell = sd['expression_sys']['host_cell'])
                   #      if 'remarks' in sd:
                   #          ce.remarks = sd['expression_sys']['remarks']
                   #      ce.save()
 
               
                   #  # create solubilization records
                   #  if 'solubilization' in sd and sd['solubilization'] and 'steps' in sd['solubilization'] and sd['solubilization']['steps']:
                   #      so = ConstructSolubilization()
                   #      so.construct = c
                   #      cl = ChemicalList.objects.create()
                   #      so.chemical_list = cl 

                   #      for step in sd['solubilization']['steps']:
                   #          if 'type' in step and 'item' in step and'concentration' in step:
                   #              chem = Chemical()
                   #              chem.chemical_type,  created = ChemicalType.objects.get_or_create(name = step['type']) 
                   #              chem.name =  step['item']
                   #              chem.save()

                   #              cc = ChemicalConc()
                   #              cc.concentration = step['concentration']
                   #              cc.chemical = chem    # since ChemicalConc has a ForeignKey to Chemical
                   #              cc.save()
                   #              cl.chemicals.add(cc)                          
                   #          else:
                   #              self.logger.error('Solubilization step incorrectly defined for {}'.format(p))                                 

                   #      if 'remarks' in sd['solubilization']:
                   #          so.remarks = sd['solubilization']['remarks']
                   #      so.save()



                   #  # create  purification records
                   #  if 'purification' in sd and sd['purification'] and sd['purification']['steps']:
                   #      pu = ConstructPurification()
                   #      pu.construct = c
                   #      if 'remarks' in sd['purification']:
                   #          pu.remarks = sd['purification']['remarks']
                   #      pu.save() 
                   #      for step in sd['purification']['steps']:
                   #          if 'type' in step and 'description' in step:
                   #              pust = PurificationStep()
                   #              pust.description = step['description']
                   #              pust.purification = pu
                   #              pust.purification_type, created = PurificationStepType.objects.get_or_create(name = step['type'] ) # 2 values returned by get_or_create
                   #              if created: 
                   #                  self.logger.info('Created purification step type {}'.format(pust.purification_type))
                   #              pust.save()

                   #          else:
                   #              self.logger.error('Purification step incorrectly defined for {}'.format(p))

                        


                   # # create crystallization records
                   #  if 'crystallization' in sd and sd['crystallization']: 
                   #      cy = ConstructCrystallization()
                   #      cy.construct = c
                   #      cyt = CrystallizationMethodTypes.objects.create()
                   #      cy.crystal_type = cyt
                   #      cy.method = sd['crystallization']['method']
                   #      cy.settings = sd['crystallization']['settings']
                   #      cy.protein_conc = sd['crystallization']['protein_conc']
                   #      cl = ChemicalList.objects.create()
                   #      cy.chemical_list = cl   

                   #      for step in sd['crystallization']['chemicallist']:
                   #          if 'type' in step and 'item' in step and'concentration' in step:
                   #              chem = Chemical()
                   #              chem.chemical_type,  created = ChemicalType.objects.get_or_create(name = step['type']) 
                                
                   #              chem.name =  step['item']
                   #              chem.save()
                   #              cc = ChemicalConc()
                   #              cc.concentration = step['concentration']
                   #              cc.chemical = chem    # since ChemicalConc has a ForeignKey to Chemical
                   #              cc.save()
                                
                   #              cl.chemicals.add(cc)                          
                   #          else:
                   #              self.logger.error('Crystallization step incorrectly defined for {}'.format(p))                        

                   #      cy.aqueous_solution_lipid_ratio = sd['crystallization']['aqueous_solution_lipid_ratio_LCP']
                   #      cy.lcp_bolus_volume = sd['crystallization']['LCP_bolus_volume']
                   #      cy.precipitant_solution_volume = sd['crystallization']['precipitant_solution_volume']
                   #      cy.temp = sd['crystallization']['temperature']
                   #      cy.ph = sd['crystallization']['ph']  


                   #      if 'remarks' in sd['crystallization']:
                   #          cy.remarks = sd['crystallization']['remarks']
                   #      cy.save()

                                     
                   #  # fusion proteins
                   #  split_segments = {}
                   #  if 'fusion_proteins' in sd and sd['fusion_proteins']:
                   #      for fp in sd['fusion_proteins']:
                   #          fp_start = Residue.objects.get(protein_conformation=ppc,
                   #              sequence_number=fp['positions'][0])
                   #          fp_end = Residue.objects.get(protein_conformation=ppc, sequence_number=fp['positions'][1])
                   #          # if the fusion protein is inserted within only one segment (the usual case), split that
                   #          # segment into two segments
                   #          if fp_start and fp_start.protein_segment == fp_end.protein_segment:
                   #              # get/create split protein segments
                   #              segment_before, created = ProteinSegment.objects.get_or_create(
                   #                  slug=fp_start.protein_segment.slug+"_1", defaults={
                   #                  'name': fp_start.protein_segment.name,
                   #                  'category': fp_start.protein_segment.category,
                   #                  'partial': True})
                   #              segment_after, created = ProteinSegment.objects.get_or_create(
                   #                  slug=fp_start.protein_segment.slug+"_2", defaults={
                   #                  'name': fp_start.protein_segment.name,
                   #                  'category': fp_start.protein_segment.category,
                   #                  'partial': True})

                   #              # keep track of  information about split segments
                   #              split_segments[fp_start.protein_segment.slug] = {
                   #                  'start': {
                   #                      'sequence_number': fp['positions'][0],
                   #                      'segment': segment_before,
                   #                  },
                   #                  'end': {
                   #                      'sequence_number': fp['positions'][1],
                   #                      'segment': segment_after,
                   #                  },
                   #              }

                   #          # get/insert fusion protein
                   #          fusion, create = ProteinFusion.objects.get_or_create(name=fp['name'], defaults={
                   #              'sequence': fp['sequence']})

                   #          # create relationship with protein
                   #          ProteinFusionProtein.objects.create(protein=p, protein_fusion=fusion,
                   #              segment_before=segment_before, segment_after=segment_after)

                    # prs = Residue.objects.filter(protein_conformation=ppc).prefetch_related(
                    #     'protein_conformation__protein', 'protein_segment', 'generic_number',
                    #     'display_generic_number__scheme', 'alternative_generic_numbers__scheme')
                    # updated_sequence = ''
                    # for pr in prs:
                    #     if pr.sequence_number not in deletions:
                    #         r = Residue()
                    #         r.protein_conformation = pc
                    #         r.generic_number = pr.generic_number
                    #         r.display_generic_number = pr.display_generic_number
                    #         r.sequence_number = pr.sequence_number
                            
                    #         # check for split segments
                    #         if pr.protein_segment.slug in split_segments:
                    #             rsns = split_segments[pr.protein_segment.slug]['start']['sequence_number']
                    #             rsne = split_segments[pr.protein_segment.slug]['end']['sequence_number']
                    #             if r.sequence_number <= rsns:
                    #                 r.protein_segment = split_segments[pr.protein_segment.slug]['start']['segment']
                    #             elif r.sequence_number >= rsne:
                    #                 r.protein_segment = split_segments[pr.protein_segment.slug]['end']['segment']
                    #         else:
                    #             r.protein_segment = pr.protein_segment

                    #         # amino acid, check for mutations
                    #         if r.sequence_number in mutations:
                    #             if mutations[r.sequence_number]['wt_res'] == pr.amino_acid:
                    #                 r.amino_acid = mutations[r.sequence_number]['mut_res']
                    #             else:
                    #                 self.logger.error('Mutation {} in construct {} does not match wild-type sequence' \
                    #                     + ' of {}'.format(mutations[r.sequence_number]['full'], pc.protein.name,
                    #                     ppc.protein.entry_name))
                    #         else:
                    #             r.amino_acid = pr.amino_acid

                    #         # save amino acid to updated sequence
                    #         updated_sequence += r.amino_acid

                    #         # save residue before populating M2M relations
                    #         r.save()

                    #         # alternative generic numbers
                    #         agns = pr.alternative_generic_numbers.all()
                    #         for agn in agns:
                    #             r.alternative_generic_numbers.add(agn)
                    
                    # # update sequence
                    # p.sequence = updated_sequence
                    # p.save()

        self.logger.info('COMPLETED CREATING CONSTRUCTS')
