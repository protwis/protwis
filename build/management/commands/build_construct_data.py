from django.core.management.base import BaseCommand, CommandError
from django.conf import settings
from django.utils.text import slugify
from django.db import IntegrityError

from protein.models import Protein, ProteinConformation
from residue.models import Residue
from structure.models import Structure
from construct.models import (Construct,Crystallization,CrystallizationLigandConc,ChemicalType,Chemical,ChemicalConc,ChemicalList,
CrystallizationMethods,CrystallizationTypes,ChemicalListName,ContributorInfo,ConstructMutation,ConstructInsertion,ConstructInsertionType,
ConstructDeletion,ConstructModification,CrystalInfo,ExpressionSystem,Solubilization,PurificationStep,Purification)

from ligand.models import Ligand, LigandType, LigandRole
from ligand.functions import get_or_make_ligand

from optparse import make_option
import logging
import csv
import os
import json
import datetime

class Command(BaseCommand):
    help = 'Build construct data'

    def add_arguments(self, parser):
        parser.add_argument('--filename', action='append', dest='filename',
            help='Filename to import. Can be used multiple times')

    logger = logging.getLogger(__name__)

        # source file directory
    construct_data_dir = os.sep.join([settings.DATA_DIR, 'structure_data','construct_data'])

    def handle(self, *args, **options):
        if options['filename']:
            filenames = options['filename']
        else:
            filenames = False
        
        try:
            self.purge_construct_data()
            self.create_construct_data(filenames)
        except Exception as msg:
            print(msg)
            self.logger.error(msg)


    def purge_construct_data(self):
        Construct.objects.all().delete()
        Crystallization.objects.all().delete()
        ChemicalConc.objects.all().delete()
        Chemical.objects.all().delete()
        ChemicalType.objects.all().delete()
        ChemicalList.objects.all().delete()
        CrystallizationLigandConc.objects.all().delete()
        CrystallizationMethods.objects.all().delete()
        CrystallizationTypes.objects.all().delete()
        ChemicalListName.objects.all().delete()
        ContributorInfo.objects.all().delete()
        ConstructMutation.objects.all().delete()
        ConstructDeletion.objects.all().delete()
        ConstructInsertion.objects.all().delete()
        ConstructInsertionType.objects.all().delete()
        ConstructModification.objects.all().delete()
        CrystalInfo.objects.all().delete()
        ExpressionSystem.objects.all().delete()
        Solubilization.objects.all().delete()
        Purification.objects.all().delete()
        PurificationStep.objects.all().delete()

    def create_construct_data(self, filenames=False):
        self.logger.info('ADDING EXPERIMENTAL CONSTRUCT DATA')

        # read source files
        if not filenames:
            filenames = os.listdir(self.construct_data_dir)

        for filename in filenames:
            if filename[-4:]!='json':
                continue
            print(filename)
            filepath = os.sep.join([self.construct_data_dir, filename])
            with open(filepath) as json_file:
                d = json.load(json_file)

                protein = Protein.objects.filter(entry_name=d['construct_crystal']['uniprot']).get()
                structure = Structure.objects.filter(pdb_code__index=d['construct_crystal']['pdb']).get()
                protein_conformation = structure.protein_conformation

                construct = Construct()
                construct.protein = protein
                construct.name = d['construct_crystal']['pdb_name']
                construct.json = d

                #CrystalInfo
                crystal = CrystalInfo()
                crystal.resolution = structure.resolution
                crystal.pdb_data = structure.pdb_data
                crystal.pdb_code = structure.pdb_code.index
                crystal.save()

                construct.crystal = crystal

                #Contact INFO
                construct.contributor, created = ContributorInfo.objects.get_or_create(name = d['contact_info']['name_cont'],
                                                               pi_email = d['contact_info']['pi_email'],
                                                               pi_name = d['contact_info']['pi_name'],
                                                               urls = d['contact_info']['url'],
                                                               date = datetime.datetime.strptime(d['contact_info']['date'], '%m/%d/%Y').strftime('%Y-%m-%d'),
                                                               address = d['contact_info']['address'])

                construct.save()
                #MUTATIONS
                for mutation in d['mutations']:
                    mut = ConstructMutation.objects.create(sequence_number=mutation['pos'],wild_type_amino_acid=mutation['wt'],mutated_amino_acid=mutation['mut'])
                    construct.mutations.add(mut)

                #DELETIONS
                insert_deletions = {}
                for deletion in d['deletions']:
                    if 'start' in deletion:
                        dele = ConstructDeletion.objects.create(start=deletion['start'],end=deletion['end'])
                    else:
                        dele = ConstructDeletion.objects.create(start=deletion['pos'],end=deletion['pos'])
                    construct.deletions.add(dele)
                    if deletion['origin']!='user':
                        id = deletion['origin'].split('_')[1]
                        insert_deletions[id] = deletion

                #INSERTIONS (AUX)
                for name,aux in d['auxiliary'].items():
                    id = name.replace('aux','')
                    aux_type,created = ConstructInsertionType.objects.get_or_create(name=aux['type'],subtype=aux['subtype'])
                    insert = ConstructInsertion.objects.create(insert_type=aux_type,presence=aux['presence'],position=aux['position']+"_"+id)

                    if insert.presence == 'YES' and insert.position.startswith('Within Receptor'):
                        #need to fetch range
                        insert.start = insert_deletions[id]['start']
                        insert.end = insert_deletions[id]['end']
                        insert.save()

                    construct.insertions.add(insert)

                #MODIFICATIONS
                for modification in d['modifications']:
                    mod = ConstructModification.objects.create(modification=modification['type'],position_type=modification['position'][0],
                                                               pos_start=modification['position'][1][0],
                                                               pos_end=modification['position'][1][1],remark=modification['remark'] )
                    construct.modifications.add(mod)


                #EXPRESSION
                construct.expression,created = ExpressionSystem.objects.get_or_create(expression_method=d['expression']['expr_method'],
                                                                    host_cell_type=d['expression']['host_cell_type'],
                                                                    host_cell=d['expression']['host_cell'],
                                                                    remarks=d['expression']['expr_remark'])


                
                #solubilization
                c_list = ChemicalList()
                list_name,created  = ChemicalListName.objects.get_or_create(name='Solubilization')
                c_list.name = list_name
                c_list.save()
                ct, created = ChemicalType.objects.get_or_create(name='detergent')
                chem, created = Chemical.objects.get_or_create(name=d['solubilization']['deterg_type'], chemical_type=ct)
                cc, created = ChemicalConc.objects.get_or_create(concentration=d['solubilization']['deterg_concentr'], concentration_unit=d['solubilization']['deterg_concentr_unit'], chemical=chem)
                c_list.chemicals.add(cc)                
                ct, created = ChemicalType.objects.get_or_create(name='additive')
                chem, created = Chemical.objects.get_or_create(name=d['solubilization']['solub_additive'], chemical_type=ct)
                cc, created = ChemicalConc.objects.get_or_create(concentration=d['solubilization']['additive_concentr'], concentration_unit=d['solubilization']['addit_concentr_unit'], chemical=chem)
                c_list.chemicals.add(cc)

                solubilization = Solubilization.objects.create(chemical_list = c_list)

                construct.solubilization = solubilization
                construct.save()

                #Purification
                purification = Purification.objects.create()
                for puri,step in d['solubilization'].items():
                    if not puri.startswith(('chem_enz_treatment','sol_remark')):
                        continue
                    else:
                        s,created = PurificationStep.objects.get_or_create(name=step)
                        purification.steps.add(s)
                        print(step)
                construct.purification = purification
                construct.save()

                #CRYSTALLIZATION 
                c = Crystallization()
                sub_name = "" if 'lcp_lipid' not in d['crystallization'] else d['crystallization']['lcp_lipid']
                c_type, created = CrystallizationTypes.objects.get_or_create(name=d['crystallization']['crystal_type'], sub_name=sub_name)
                c_method, created = CrystallizationMethods.objects.get_or_create(name=d['crystallization']['crystal_method'])
            
                c.crystal_type = c_type
                c.crystal_method = c_method
                c.remarks = d['crystallization']['crystal_remark']
                c.temp = d['crystallization']['temperature']

                if d['crystallization']['ph']=='single_ph':
                    c.ph_start = d['crystallization']['ph_single']
                    c.ph_end = d['crystallization']['ph_single']
                else:
                    c.ph_start = d['crystallization']['ph_range_one']
                    c.ph_end = d['crystallization']['ph_range_two']


                c.protein_conc = d['crystallization']['protein_concentr']
                c.protein_conc_unit = d['crystallization']['protein_conc_unit']

                c.json = d
                c.save()

                #MAKE LISTS
                c_list = ChemicalList()
                list_name,created  = ChemicalListName.objects.get_or_create(name='crystallization_chemical_components')
                c_list.name = list_name
                c_list.save()
                for chemical in d['crystallization']['chemical_components']:
                    ct, created = ChemicalType.objects.get_or_create(name='crystallization_chemical_components')
                    chem, created = Chemical.objects.get_or_create(name=chemical['component'], chemical_type=ct)
                    cc, created = ChemicalConc.objects.get_or_create(concentration=chemical['value'], concentration_unit=chemical['unit'], chemical=chem)
                    c_list.chemicals.add(cc)
                c.chemical_lists.add(c_list)

                if d['crystallization']['crystal_type']=='lipidic cubic phase': #make list of LCP stuff
                    c_list = ChemicalList()
                    # c_list.name = d['crystallization']['lcp_lipid']
                    list_name,created  = ChemicalListName.objects.get_or_create(name='LCP')
                    c_list.name = list_name
                    c_list.save()
                    ct, created = ChemicalType.objects.get_or_create(name='LCP Lipid additive')
                    chem, created = Chemical.objects.get_or_create(name=d['crystallization']['lcp_add'], chemical_type=ct)
                    cc, created = ChemicalConc.objects.get_or_create(concentration=d['crystallization']['lcp_conc'], concentration_unit=d['crystallization']['lcp_conc_unit'], chemical=chem)
                    c_list.chemicals.add(cc)
                    c.chemical_lists.add(c_list)

                #DETERGENT
                c_list = ChemicalList()
                list_name,created  = ChemicalListName.objects.get_or_create(name='Detergent')
                c_list.name = list_name
                c_list.save()
                ct, created = ChemicalType.objects.get_or_create(name='detergent')
                chem, created = Chemical.objects.get_or_create(name=d['crystallization']['detergent'], chemical_type=ct)
                cc, created = ChemicalConc.objects.get_or_create(concentration=d['crystallization']['deterg_conc'], concentration_unit=d['crystallization']['deterg_conc_unit'], chemical=chem)
                c_list.chemicals.add(cc)
                c.chemical_lists.add(c_list)

                #LIPID
                c_list = ChemicalList()
                list_name,created  = ChemicalListName.objects.get_or_create(name='Lipid')
                c_list.name = list_name
                c_list.save()
                ct, created = ChemicalType.objects.get_or_create(name='lipid')
                chem, created = Chemical.objects.get_or_create(name=d['crystallization']['lipid'], chemical_type=ct)
                cc, created = ChemicalConc.objects.get_or_create(concentration=d['crystallization']['lipid_concentr'], concentration_unit=d['crystallization']['lipid_concentr_unit'], chemical=chem)
                c_list.chemicals.add(cc)
                c.chemical_lists.add(c_list)



                #Use ligand function to get ligand if it exists or otherwise create. Lots of checks for inchi/smiles/name
                ligand = get_or_make_ligand(d['construct_crystal']['ligand_id'],d['construct_crystal']['ligand_id_type'],d['construct_crystal']['ligand_name'])

                if ligand and 'ligand_activity' in d['construct_crystal']:
                    role_slug = slugify(d['construct_crystal']['ligand_activity'])
                    try:
                        lr, created = LigandRole.objects.get_or_create(slug=role_slug,
                        defaults={'name': d['construct_crystal']['ligand_activity']})
                    except IntegrityError:
                        lr = LigandRole.objects.get(slug=role_slug)

                ligand_c = CrystallizationLigandConc()
                ligand_c.construct_crystallization = c
                ligand_c.ligand = ligand
                ligand_c.ligand_role = lr
                if 'ligand_conc' in d['construct_crystal']:
                    ligand_c.ligand_conc = d['construct_crystal']['ligand_conc']
                if 'ligand_conc_unit' in d['construct_crystal']:
                    ligand_c.ligand_conc_unit = d['construct_crystal']['ligand_conc_unit']
                ligand_c.save()

                c.ligands.add(ligand_c)

                construct.crystallization = c


                construct.save()



                #print(json_data)


        self.logger.info('COMPLETED CREATING EXPERIMENTAL CONSTRUCT DATA')