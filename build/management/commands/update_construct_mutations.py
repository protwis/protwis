from django.utils.text import slugify
from django.db import IntegrityError
from django.core.management.base import BaseCommand, CommandError
from django.core.management import call_command
from django.conf import settings
from django.db import connection
from django.core.cache import cache
from common.tools import test_model_updates
from build.management.commands.build_ligand_functions import get_or_create_ligand

from structure.models import Structure
from construct.functions import  fetch_pdb_info
from construct.models import *
from residue.models import Residue

from ligand.models import Ligand, LigandType, LigandRole

from operator import itemgetter
from itertools import groupby
import datetime
import logging
from optparse import make_option
import os
import xlrd
import json
import django.apps

import time
from collections import OrderedDict


class Command(BaseCommand):
    help = 'Update construct mutations from excel file'

    logger = logging.getLogger(__name__)

    path = os.sep.join([settings.DATA_DIR, 'structure_data', 'construct_data', 'Stabilising_Mutations_In_Xtal_Constructs.xlsx'])
    annotation_file = os.sep.join([settings.DATA_DIR, 'structure_data', 'construct_data', 'construct_annotations.xlsx'])
    tracker = {}
    all_models = django.apps.apps.get_models()[6:]
    test_model_updates(all_models, tracker, initialize=True)

    def handle(self, *args, **options):

        ## DEBUG ITEMS
        # self.match_all_with_uniprot_mutations()
        # Simply check deletions on record vs newest pdb
        # self.check_deletions()
        # # Make sure json file is correct
        # # self.json_check_for_mutations_deletions()

        ### EXPORTS ###
        # # Export current inserts to file
        # self.export_inserts()
        ########################

        self.all_pdbs = []

        # Restart constructs
        self.rebuild_constructs()

        ## MUST HAVE ##
        self.excel_mutations = self.parse_excel(self.path,'Mutation_Data')
        self.check_mutations()
        # changes deletions to match PDB
        # Custom rules exist in the function
        self.replace_deletions()

        ## IMPORTS ###
        self.import_inserts()
        self.import_expression()
        self.import_solub()
        self.import_puri()
        self.import_xtal()
        test_model_updates(self.all_models, self.tracker, check=True)

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

    def rebuild_constructs(self):
        self.purge_construct_data()
        structures = Structure.objects.all().exclude(structure_type__slug__startswith='af-')
        for s in structures:
            pdbname = str(s)
            cache.delete(pdbname+"_auto_d")
            self.all_pdbs.append(pdbname)
            protein_conformation = s.protein_conformation

            construct = Construct()
            construct.protein = protein_conformation.protein.parent
            construct.name = pdbname
            construct.json = ''
            construct.structure = s

            #CrystalInfo
            crystal = CrystalInfo()
            crystal.resolution = s.resolution
            crystal.pdb_data = s.pdb_data
            crystal.pdb_code = s.pdb_code.index
            crystal.save()

            construct.crystal = crystal

            d = {}
            d['contact_info'] = {}
            d['contact_info']['name_cont'] = 'gpcrdb'
            d['contact_info']['pi_email'] = 'info@gpcrdb.org'
            d['contact_info']['pi_name'] = 'gpcrdb'
            d['contact_info']['url'] = 'gpcrdb.org'
            d['contact_info']['date'] = time.strftime('%m/%d/%Y')
            d['contact_info']['address'] = ''
            construct.contributor, created = ContributorInfo.objects.get_or_create(name = d['contact_info']['name_cont'],
                                                       pi_email = d['contact_info']['pi_email'],
                                                       pi_name = d['contact_info']['pi_name'],
                                                       urls = d['contact_info']['url'],
                                                       date = datetime.datetime.strptime(d['contact_info']['date'], '%m/%d/%Y').strftime('%Y-%m-%d'),
                                                       address = d['contact_info']['address'])

            construct.save()

    def import_expression(self):
        expressions = self.parse_excel(self.annotation_file,'Expression')
        exp_list = {}
        for e in expressions:
            if e[1] not in exp_list:
                exp_list[e[1]] = []
            exp_list[e[1]].append(e)
            # print("#####",e[1])
            for construct in Construct.objects.filter(structure__pdb_code__index=e[1]).all():
                ce = construct.expression
                 #print(construct.name)
                # if ce:
                #     print(ce.expression_method,ce.host_cell_type,ce.host_cell,ce.expression_time,ce.remarks)
                #     print(e)

                new_e = ExpressionSystem()
                new_e.expression_method = e[4]
                new_e.host_cell_type = e[2]
                new_e.host_cell = e[3]
                new_e.expression_time = e[5]
                new_e.remarks = e[6]

                # if ce:
                #     ce.delete()

                new_e.save()
                construct.expression = new_e
                construct.save()

        missing = list(set(self.all_pdbs) - set(exp_list.keys()))
        print(sorted(missing)," do not have any Expression annotated -- add them to sheet with NONE if they have none")


    def import_solub(self):
        solubs = self.parse_excel(self.annotation_file,'Solubzn-Detergent')
        solub_list = {}
        for s in solubs:
            if s[1] not in solub_list:
                solub_list[s[1]] = []
            solub_list[s[1]].append(s)
        for pdb,chems in solub_list.items():
            c_list = ChemicalList()
            list_name,created  = ChemicalListName.objects.get_or_create(name='Solubilization')
            c_list.name = list_name
            c_list.save()
            for c in chems:
                ct, created = ChemicalType.objects.get_or_create(name=c[3])
                chem, created = Chemical.objects.get_or_create(name=c[2], chemical_type=ct)
                cc, created = ChemicalConc.objects.get_or_create(concentration=c[4], concentration_unit=c[5], chemical=chem)
                c_list.chemicals.add(cc)

            solubilization = Solubilization.objects.create(chemical_list = c_list)
            try:
                construct = Construct.objects.get(structure__pdb_code__index=c[1].upper())
                construct.solubilization = solubilization
                construct.save()
            except:
                print(c[1],'cannot find pdb construct')

        missing = list(set(self.all_pdbs) - set(solub_list.keys()))
        print(sorted(missing)," do not have any Solubzn-Detergent annotated -- add them to sheet with NONE if they have none")


    def import_puri(self):
        puris = self.parse_excel(self.annotation_file,'Purification-Treatment')
        puri_list = {}

        for p in puris:
            if p[1] not in puri_list:
                puri_list[p[1]] = []
            puri_list[p[1]].append(p)

        for pdb,treat in puri_list.items():
            try:
                construct = Construct.objects.get(structure__pdb_code__index=pdb.upper())
            except:
                print(pdb,'cannot find pdb construct')
                continue

            purification = Purification.objects.create()
            # print('purification made')
            for t in treat:
                s,created = PurificationStep.objects.get_or_create(name=t[2])
                # print(t,s,c)
                purification.steps.add(s)

            construct.purification = purification
            construct.save()

        missing = list(set(self.all_pdbs) - set(puri_list.keys()))
        print(sorted(missing)," do not have any Purification-Treatment annotated -- add them to sheet with NONE if they have none")

    def import_xtal(self):
        xtals = self.parse_excel(self.annotation_file,'Xtal-methods')
        xtals_list = {}
        for x in xtals:
            if x[1] in xtals_list:
                print('pdbcode duplicate?',x[1])
            xtals_list[x[1]] = x

        xtal_chems = self.parse_excel(self.annotation_file,'Xtal-Chemicals')
        xtals_chems_list = {}
        for x in xtal_chems:
            if x[1] not in xtals_chems_list:
                xtals_chems_list[x[1]] = []
            xtals_chems_list[x[1]].append(x)

        xtal_ligands = self.parse_excel(self.annotation_file,'PDB-ligand-complex')
        xtal_ligands_list = {}
        for x in xtal_ligands:
            if x[1] not in xtal_ligands_list:
                xtal_ligands_list[x[1]] = []
            xtal_ligands_list[x[1]].append(x)


        missing = list(set(self.all_pdbs) - set(xtals_list.keys()))
        print(sorted(missing)," do not have any Xtal-methods annotated -- add them to sheet with NONE if they have none")
        missing = list(set(self.all_pdbs) - set(xtals_chems_list.keys()))
        print(sorted(missing)," do not have any Xtal-Chemicals annotated -- add them to sheet with NONE if they have none")
        missing = list(set(self.all_pdbs) - set(xtal_ligands_list.keys()))
        print(sorted(missing)," do not have any PDB-ligand-complex annotated -- add them to sheet with NONE if they have none")

        for pdb,x in xtals_list.items():
            try:
                construct = Construct.objects.get(structure__pdb_code__index=pdb.upper())
            except:
                print(pdb,'cannot find pdb construct')
                continue
            c = Crystallization()
            c_type, created = CrystallizationTypes.objects.get_or_create(name=x[3], sub_name=x[4])
            c_method, created = CrystallizationMethods.objects.get_or_create(name=x[2])

            c.crystal_type = c_type
            c.crystal_method = c_method
            c.remarks = x[10]
            c.temp = x[7]

            # Some entries have it wrong here
            try:
                c.ph_start = float(x[8])
                c.ph_end = float(x[9])
            except:
                c.ph_start = 0
                c.ph_end = 0

            c.protein_conc = x[5]
            c.protein_conc_unit = x[6]
            c.save()

            if pdb in xtals_chems_list:
                chems = xtals_chems_list[pdb]
                chem_types = {}
                for chem in chems:
                    ctype = chem[2]
                    if ctype not in chem_types:
                        chem_types[ctype] = []
                    chem_types[ctype].append(chem)

                for ctype,chems in chem_types.items():
                    c_list = ChemicalList()
                    list_name,created  = ChemicalListName.objects.get_or_create(name=ctype)
                    c_list.name = list_name
                    c_list.save()
                    for ch in chems:
                        ct, created = ChemicalType.objects.get_or_create(name=ch[4])
                        chem, created = Chemical.objects.get_or_create(name=ch[3], chemical_type=ct)
                        cc, created = ChemicalConc.objects.get_or_create(concentration=ch[5], concentration_unit=ch[6], chemical=chem)
                        c_list.chemicals.add(cc)
                    c.chemical_lists.add(c_list)
            else:
                print('no chems for ',pdb)

            if pdb in xtal_ligands_list:
                l = xtal_ligands_list[pdb][0]
                ligand = get_or_create_ligand(l[2])
                # ligand = get_or_create_ligand(l[7],l[6],l[2])
                role_slug = slugify(l[3])
                try:
                    lr, created = LigandRole.objects.get_or_create(slug=role_slug,
                    defaults={'name': l[3]})
                except IntegrityError:
                    lr = LigandRole.objects.get(slug=role_slug)
                if ligand:
                    ligand_c = CrystallizationLigandConc()
                    ligand_c.construct_crystallization = c
                    ligand_c.ligand = ligand
                    if lr:
                        ligand_c.ligand_role = lr
                    if l[4]:
                        ligand_c.ligand_conc = l[4]
                    if l[5]:
                        ligand_c.ligand_conc_unit = l[5]
                    ligand_c.save()

                    c.ligands.add(ligand_c)

            construct.crystallization = c
            construct.save()


    def import_inserts(self):

        # Delete current
        ConstructInsertion.objects.all().delete()

        inserts = self.parse_excel(self.annotation_file,'inserts')

        tracked_pdbs = []
        for i in inserts:

            if i[1] not in tracked_pdbs:
                tracked_pdbs.append(i[1])

            if i[2]=='NONE':
                # SKip those that are entries just to show there is nothing
                continue
            if i[3]=='?':
                continue
            print(i)
            aux_type, created = ConstructInsertionType.objects.get_or_create(name=i[5],subtype=i[6])
            for construct in Construct.objects.filter(structure__pdb_code__index=i[1].upper()):
                try:
                    if i[3]=='':
                        i[3] = 0
                    insert = ConstructInsertion.objects.create(construct=construct, insert_type=aux_type,presence=i[7],position=i[2]+"_"+str(int(i[3])))
                except Exception as e:
                    print('Error with insert! FIXIT',i,str(e))
                if i[4]:
                    i[4] = str(i[4])
                    #if position information add that
                    if len(i[4].split(":"))>1:
                        insert.start = i[4].split(":")[0]
                        insert.end = i[4].split(":")[1]
                    else:
                        insert.start = int(i[4].split('.')[0])
                        insert.end = int(i[4].split('.')[0])
                insert.save()
                construct.invalidate_schematics()

        missing = list(set(self.all_pdbs) - set(tracked_pdbs))
        print(sorted(missing)," do not have any inserts annotated -- add them to sheet with NONE if they have none")


    def export_inserts(self):

        ## Get a list of annotated structures
        inserts = self.parse_excel(self.annotation_file,'inserts')
        tracked_pdbs = []
        for i in inserts:
            if i[1] not in tracked_pdbs:
                tracked_pdbs.append(i[1])

        checked_pdbcode_without_fusions = ['4AMI']
        checked_pdb_code_with_fusions = {'3PDS' : 'icl3'}
        # 3pds has P00720 LYS

        constructs = Construct.objects.all().order_by('protein__entry_name')
        list_of_comfirmed_fusion = ['C8TP59','Q0SXH8','Q9V2J8','Soluble cytochrome b562','Endolysin','Rubredoxin','Lysozyme','Flavodoxin','GlgA glycogen synthase']
        csv_rows = []

        csv_rows_expression = []
        csv_rows_solub = []
        csv_rows_puri = []
        csv_rows_xtal = []
        csv_rows_xtal_chems = []
        csv_rows_xtal_ligands = []
        for c in constructs:
            pdbname = c.structure.pdb_code.index
            if pdbname in tracked_pdbs:
                # print("skip",pdbname,"already annotated.")
                continue
            else:
                print("dump auto annotation for",pdbname," (considering adding these to excel sheet)")

            fusion_position, fusions, linkers = c.fusion()
            protein = Protein.objects.filter(entry_name=pdbname.lower()).get()
            uniprot = protein.parent.entry_name
            #print(c.name,fusion_position)

            # if pdbname != '4UHR':
            #     continue

            d = json.loads(c.json)
            # print(d)
            # print(pdbname)

            if c.crystallization:
                if c.crystallization.crystal_type.sub_name == 'other [See next field]':
                    print(pdbname,'Has other subname',c.crystallization.crystal_type.sub_name)
                    print(d['crystallization']['lcp_lipid'],d['crystallization']['other_lcp_lipid'])
                    c_type, created = CrystallizationTypes.objects.get_or_create(name=d['crystallization']['crystal_type'], sub_name=d['crystallization']['other_lcp_lipid'])
                    c.crystallization.crystal_type = c_type
                    c.crystallization.save()

            if 'crystallization' in d and 'crystal_type' in d['crystallization'] and (d['crystallization']['crystal_type']=='lipidic cubic phase' or d['crystallization']['crystal_type']=='lipidic cubic phase (LCP)'): #make list of LCP stuff
                if not c.crystallization.chemical_lists.filter(name__name='LCP').exists() or (c.crystallization.chemical_lists.filter(name__name='LCP').exists() and c.crystallization.chemical_lists.filter(name__name='LCP').first().chemicals.count() == 0) :
                    print(pdbname,"has LCP stuff to do!")
                    if not 'lcp_add' in d['crystallization']:
                        print('lcp_add missing')
                    if not 'lcp_conc' in d['crystallization']:
                        print('lcp_conc missing')
                    if not 'lcp_conc_unit' in d['crystallization']:
                        print('lcp_conc_unit missing')
                    try:
                        c.crystallization.chemical_lists.filter(name__name='LCP').delete()
                        c_list = ChemicalList()
                        # c_list.name = d['crystallization']['lcp_lipid']
                        list_name,created  = ChemicalListName.objects.get_or_create(name='LCP')
                        c_list.name = list_name
                        c_list.save()
                        ct, created = ChemicalType.objects.get_or_create(name='LCP Lipid additive')
                        chem, created = Chemical.objects.get_or_create(name=d['crystallization']['lcp_add'], chemical_type=ct)
                        cc, created = ChemicalConc.objects.get_or_create(concentration=d['crystallization']['lcp_conc'], concentration_unit=d['crystallization']['lcp_conc_unit'], chemical=chem)
                        c_list.chemicals.add(cc)
                        c.crystallization.chemical_lists.add(c_list)
                    except:
                        print('error for ',pdbname)

            if 'crystallization' in d and 'detergent' in d['crystallization']:
                if d['crystallization']['detergent'] == 'other [See next field]':
                    print(pdbname,'Has other detergent!',d['crystallization']['other_deterg'])
                    c.crystallization.chemical_lists.filter(name__name='Detergent').delete()
                    c_list = ChemicalList()
                    list_name,created  = ChemicalListName.objects.get_or_create(name='Detergent')
                    c_list.name = list_name
                    c_list.save()
                    ct, created = ChemicalType.objects.get_or_create(name='detergent')
                    chem, created = Chemical.objects.get_or_create(name=d['crystallization']['other_deterg'], chemical_type=ct)
                    cc, created = ChemicalConc.objects.get_or_create(concentration=d['crystallization']['deterg_conc'], concentration_unit=d['crystallization']['deterg_conc_unit'], chemical=chem)
                    c_list.chemicals.add(cc)
                    c.crystallization.chemical_lists.add(c_list)

            if 'crystallization' in d and 'lipid' in d['crystallization']:
                if d['crystallization']['lipid'] == 'other [See next field]':
                    print(pdbname,'Has other Lipid!',d['crystallization']['other_lipid'])
                    c.crystallization.chemical_lists.filter(name__name='Lipid').delete()
                    c_list = ChemicalList()
                    list_name,created  = ChemicalListName.objects.get_or_create(name='Lipid')
                    c_list.name = list_name
                    c_list.save()
                    ct, created = ChemicalType.objects.get_or_create(name='lipid')
                    chem, created = Chemical.objects.get_or_create(name=d['crystallization']['other_lipid'], chemical_type=ct)
                    cc, created = ChemicalConc.objects.get_or_create(concentration=d['crystallization']['lipid_concentr'], concentration_unit=d['crystallization']['lipid_concentr_unit'], chemical=chem)
                    c_list.chemicals.add(cc)
                    c.crystallization.chemical_lists.add(c_list)
            d = cache.get(pdbname+"_deletions")
            # d = None
            if not d:
                d = fetch_pdb_info(pdbname,protein,ignore_gasper_annotation=True)
                cache.set(pdbname+"_deletions",d,60*60*24)
            #print(d['auxiliary'])
            found = False
            found_where = None
            found_type = None
            for aux, v in d['auxiliary'].items():
                # print(v['subtype'])
                if v['subtype'] in list_of_comfirmed_fusion:
                    found = True
                    # print(v['subtype'],v['position'])
                    found_where = v['position']
                    found_type = v['subtype']
            # print(d['construct_sequences'])
            for aux, v in d['construct_sequences'].items():
                if aux in list_of_comfirmed_fusion:
                    found = True
                    # print(aux,v['where'])
                    found_where = v['where']
                    found_type = aux

            if not found and fusion_position:
                # print(fusion_position, "not found in", c.name)
                #print(d['construct_sequences'])
                pass


            for i in c.insertions.all().order_by('position'):
                #print(i)
                position = i.position.split("_")
                seq_pos = ''

                if i.insert_type.subtype == 'Soluble cytochrome b562':
                    i.insert_type.subtype = 'BRIL (b562RIL)'
                    i.insert_type.name = 'fusion'
                if i.insert_type.subtype == 'Endolysin':
                    i.insert_type.subtype = 'T4 Lysozyme (T4L)'
                    i.insert_type.name = 'fusion'
                if i.insert_type.subtype == 'Rubredoxin':
                    i.insert_type.name = 'fusion'
                if i.insert_type.subtype == 'Flavodoxin':
                    i.insert_type.name = 'fusion'
                if i.insert_type.subtype == 'GlgA glycogen synthase':
                    i.insert_type.subtype = 'PGS (Pyrococcus abyssi glycogen synthase)'
                    i.insert_type.name = 'fusion'

                if i.start:
                    seq_pos = '%s:%s' % (i.start,i.end)
                if i.insert_type.name=='fusion' or i.insert_type.subtype in list_of_comfirmed_fusion:
                    insert = [uniprot,pdbname,position[0],position[1],seq_pos,i.insert_type.name,i.insert_type.subtype,i.presence,found_where,found_type]
                else:
                    insert = [uniprot,pdbname,position[0],position[1],seq_pos,i.insert_type.name,i.insert_type.subtype,i.presence]
                #print(insert)
                csv_rows.append(insert)

            if c.insertions.count()==0:
                insert = [uniprot,pdbname,'NONE','','','','','',found_where,found_type]
                #print(insert)
                csv_rows.append(insert)

            # # EXPRESSION #
            # if c.expression:
            #     csv_rows_expression.append([uniprot,pdbname,c.expression.host_cell_type,c.expression.host_cell,c.expression.expression_method,c.expression.expression_time,c.expression.remarks.strip()])
            # else:
            #     csv_rows_expression.append([uniprot,pdbname,'MISSING'])

            # SOLUB #
            if c.solubilization:
                remarks = c.solubilization.remarks
                for chem in c.solubilization.chemical_list.chemicals.all():
                    chem_name = chem.chemical.name
                    chem_type = chem.chemical.chemical_type.name
                    chem_conc = chem.concentration
                    chem_unit = chem.concentration_unit
                    csv_rows_solub.append([uniprot,pdbname,chem_name,chem_type,chem_conc,chem_unit,remarks])
                if c.solubilization.chemical_list.chemicals.count()==0:
                    csv_rows_solub.append([uniprot,pdbname,'MISSING'])
            else:
                csv_rows_solub.append([uniprot,pdbname,'MISSING'])

            # PURI #
            if c.purification:
                remarks = c.purification.remarks
                for step in c.purification.steps.all():
                    step_name = step.name
                    step_description = step.description
                    csv_rows_puri.append([uniprot,pdbname,step_name,step_description,remarks])
                if c.solubilization.chemical_list.chemicals.count()==0:
                    csv_rows_puri.append([uniprot,pdbname,'MISSING'])
            else:
                csv_rows_puri.append([uniprot,pdbname,'MISSING'])

            # XTAL #
            if c.crystallization:

                crystal_type_name = c.crystallization.crystal_type.name
                crystal_type_sub_name = c.crystallization.crystal_type.sub_name

                crystal_method = c.crystallization.crystal_method.name

                remarks = c.crystallization.remarks
                if remarks:
                    remars = remarks.encode('utf8')
                    remarks = remarks.replace("\n","")
                    remarks = remarks.replace("\r","")
                    remarks = remarks.replace("\t","")
                    import re
                    remarks = re.sub(' +',' ',remarks)

                protein_conc = c.crystallization.protein_conc
                protein_conc_unit = c.crystallization.protein_conc_unit

                temp = c.crystallization.temp

                ph_start = c.crystallization.ph_start
                ph_end = c.crystallization.ph_end
                csv_rows_xtal.append([uniprot,pdbname,crystal_type_name,crystal_type_sub_name,crystal_method,protein_conc,protein_conc_unit,temp,ph_start,ph_end,remarks])

                # XTAL CHEMS #
                if c.crystallization.chemical_lists:
                    for clist in c.crystallization.chemical_lists.all():
                        list_name = clist.name
                        for chem in clist.chemicals.all():
                            chem_name = chem.chemical.name
                            chem_type = chem.chemical.chemical_type.name
                            chem_conc = chem.concentration
                            chem_unit = chem.concentration_unit
                            csv_rows_xtal_chems.append([uniprot,pdbname,list_name,chem_name,chem_type,chem_conc,chem_unit])
                else:
                    csv_rows_xtal_chems.append([uniprot,pdbname,'MISSING'])

                # XTAL LIGAND #
                if c.crystallization.ligands:
                    for lig in c.crystallization.ligands.all():
                        #print('ligand',lig)
                        l = lig.ligand
                        l_name = l.name
                        smiles = l.smiles
                        inchi = l.inchikey

                        id_type = None
                        id_index = None
                        for links in l.ids.all():
                            #Just use the first link
                            id_type = links.web_resource.slug
                            id_index = links.index
                            break


                        l_role = lig.ligand_role.name
                        l_conc = lig.ligand_conc
                        l_conc_unit = lig.ligand_conc_unit
                        if not id_type:
                            if inchi:
                                id_type = 'inchikey'
                                id_index = inchi
                            elif smiles:
                                id_type = 'smiles'
                                id_index = smiles

                        csv_rows_xtal_ligands.append([uniprot,pdbname,l_name,l_role,l_conc,l_conc_unit,id_type,id_index])
                else:
                    csv_rows_xtal_ligands.append([uniprot,pdbname,'MISSING'])


            else:
                csv_rows_xtal.append([uniprot,pdbname,'MISSING'])
                csv_rows_xtal_chems.append([uniprot,pdbname,'MISSING'])
                csv_rows_xtal_ligands.append([uniprot,pdbname,'MISSING'])




        import csv
        with open('construct_inserts.csv', 'w') as f:
            writer = csv.writer(f, delimiter = '\t')
            writer.writerows(csv_rows)
        with open('construct_expression.csv', 'w') as f:
            writer = csv.writer(f, delimiter = '\t')
            writer.writerows(csv_rows_expression)
        with open('construct_solub.csv', 'w') as f:
            writer = csv.writer(f, delimiter = '\t')
            writer.writerows(csv_rows_solub)
        with open('construct_puri.csv', 'w') as f:
            writer = csv.writer(f, delimiter = '\t')
            writer.writerows(csv_rows_puri)
        with open('construct_xtal.csv', 'w') as f:
            writer = csv.writer(f, delimiter = '\t')
            writer.writerows(csv_rows_xtal)
        with open('construct_xtal_chems.csv', 'w') as f:
            writer = csv.writer(f, delimiter = '\t')
            writer.writerows(csv_rows_xtal_chems)
        with open('construct_xtal_ligs.csv', 'w') as f:
            writer = csv.writer(f, delimiter = '\t')
            writer.writerows(csv_rows_xtal_ligands)

    def parse_excel(self,path, sheet = None):
        workbook = xlrd.open_workbook(path)
        worksheets = workbook.sheet_names()
        d = []
        for worksheet_name in worksheets:
            if worksheet_name in d:
                print('Error, worksheet with this name already loaded')
                continue
            if sheet and worksheet_name != sheet:
                # Only run this sheet
                continue

            #d[worksheet_name] = OrderedDict()
            worksheet = workbook.sheet_by_name(worksheet_name)

            num_rows = worksheet.nrows - 1
            num_cells = worksheet.ncols
            curr_row = 0 #skip first, otherwise -1

            headers = []
            for i in range(num_cells):
                h = worksheet.cell_value(0, i)
                if h=="":
                    h = "i_"+str(i)
                if h in headers:
                    h += "_"+str(i)
                headers.append(worksheet.cell_value(0, i))
            for curr_row in range(1,num_rows+1):
                row = worksheet.row(curr_row)
                key = worksheet.cell_value(curr_row, 0)

                # if key=='':
                #     continue
                # if key not in d:
                #     d[key] = []
                temprow = OrderedDict()
                temprow = []
                for curr_cell in range(num_cells):
                    cell_value = worksheet.cell_value(curr_row, curr_cell)
                    temprow.append(cell_value)
                    # if headers[curr_cell] not in temprow:
                    #     temprow[headers[curr_cell]] = cell_value
                d.append(temprow)
        return d



    def fix_wrong_mutations(self):
        muts = ConstructMutation.objects.all()
        for mut in muts:
            # mut.save()
            if mut.residue == None:
                print(mut,'no link')

    def check_deletions(self):
        constructs = Construct.objects.all()
        csv_rows = []
        for c in constructs:
            issues = []
            pdbname = c.structure.pdb_code.index
            cname = c.name
            # if pdbname!='4GPO':
            #     continue
            protein = Protein.objects.filter(entry_name=pdbname.lower()).get()
            uniprot = protein.parent.entry_name
            d = cache.get(pdbname+"_auto_d")
            # d = None
            if not d:
                d = fetch_pdb_info(c_pdb,protein,ignore_gasper_annotation=True)
                cache.set(pdbname+"_auto_d",d,60*60*24)
            pdb_deletions = []
            for d in d['deletions']:
                pdb_deletions += range(d['start'],d['end']+1)
            cons_dels = c.deletions.all()
            db_deletions = []
            for d in cons_dels:
                db_deletions += range(d.start,d.end+1)

            present_in_pdb_only = set(pdb_deletions)-set(db_deletions)
            present_in_pdb_only_list = []
            for k, g in groupby(enumerate(present_in_pdb_only), lambda x:x[0]-x[1]):
                group = list(map(itemgetter(1), g))
                present_in_pdb_only_list.append([group[0], group[-1]])

            present_in_db_only = set(db_deletions)-set(pdb_deletions)
            present_in_db_only_list = []
            for k, g in groupby(enumerate(present_in_db_only), lambda x:x[0]-x[1]):
                group = list(map(itemgetter(1), g))
                present_in_db_only_list.append([group[0], group[-1]])

            if present_in_pdb_only or present_in_db_only:
                print(pdbname)
            if present_in_pdb_only: print("PDBONLY",present_in_pdb_only)
            if present_in_db_only: print("DBONLY",present_in_db_only)
            csv_rows.append([pdbname,uniprot,cname,present_in_db_only_list,present_in_pdb_only_list,''])
        import csv
        with open('construct_del_issues.csv', 'w') as f:
            writer = csv.writer(f, delimiter = '\t')
            writer.writerows(csv_rows)

    def replace_deletions(self):
        # delete alle deletions
        # ConstructDeletion.objects.all().delete()
        for c in Construct.objects.all():

            pdbname = c.structure.pdb_code.index

            # if not pdbname in ['5F8U','2VT4']:
            #     continue
            # print(pdbname)

            #reset caches
            c.schematics = None
            c.snakecache = None
            c.save()

            c.deletions.all().delete()

            pdbname = c.structure.pdb_code.index
            cname = c.name
            protein = Protein.objects.filter(entry_name=pdbname.lower()).get()
            uniprot = protein.parent.entry_name
            d = cache.get(pdbname+"_auto_d")
            # d = None
            if d == None or (not d and 'deletions' in d):
                d = fetch_pdb_info(c_pdb,protein,ignore_gasper_annotation=True)
                cache.set(pdbname+"_auto_d",d,60*60*24)
            if 'deletions' in d:
                for d in d['deletions']:
                    dele, created = ConstructDeletion.objects.get_or_create(construct=c, start=d['start'],end=d['end'])
            else:
                print('No deletions in d[]',pdbname)

    def json_check_for_mutations_deletions(self):
        for c in Construct.objects.all():
            print(c.name)
            d = json.loads(c.json)

            items_to_delete = ['mut_aa_','wt_aa_','aa_no_','mut_type_','delet_start_','delet_origin_','deletion_']

            #clean json
            if 'raw_data' not in d:
                d['raw_data'] = {}
            else:
                # remove points associated with mutations and deletions
                # FIXME also do inserts / modifications
                delete = []
                for k,v in d['raw_data'].items():
                    if k.startswith(tuple(items_to_delete)):
                        print('delete',k)
                        delete.append(k)
                for k in delete:
                    del(d['raw_data'][k])

            # reset these lists to recreate them
            print(d['mutations'])
            d['mutations'] = []
            d['deletions'] = []

            cons_muts = ConstructMutation.objects.filter(construct = c)
            i = 2
            for m in cons_muts:
                d['mutations'].append({'mut':m.mutated_amino_acid,'wt':m.wild_type_amino_acid,'pos':m.sequence_number,'type':str(m.effects.all())})
                d['raw_data']['mut_aa_'+str(i)] = m.mutated_amino_acid
                d['raw_data']['wt_aa_'+str(i)] = m.wild_type_amino_acid
                d['raw_data']['aa_no_'+str(i)] = m.sequence_number
                d['raw_data']['mut_type_'+str(i)] = str(m.effects.all())
                i+=1
            print(d['mutations'])
    def match_all_with_uniprot_mutations(self):
        constructs = Construct.objects.all()
        csv_rows = [['reference','pdb','construct name','class','lig type','rec fam','uniprot','segment','gpcrdb#','AA no.','WT aa','Mut aa','','','','','','Remark']]
        for c in constructs:
            issues = []
            pdbname = c.structure.pdb_code.index
            # if pdbname!='4GPO':
            #     continue
            protein = Protein.objects.filter(entry_name=pdbname.lower()).get()
            uniprot = protein.parent.entry_name
            # if uniprot !='glp1r_human':
            #     continue
            d = cache.get(pdbname+"_auto_d")
            # d = None
            if not d:
                d = fetch_pdb_info(c_pdb,protein,ignore_gasper_annotation=True)
                cache.set(pdbname+"_auto_d",d,60*60*24)
            # print('pdb',d['mutations'])
            cons_muts = ConstructMutation.objects.filter(construct = c)
            for m in cons_muts:
                seq_pos = m.sequence_number
                found = False
                for pdb_m in d['mutations']:
                    if int(pdb_m['pos']) == seq_pos and pdb_m['wt']==m.wild_type_amino_acid and (pdb_m['mut']==m.mutated_amino_acid):
                        found = True
                        break
                if not found:
                    ignore = False
                    for m_xlx in self.excel_mutations:
                        if m_xlx[6]==uniprot and int(m_xlx[9])==seq_pos and m_xlx[11]==m.mutated_amino_acid and m_xlx[10]==m.wild_type_amino_acid:
                            found = False
                            if pdbname in m_xlx[1] or m_xlx[1]=='':
                                found = True
                            if '%'+pdbname in m_xlx[1]:
                                found = False
                            if found:
                                if  m_xlx[16]!='Non-receptor' and m_xlx[16]!='Wrong annotation - remove!':
                                    ignore = True
                    if ignore:
                        issues.append(('In excel but missing in pdb?',seq_pos,m.mutated_amino_acid))
                    else:
                        issues.append(('Not in excel nor pdb, deleting',seq_pos,m.mutated_amino_acid))
                        m.delete()
                    print(issues)
                    continue

                    issues.append(('missing in pdb?',seq_pos,m.mutated_amino_acid))
                    mut_aa = m.mutated_amino_acid
                    pos = m.sequence_number
                    wt_aa = m.wild_type_amino_acid
                    annotated_effect = "Not identified in PDB -- perhaps delete?"

                    res = Residue.objects.get(protein_conformation__protein=protein.parent, sequence_number=pos)
                    seg = res.protein_segment.slug
                    if res.generic_number:
                        gn = res.generic_number.label
                    else:
                        gn = ''
                    csv_rows.append(['',pdbname,'','','','',protein.parent.entry_name,seg,gn,pos,wt_aa,mut_aa,'','','','','',annotated_effect])

            for m in d['mutations']:
                cons_muts = ConstructMutation.objects.filter(construct = c, sequence_number = m['pos'], mutated_amino_acid = m['mut'], wild_type_amino_acid = m['wt'])
                if not cons_muts.exists():
                    # print('missing',m)
                    ignore = False
                    for m_xlx in self.excel_mutations:
                        if m_xlx[6]==uniprot and int(m_xlx[9])==m['pos']:
                            found = False
                            if pdbname in m_xlx[1] or m_xlx[1]=='':
                                found = True
                            if '%'+pdbname in m_xlx[1]:
                                found = False
                            if found:
                                if  m_xlx[16]=='Non-receptor' or m_xlx[16]=='Wrong annotation - remove!':
                                    ignore = True
                    if ignore:
                        continue
                    issues.append(('{}{}{} ({})'.format(m['wt'],m['pos'],m['mut'],m['type']),' not in db, nor to be ignored in excel'))
                    mut_aa = m['mut']
                    pos = m['pos']
                    wt_aa = m['wt']
                    annotated_effect = m['type']
                    res = Residue.objects.get(protein_conformation__protein=protein.parent, sequence_number=pos)
                    seg = res.protein_segment.slug
                    if res.generic_number:
                        gn = res.generic_number.label
                    else:
                        gn = ''
                    csv_rows.append(['',pdbname,'','','','',protein.parent.entry_name,seg,gn,pos,wt_aa,mut_aa,'','','','','',annotated_effect])
            if issues:
                print(pdbname)
                for i in issues:
                    print(i)
        import csv
        with open('construct_mut_issues.csv', 'w') as f:
            writer = csv.writer(f, delimiter = '\t')
            writer.writerows(csv_rows)


    def check_mutations(self):
        track_annotated_mutations = []
        cached_mutations = {}
        mut_pdb_list = {}
        for i,mut in enumerate(self.excel_mutations):
            # print("Progress ",i,len(self.excel_mutations))
            # continue
            #print(mut)
            m = {}
            m['gn'] = mut[8]
            m['mut_aa'] = mut[11]
            m['wt_aa'] = mut[10]
            m['entry_name'] = mut[6]
            m['pos'] = int(mut[9])
            m['pdb'] = mut[1]
            m['thermo_effect'] = mut[12]
            m['expression_effect'] = mut[13]
            m['site_effect'] = mut[14]
            m['site_effect_type'] = mut[15]
            m['other_effect'] = mut[16]
            # if m['entry_name']!='glp1r_human':
            #     continue


            # print(m)
            if m['pdb'] and m['pdb'][0]!='%':
                pdbs = m['pdb'].split(',')
                # print(pdbs)
                for pdb in pdbs:
                    if pdb not in mut_pdb_list:
                        mut_pdb_list[pdb] = []
                if len(pdbs)==0:
                    cons = Construct.objects.filter(structure__pdb_code__index = m['pdb'])
                else:
                    cons = Construct.objects.filter(structure__pdb_code__index__in = pdbs)
            else:
                cons = Construct.objects.filter(structure__protein_conformation__protein__parent__entry_name = m['entry_name'])

            pdbs_has = []
            pdbs_hasnot = []
            for c in cons:
                c_pdb = c.structure.pdb_code.index
                if c_pdb not in mut_pdb_list:
                    mut_pdb_list[c_pdb] = []
                not_to_check = None
                if m['pdb'] and m['pdb'][0]=='%':
                    # if there are some pdbs not to check on this uniport
                    not_to_check = m['pdb'].replace("%","").split(",")
                    for pdb in not_to_check:
                        if pdb not in mut_pdb_list:
                            mut_pdb_list[pdb] = []
                    if c_pdb in not_to_check:
                        continue
                # if not_to_check:
                #     print(c_pdb,not_to_check,m['pdb'][0])
                protein = Protein.objects.filter(entry_name=c_pdb.lower()).get()
                # print(c_pdb,protein)
                if c_pdb in cached_mutations:
                    d = cached_mutations[c_pdb]
                else:
                    d = cache.get(c_pdb+"_auto_d")
                    if not d:
                        d = fetch_pdb_info(c_pdb,protein,ignore_gasper_annotation=True)
                        cache.set(c_pdb+"_auto_d",d,60*60*24)
                    cached_mutations[c_pdb] = d
                # Find construct mutation
                cons_muts = ConstructMutation.objects.filter(construct=c, sequence_number = m['pos'], mutated_amino_acid = m['mut_aa'], wild_type_amino_acid = m['wt_aa'])

                if not cons_muts.exists() and m['other_effect']!='Non-receptor' and m['other_effect']!='Wrong annotation - remove!':
                    # If no hits something is odd
                    # print(c.structure.pdb_code.index,' do not have following mutation:',mut)
                    found = False
                    for pdb_m in d['mutations']:
                        if int(pdb_m['pos']) == m['pos']:
                            found = True
                            break
                    if found:
                        # print('It was however found in pdb! ADDING')
                        res_wt = Residue.objects.get(protein_conformation__protein=protein.parent, sequence_number=m['pos'])
                        mut = ConstructMutation.objects.create(construct=c, sequence_number=m['pos'],wild_type_amino_acid= m['wt_aa'],mutated_amino_acid=m['mut_aa'], residue=res_wt)
                        pdbs_has.append(c_pdb)
                    else:
                        # print('Was also not found in pdb!')
                        pdbs_hasnot.append("%"+c_pdb)
                        cons_muts_odd = ConstructMutation.objects.filter(construct=c, sequence_number = m['pos'])
                        for cons_mut in cons_muts_odd:
                            print(c_pdb,cons_mut)
                else:
                    # print(c.structure.pdb_code.index,' HAS following mutation:',mut)
                    pdbs_has.append(c_pdb)

                cons_muts = ConstructMutation.objects.filter(construct=c, sequence_number = m['pos'], mutated_amino_acid = m['mut_aa'], wild_type_amino_acid = m['wt_aa'])
                for cons_mut in cons_muts:
                    if m['other_effect']=='Non-receptor' or m['other_effect']=='Wrong annotation - remove!':
                        # print('Delete!',cons_mut.construct.structure.pdb_code.index,cons_mut)
                        cons_mut.delete()
                        continue
                    # Clear existing to replace with current
                    cons_mut.effects.clear()

                    if m['thermo_effect']:
                        mutation_type, created = ConstructMutationType.objects.get_or_create(slug=slugify('Thermostabilising'),name='Thermostabilising', effect=m['thermo_effect'])
                        cons_mut.effects.add(mutation_type)

                    if m['expression_effect']:
                        mutation_type, created = ConstructMutationType.objects.get_or_create(slug=slugify('Receptor Expression'),name='Receptor Expression', effect=m['expression_effect'])
                        cons_mut.effects.add(mutation_type)

                    if m['site_effect']:
                        mutation_type, created = ConstructMutationType.objects.get_or_create(slug=slugify(m['site_effect']),name=m['site_effect'], effect=m['site_effect_type'])
                        cons_mut.effects.add(mutation_type)

                    if m['other_effect']:
                        # print(m['other_effect'])
                        mutation_type, created = ConstructMutationType.objects.get_or_create(slug=slugify('Other effect'),name='Other effect', effect=m['other_effect'])
                        cons_mut.effects.add(mutation_type)

                    track_annotated_mutations.append(cons_mut.pk)
            # if not m['pdb'] and len(pdbs_hasnot) and len(pdbs_has):
            #     print(m['entry_name'],m['wt_aa']+str(m['pos'])+m['mut_aa'])
            #     print("has",",".join(pdbs_has))
            #     print("hasnot",",".join(pdbs_hasnot))
            # if not len(pdbs_has):
            #     print('NOONE HAS',m['entry_name'],m['wt_aa']+str(m['pos'])+m['mut_aa'])
        print(len(track_annotated_mutations),'annotated mutations')

        non_annotated_muts = ConstructMutation.objects.all().exclude(pk__in=track_annotated_mutations).order_by('construct__protein__entry_name','sequence_number')
        print(len(non_annotated_muts),'non-annotated mutations')
        csv_rows = [['reference','pdb','construct name','class','lig type','rec fam','uniprot','segment','gpcrdb#','AA no.','WT aa','Mut aa','','','','','','Remark']]


        missing = list(set(self.all_pdbs) - set(mut_pdb_list.keys()))
        print(sorted(missing)," do not have any mutations annotated -- add them to sheet with NONE if they have none")
        missing2 = []
        for c_pdb in missing:
            d = cache.get(c_pdb+"_auto_d")
            if not d:
                protein = Protein.objects.filter(entry_name=c_pdb.lower()).get()
                d = fetch_pdb_info(c_pdb,protein,ignore_gasper_annotation=True)
                cache.set(c_pdb+"_auto_d",d,60*60*24)
            if len(d['mutations']):
                missing2.append(c_pdb)
        print(sorted(missing2)," do not have any mutations annotated (but auto has them having) -- add them to sheet with NONE if they have none")

        for mut in non_annotated_muts:
            pdb = mut.construct.structure.pdb_code.index
            uniprot = mut.construct.protein.entry_name
            seg = mut.residue.protein_segment.slug
            if mut.residue.generic_number:
                gn = mut.residue.generic_number.label
            else:
                gn = ''
            pos = mut.sequence_number
            wt_aa = mut.wild_type_amino_acid
            mut_aa = mut.mutated_amino_acid

            annotated_effect = [e.slug for e in mut.effects.all()]

            # print(annotated_effect)

            csv_rows.append(['',pdb,'','','','',uniprot,seg,gn,pos,wt_aa,mut_aa,'','','','','',','.join(annotated_effect)])
            # print(csv_rows[-1])

       #  print(csv_rows)
        import csv
        with open('construct_mut_missing.csv', 'w') as f:
            writer = csv.writer(f, delimiter = '\t')
            writer.writerows(csv_rows)
