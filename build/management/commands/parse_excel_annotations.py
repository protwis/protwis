from django.core.management.base import BaseCommand, CommandError
from django.core.management import call_command
from django.conf import settings
from django.db import connection

from common.alignment import Alignment, ClosestReceptorHomolog
from protein.models import Protein, ProteinSegment
from structure.models import Structure

import datetime
import logging
from optparse import make_option
import os
import shutil
import xlrd
import yaml
from collections import OrderedDict
import pprint


_mapping_tag = yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG
def dict_constructor(loader, node):
    return OrderedDict(loader.construct_pairs(node))

def represent_ordereddict(dumper, data):
    value = []

    for item_key, item_value in data.items():
        node_key = dumper.represent_data(item_key)
        node_value = dumper.represent_data(item_value)

        value.append((node_key, node_value))

    return yaml.nodes.MappingNode(u'tag:yaml.org,2002:map', value)

yaml.add_representer(OrderedDict, represent_ordereddict)
yaml.add_constructor(_mapping_tag, dict_constructor)


class Command(BaseCommand):
    help = 'Basic functions for build scrips'

    logger = logging.getLogger(__name__)

    def add_arguments(self, parser):
        parser.add_argument('-f', '--filename',
            action='store',
            dest='filename',
            help='Path to Uniprot text file')
        parser.add_argument('-m',
            action='store_true',
            default=False,
            help='Run main template search. Updates Xtal_Templ.csv with closest receptor homologs')

    annotation_source_file = os.sep.join([settings.DATA_DIR, 'structure_data', 'Structural_Annotation.xlsx'])
    xtal_seg_end_file = os.sep.join([settings.DATA_DIR, 'structure_data', 'annotation', 'xtal_segends.yaml'])
    mod_xtal_seg_end_file = os.sep.join([settings.DATA_DIR, 'structure_data', 'annotation', 'mod_xtal_segends.yaml'])
    xtal_seg_end_bw_file = os.sep.join([settings.DATA_DIR, 'structure_data', 'annotation', 'xtal_segends_bw.yaml'])
    ECD_annotation_source_file = os.sep.join([settings.DATA_DIR, 'structure_data', 'ECD_annotation.xlsx'])
    ClassD_annotation_source_file = os.sep.join([settings.DATA_DIR, 'structure_data', 'Class_D_Annotation.xlsx'])
    GPCRdb_structure_info = os.sep.join([settings.DATA_DIR, 'structure_data', 'annotation', 'GPCRdb_structure_info.xlsx'])

    non_xtal_seg_end_file = os.sep.join([settings.DATA_DIR, 'structure_data', 'annotation', 'non_xtal_segends.yaml'])
    non_xtal_seg_end_bw_file = os.sep.join([settings.DATA_DIR, 'structure_data', 'annotation', 'non_xtal_segends_bw.yaml'])

    all_anomalities_file = os.sep.join([settings.DATA_DIR, 'structure_data', 'annotation', 'all_anomalities.yaml'])
    xtal_anomalities_file = os.sep.join([settings.DATA_DIR, 'structure_data', 'annotation', 'xtal_anomalities.yaml'])

    sequence_file = os.sep.join([settings.DATA_DIR, 'structure_data', 'annotation', 'sequences.yaml'])

    ECD_wt_file = os.sep.join([settings.DATA_DIR, 'structure_data', 'annotation', 'ECD_wt.yaml'])
    ECD_anomalies_file = os.sep.join([settings.DATA_DIR, 'structure_data', 'annotation', 'ECD_anomalies.yaml'])


    if not os.path.exists(os.sep.join([settings.DATA_DIR, 'structure_data', 'annotation'])):
        os.makedirs(os.sep.join([settings.DATA_DIR, 'structure_data', 'annotation']))

    def handle(self, *args, **options):
        self.data = self.parse_excel(self.annotation_source_file)
        # self.dump_files()
        self.ECD_data = self.parse_excel(self.ECD_annotation_source_file)
        self.dump_ECD_files()
        self.ClassD_data = self.parse_excel(self.ClassD_annotation_source_file)
        # self.dump_ClassD_data()
        self.GPCRdb_structure_data = self.parse_excel(self.GPCRdb_structure_info)
        self.dump_GPCRdb_structure_data()
        # self.analyse_annotation_consistency()
        self.find_representatives()
        if options['m']:
            self.main_template_search()

    def dict_to_csv(self, file_name, tab_name):
        with open(os.sep.join([settings.DATA_DIR, 'structure_data', 'annotation', file_name+'.csv']), 'w') as f1:
            c = 0
            for key, val in self.GPCRdb_structure_data[tab_name].items():
                if c==0:
                    header = '\t'.join([i for i in val])+'\n'
                    f1.write(header)
                line = '\t'.join([str(j) for i,j in val.items()])+'\n'
                f1.write(line)
                c+=1

    def dump_GPCRdb_structure_data(self):
        structures, ligands, nanobodies, fusion_proteins, g_proteins, arrestins, ramps, grks = OrderedDict(),OrderedDict(),OrderedDict(),OrderedDict(),OrderedDict(),OrderedDict(),OrderedDict(),OrderedDict()
        for i,j in [('structures', 'Structures'), ('ligands','Ligands'), ('nanobodies','Nanobodies'), ('fusion_proteins','Fusion proteins'), ('g_proteins','G protein'), ('arrestins','Arrestin'), ('ramp','RAMP'), ('grk', 'GRK')]:
            self.dict_to_csv(i,j)

    def dump_ECD_files(self):
        data_dict = OrderedDict()
        for key, val in self.ECD_data['wt'].items():
            if val['H1x50']=='':
                continue
            entry_name = val['UniProt']
            del val['Key']
            del val['UniProt']
            data_dict[entry_name] = val
        with open(self.ECD_wt_file, 'w') as outfile:
            yaml.dump(data_dict, outfile, indent=4)
        anomalies = OrderedDict()
        for key, val in self.ECD_data['anomalies'].items():
            entry_name = val['protein']
            del val['protein']
            anomalies[entry_name] = val
        with open(self.ECD_anomalies_file, 'w') as outfile:
            yaml.dump(anomalies, outfile, indent=4)


    def dump_ClassD_data(self):
        data_dict1, data_dict2 = OrderedDict(), OrderedDict()
        for key, val in self.ClassD_data['SegEnds_NonXtal_Prot#'].items():
            entry_name = val['UniProt'].lower()
            del val['Key']
            del val['UniProt']
            del val['']
            data_dict1[entry_name] = val
        with open(self.non_xtal_seg_end_file, 'a') as outfile:
            yaml.dump(data_dict1, outfile, indent=4)
        for key, val in self.ClassD_data['SegEnds_NonXtal_BW#'].items():
            entry_name = val['UniProt'].lower()
            del val['UniProt']
            data_dict2[entry_name] = val
        with open(self.non_xtal_seg_end_bw_file, 'a') as outfile:
            yaml.dump(data_dict2, outfile, indent=4)
        data = self.ClassD_data["Bulges_Constrictions"]
        NonXtal_Bulges_Constr_GPCRdb = {}
        for structure,vals in data.items():
            entry = structure.lower()
            NonXtal_Bulges_Constr_GPCRdb[entry] = OrderedDict()
            for key,val in vals.items():
                if not key:
                    continue
                NonXtal_Bulges_Constr_GPCRdb[entry][key] = val
        NonXtal_Bulges_Constr_GPCRdb = OrderedDict(sorted(NonXtal_Bulges_Constr_GPCRdb.items()))
        with open(self.all_anomalities_file, 'a') as outfile:
            yaml.dump(NonXtal_Bulges_Constr_GPCRdb, outfile, indent=4)
        data = self.ClassD_data["Seqs"]
        Seqs = {}
        for structure,vals in data.items():
            entry = structure.lower()
            Seqs[entry] = OrderedDict()
            for key,val in vals.items():
                if not key:
                    continue
                Seqs[entry][key] = val
        Seqs = OrderedDict(sorted(Seqs.items()))
        with open(self.sequence_file, 'a') as outfile:
            yaml.dump(Seqs, outfile, indent=4)

        structures = self.ClassD_data["SegEnds_Xtal_Prot#"]
        pdb_info = {}
        pdb_info_all = {}
        for structure,vals in structures.items():
            if structure.split("_")[-1] == "wt":
                continue
            if structure.split("_")[-1] == "dist":
                continue
            #print(structure)
            pdb_id = structure.split("_")[-1]
            pdb_info[pdb_id] = OrderedDict()
            for key,val in vals.items():
                if len(key)>3:
                    continue
                if not key:
                    continue
                if key[-1]!="b" and key[-1]!="e":
                    continue
                pdb_info[pdb_id][key] = val

        for structure,vals in structures.items():
            entry = structure
            pdb_info_all[entry] = OrderedDict()
            for key,val in vals.items():
                if len(key)>3:
                    continue
                if not key:
                    continue
                if key[-1]!="b" and key[-1]!="e":
                    continue
                pdb_info_all[entry][key] = val

        pdb_info = OrderedDict(sorted(pdb_info.items()))
        with open(self.mod_xtal_seg_end_file, 'a') as outfile:
            yaml.dump(pdb_info, outfile, indent=4)
        pdb_info_all = OrderedDict(sorted(pdb_info_all.items()))
        with open(self.xtal_seg_end_file, 'a') as outfile:
            yaml.dump(pdb_info_all, outfile, indent=4)

    def parse_excel(self,path,remove_header_linebreak=False):
        self.workbook = xlrd.open_workbook(path)
        worksheets = self.workbook.sheet_names()
        d = {}
        for worksheet_name in worksheets:
            if worksheet_name in d:
                print('Error, worksheet with this name already loaded')
                continue

            d[worksheet_name] = OrderedDict()
            worksheet = self.workbook.sheet_by_name(worksheet_name)

            num_rows = worksheet.nrows - 1
            num_cells = worksheet.ncols - 1
            curr_row = 0 #skip first, otherwise -1

            headers = []
            for i in range(num_cells):
                h = worksheet.cell_value(0, i)
                if h=="":
                    #replace header with index if empty
                    h = "i_"+str(i)
                if h in headers:
                    # print('already have ',h)
                    h += "_"+str(i)
                if remove_header_linebreak:
                    h = h.replace('\n', ' ')
                    headers.append(h)
                else:
                    headers.append(worksheet.cell_value(0, i))
            for curr_row in range(1,num_rows+1):
                row = worksheet.row(curr_row)
                key = worksheet.cell_value(curr_row, 0)
                # Workaround for structures as there can be multiple entries per structure
                # TODO - discuss if this shouldn't be the default?
                if "GPCRdb_structure_info.xlsx" in path:
                    key = curr_row


                if key=='':
                    #in case there is no key for whatever reason
                    # print("no key!")
                    continue
                # if key in d[worksheet_name]:
                #     print(key, "already in",worksheet_name)
                d[worksheet_name][key] = OrderedDict()
                temprow = {}
                for curr_cell in range(num_cells):
                    # cell_type = worksheet.cell_type(curr_row, curr_cell)
                    cell_value = worksheet.cell_value(curr_row, curr_cell)
                    if "GPCRdb_structure_info.xlsx" in path and headers[curr_cell]=='Date':
                        cell_value = str(datetime.datetime(*xlrd.xldate_as_tuple(cell_value, self.workbook.datemode)))[:10]

                    # temprow.append(cell_value)
                    if headers[curr_cell] not in d[worksheet_name][key]:
                        #do not overwrite
                        d[worksheet_name][key][headers[curr_cell]] = cell_value
                # if curr_row>2: break
        return d

    def analyse_annotation_consistency(self):
        NonXtal = self.data["Bulges_Constr_NonXtal_GPCRdb#"]
        Xtal = self.data["Bulges_Constr_Xtal_GPCRdb#"]
        output = {}
        counting_xtal = {}
        counting_non_xtal = {}
        for entry_protein,vals in NonXtal.items():
            anomalies=[]
            anomalies_xtal=[]
            for key,val in vals.items():
                if "x" in val and "_" not in val:
                    if val.index("x") in [1,2]:
                        anomalies.append(val)
            if vals['Xtal Templ'] in Xtal:
                #print(Xtal[vals['Xtal Templ']])
                for key,val in Xtal[vals['Xtal Templ']].items():
                    if "x" in val and "_" not in val:
                        if val.index("x") in [1,2]:
                            anomalies_xtal.append(val)

            if entry_protein==vals['Xtal Templ']:
                list1 = list(set(anomalies) - set(anomalies_xtal))
                list2 = list(set(anomalies_xtal) - set(anomalies))
                if list1 or list2:
                    for i in list1:
                        if i not in counting_non_xtal:
                            counting_non_xtal[i] = 0
                        counting_non_xtal[i] += 1
                    for i in list2:
                        if i not in counting_xtal:
                            counting_xtal[i] = 0
                        counting_xtal[i] += 1
                    #print("ISSUE!")
                    #print(entry_protein)
                    #print("NonXtal_anomalies",anomalies,"Xtal_anomalies",anomalies_xtal)
                    if list1: print(entry_protein,vals['Xtal Templ'],"Present in non-xtal, but not xtal",list1)
                    if list2: print(entry_protein,vals['Xtal Templ'],"Present in xtal, but not non-xtal",list2)

        print("Overall")
        print("Present in non-xtal, but not xtal",counting_xtal)
        print("Present in xtal, but not non-xtal",counting_non_xtal)

        structures = self.data["SegEnds_Xtal_Prot#"]
        structures_non_xtal = self.data["SegEnds_NonXtal_Prot#"]
        info = {}
        for structure,vals in structures.items():
            if structure.split("_")[-1] == "wt":
                # print(structure)
                entry = vals['UniProt']
                info[entry] = {}
                for key,val in vals.items():
                    # print(val,key)
                    if len(key)>3:
                        continue
                    if not key:
                        continue
                    if key[-1]!="b" and key[-1]!="e":
                        continue

                    info[entry][key] = val

                    if structures_non_xtal[entry][key]!=val:
                        print("error with ",entry,key,"Xtal sheet:",val,"NonXtal sheet:",structures_non_xtal[entry][key])
                        print(structures_non_xtal[entry])
                        print(vals)
                #print(structure,info)

        # with open(self.xtal_seg_end_file, 'w') as outfile:
        #     yaml.dump(pdb_info, outfile)

    def dump_files(self):
        structures = self.data["SegEnds_Xtal_Prot#"]
        pdb_info = {}
        pdb_info_all = {}
        for structure,vals in structures.items():
            if structure.split("_")[-1] == "wt":
                continue
            if structure.split("_")[-1] == "dist":
                continue
            #print(structure)
            pdb_id = structure.split("_")[-1]
            pdb_info[pdb_id] = OrderedDict()
            for key,val in vals.items():
                if len(key)>3:
                    continue
                if not key:
                    continue
                if key[-1]!="b" and key[-1]!="e":
                    continue
                pdb_info[pdb_id][key] = val

        for structure,vals in structures.items():
            entry = structure
            pdb_info_all[entry] = OrderedDict()
            for key,val in vals.items():
                if len(key)>3:
                    continue
                if not key:
                    continue
                if key[-1]!="b" and key[-1]!="e":
                    continue
                pdb_info_all[entry][key] = val

        data = self.data["SegEnds_Xtal_BW#"]
        Xtal_SegEnds_BW = {}
        for structure,vals in data.items():
            entry = structure
            Xtal_SegEnds_BW[entry] = OrderedDict()
            for key,val in vals.items():
                if not key:
                    continue
                if len(key)>3 and key[-1]!="b" and key[-1]!="e":
                    continue
                Xtal_SegEnds_BW[entry][key] = val

        data = self.data["SegEnds_NonXtal_BW#"]
        NonXtal_SegEnds_BW = {}
        for structure,vals in data.items():
            entry = structure
            NonXtal_SegEnds_BW[entry] = OrderedDict()
            for key,val in vals.items():
                if not key:
                    continue
                if len(key)>3 and key[-1]!="b" and key[-1]!="e" and key!="XtalTempl":
                    continue
                NonXtal_SegEnds_BW[entry][key] = val

        data = self.data["SegEnds_NonXtal_Prot#"]
        NonXtal_SegEnds_Prot = {}
        for structure,vals in data.items():
            entry = structure
            NonXtal_SegEnds_Prot[entry] = OrderedDict()
            for key,val in vals.items():
                if not key:
                    continue
                if len(key)>3 and key[-1]!="b" and key[-1]!="e" and key!="Xtal Templ":
                    continue
                NonXtal_SegEnds_Prot[entry][key] = val


        # data = self.data["Bulges_Constr_Xtal_GPCRdb#"]
        # Xtal_Bulges_Constr_GPCRdb = {}
        # for structure,vals in data.items():
        #     entry = structure
        #     Xtal_Bulges_Constr_GPCRdb[entry] = OrderedDict()
        #     for key,val in vals.items():
        #         if not key:
        #             continue
        #         Xtal_Bulges_Constr_GPCRdb[entry][key] = val

        data = self.data["Bulges_Constr_NonXtal_GPCRdb#"]
        NonXtal_Bulges_Constr_GPCRdb = {}
        for structure,vals in data.items():
            entry = structure
            NonXtal_Bulges_Constr_GPCRdb[entry] = OrderedDict()
            for key,val in vals.items():
                if not key:
                    continue
                NonXtal_Bulges_Constr_GPCRdb[entry][key] = val


        data = self.data["Seqs"]
        Seqs = {}
        for structure,vals in data.items():
            entry = structure
            Seqs[entry] = OrderedDict()
            for key,val in vals.items():
                if not key:
                    continue
                Seqs[entry][key] = val


        pdb_info = OrderedDict(sorted(pdb_info.items()))
        with open(self.mod_xtal_seg_end_file, 'w') as outfile:
            yaml.dump(pdb_info, outfile, default_flow_style=False)

        pdb_info_all = OrderedDict(sorted(pdb_info_all.items()))
        with open(self.xtal_seg_end_file, 'w') as outfile:
            yaml.dump(pdb_info_all, outfile, indent=4)

        Xtal_SegEnds_BW = OrderedDict(sorted(Xtal_SegEnds_BW.items()))
        with open(self.xtal_seg_end_bw_file, 'w') as outfile:
            yaml.dump(Xtal_SegEnds_BW, outfile, indent=4)

        NonXtal_SegEnds_BW = OrderedDict(sorted(NonXtal_SegEnds_BW.items()))
        with open(self.non_xtal_seg_end_bw_file, 'w') as outfile:
            yaml.dump(NonXtal_SegEnds_BW, outfile, indent=4)

        NonXtal_SegEnds_Prot = OrderedDict(sorted(NonXtal_SegEnds_Prot.items()))
        with open(self.non_xtal_seg_end_file, 'w') as outfile:
            yaml.dump(NonXtal_SegEnds_Prot, outfile, indent=4)

        # Xtal_Bulges_Constr_GPCRdb = OrderedDict(sorted(Xtal_Bulges_Constr_GPCRdb.items()))
        # with open(self.xtal_anomalities_file, 'w') as outfile:
        #     yaml.dump(Xtal_Bulges_Constr_GPCRdb, outfile, indent=4)

        NonXtal_Bulges_Constr_GPCRdb = OrderedDict(sorted(NonXtal_Bulges_Constr_GPCRdb.items()))
        with open(self.all_anomalities_file, 'w') as outfile:
            yaml.dump(NonXtal_Bulges_Constr_GPCRdb, outfile, indent=4)

        Seqs = OrderedDict(sorted(Seqs.items()))
        with open(self.sequence_file, 'w') as outfile:
            yaml.dump(Seqs, outfile, indent=4)

    def main_template_search(self):
        output_csv = ''
        changes = {}
        counter = 0
        for protein, values in self.data['Xtal_Templ'].items():
            values = self.data['Xtal_Templ'][protein]
            crh = ClosestReceptorHomolog(protein)
            closest_hom = crh.find_closest_receptor_homolog()
            if values['Template']!=closest_hom.entry_name:
                changes[protein] = [values['Template'], closest_hom.entry_name]
            output_csv+='{},{}\n'.format(protein, closest_hom.entry_name)
            counter+=1
        with open(os.sep.join([settings.DATA_DIR,'structure_data','annotation','xtal_templates.csv']),'w') as f:
            f.write(output_csv)
        if len(changes)>0:
            print('Changed {} entries out of {} (reference: [changed_from, changed_to]):'.format(len(changes), counter))
            print(changes)
            print('INFO: xtal_templates.csv file updated. Please update Structural_Annotation.xlsx Xtal_Templ sheet with this .csv')
        return changes

    def find_representatives(self):
        grouped = {}
        counter = 0
        xtals, nums, states, resolutions = [], [], [], []
        out = OrderedDict()
        exceptions = ['4L6R']
        with open(os.sep.join([settings.DATA_DIR,'structure_data','annotation','xtal_representatives.yaml']), 'w') as outfile:
            for key, values in self.data['SegEnds_Xtal_Prot#'].items():
                if counter==0:
                    prev_rec = values['UniProt']
                counter+=1
                if values['PDB']=='_wt' or 'dist' in key:
                    continue
                if values['Repr']!='-':
                    if values['Repr']=='Repr_Act':
                        actstat = 'Active'
                    elif values['Repr']=='Repr_Inter':
                        actstat = 'Intermediate'
                    elif values['Repr']=='Repr_Inact':
                        actstat = 'Inactive'
                    out[values['PDB']] = actstat
            yaml.dump(out, outfile, default_flow_style=False)
                # if prev_rec!=values['UniProt'] or counter==len(self.data['SegEnds_Xtal_Prot#']):
                #     if counter==len(self.data['SegEnds_Xtal_Prot#']):
                #         xtals.append(key)
                #         nums.append(values['#Res'])
                #         states.append(values['State'])
                #         resolutions.append(values['Resolution'])
                #     if len(xtals)>0:
                #         max_num_ia, max_x_ia, max_num_a, max_x_a, ia_count, a_count = 0, 0, 0, 0, 0, 0
                #         for x, n, s, r in zip(xtals, nums, states, resolutions):
                #             if s=='Inact':
                #                 if ia_count==0:
                #                     max_res_ia = r
                #                 if n>max_num_ia and x[-4:] not in exceptions:
                #                     max_num_ia = n
                #                     max_x_ia = x
                #                     max_res_ia = r
                #                 elif n==max_num_ia and x[-4:] not in exceptions:
                #                     if r<max_res_ia:
                #                         max_num_ia = n
                #                         max_x_ia = x
                #                         max_res_ia = r
                #                 ia_count+=1
                #             elif s=='Act':
                #                 if a_count==0:
                #                     max_res_a = r
                #                 if n>max_num_a and x[-4:] not in exceptions:
                #                     max_num_a = n
                #                     max_x_a = x
                #                 elif n==max_num_a and x[-4:] not in exceptions:
                #                     if r<max_res_a:
                #                         max_num_a = n
                #                         max_x_a = x
                #                         max_res_a = r
                #                 a_count+=1
                #         for x, n in zip(xtals, nums):
                #             if x==max_x_ia:
                #                 out[x] = 'Repr_Inact'
                #             elif x==max_x_a:
                #                 out[x] = 'Repr_Act'
                #             else:
                #                 out[x] = '-'
                #         yaml.dump(out, outfile, indent=4)
                #         xtals, nums, states, resolutions = [], [], [], []
                #         out = OrderedDict()
                #         xtals.append(key)
                #         nums.append(values['#Res'])
                #         states.append(values['State'])
                #         resolutions.append(values['Resolution'])
                # else:
                #     xtals.append(key)
                #     nums.append(values['#Res'])
                #     states.append(values['State'])
                #     resolutions.append(values['Resolution'])
                # prev_rec = values['UniProt']
