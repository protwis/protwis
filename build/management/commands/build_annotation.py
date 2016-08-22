from django.core.management.base import BaseCommand, CommandError
from django.conf import settings
from django.db import connection
from django.utils.text import slugify


from residue.models import Residue

import logging
import os
import sys
import re
from datetime import datetime
from collections import OrderedDict
from urllib.request import urlopen, quote
import xlrd
import operator
import traceback

class Command(BaseCommand):
    help = 'Reads source data and creates annotations'
    
    def add_arguments(self, parser):
        parser.add_argument('-f', '--filename',
            action='append',
            dest='filename',
            help='Filename to import. Can be used multiple times')
        parser.add_argument('-u', '--purge',
            action='store_true',
            dest='purge',
            default=False,
            help='Purge existing mutations records')

    logger = logging.getLogger(__name__)

    # source file directory
    annotation_source_file = os.sep.join([settings.DATA_DIR, 'structure_data', 'Structural_Annotation.xlsx'])

    def handle(self, *args, **options):
        # delete any existing structure data
        if options['purge']:
            try:
                # self.purge_mutants()
                pass
            except Exception as msg:
                print(msg)
                self.logger.error(msg)
        # import the structure data
        try:
            self.create_annotations(options['filename'])
        except Exception as msg:
            print(msg)
            traceback.print_exc()
            self.logger.error(msg)

    def parse_excel(self,path):
        workbook = xlrd.open_workbook(path)
        worksheets = workbook.sheet_names()
        d = {}
        for worksheet_name in worksheets:
            if worksheet_name in d:
                print('Error, worksheet with this name already loaded')
                continue

            d[worksheet_name] = {}  
            worksheet = workbook.sheet_by_name(worksheet_name)

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
                # print(h)
                headers.append(worksheet.cell_value(0, i))

            for curr_row in range(1,num_rows+1):
                row = worksheet.row(curr_row)
                key = worksheet.cell_value(curr_row, 0)

                if key=='':
                    #in case there is no key for whatever reason
                    # print("no key!")
                    continue

                d[worksheet_name][key] = {}
                temprow = {}
                for curr_cell in range(num_cells):
                    # cell_type = worksheet.cell_type(curr_row, curr_cell)
                    cell_value = worksheet.cell_value(curr_row, curr_cell)
                    # temprow.append(cell_value)
                    d[worksheet_name][key][headers[curr_cell]] = cell_value

                # if curr_row>2: break
        return d


    def generate_bw(self, i, v, aa):
        #return dict
        a = {'aa':aa, 'pos':i, 's':'', 'bw':''}
        if i<int(v['1b']):
            a['s'] = 'N-term'
        elif i<=int(v['1e']):
            a['s'] = 'TM1'
            a['bw'] = '1.'+str(50+i-int(v['1x']))
        elif v['i1x']!="-" and v['i1e']!="-" and i<int(v['2b']):
            if i<int(v['i1b']):
                a['s'] = 'ICL1'
            elif i<=int(v['i1e']):
                a['s'] = 'ICL1'
                a['bw'] = '12.'+str(50+i-int(v['i1x']))
            else:
                a['s'] = 'ICL1'
        elif i<int(v['2b']):
            a['s'] = 'ICL1'
        elif i<=int(v['2e']):
            a['s'] = 'TM2'
            a['bw'] = '2.'+str(50+i-int(v['2x']))
        elif v['e1x']!="-" and i<int(v['3b']):
            if i<int(v['e1b']):
                a['s'] = 'ECL1'
            elif i<=int(v['e1e']):
                a['s'] = 'ECL1'
                a['bw'] = '23.'+str(50+i-int(v['e1x']))
            else:
                a['s'] = 'ECL1'
        elif i<int(v['3b']):
            a['s'] = 'ECL1'
        elif i<=int(v['3e']):
            a['s'] = 'TM3'
            a['bw'] = '3.'+str(50+i-int(v['3x']))
        elif v['i2x']!="-" and i<int(v['4b']):
            if i<int(v['i2b']):
                a['s'] = 'ICL2'
            elif i<=int(v['i2e']):
                a['s'] = 'ICL2'
                a['bw'] = '34.'+str(50+i-int(v['i2x']))
            else:
                a['s'] = 'ICL2'
        elif i<int(v['4b']):
            a['s'] = 'ICL2'
        elif i<=int(v['4e']):
            a['s'] = 'TM4'
            a['bw'] = '4.'+str(50+i-int(v['4x']))
        elif v['e2x']!="-" and i<int(v['5b']) and v['e2b']!="-":
            if i<int(v['e2b']):
                a['s'] = 'ECL2'
            elif i<=int(v['e2e']):
                a['s'] = 'ECL2'
                a['bw'] = '45.'+str(50+i-int(v['e2x']))
            else:
                a['s'] = 'ECL2'
        elif i<int(v['5b']):
            a['s'] = 'ECL2'
        elif i<=int(v['5e']):
            a['s'] = 'TM5'
            a['bw'] = '5.'+str(50+i-int(v['5x']))
        elif i<int(v['6b']):
            a['s'] = 'ICL3'
        elif i<=int(v['6e']):
            a['s'] = 'TM6'
            a['bw'] = '6.'+str(50+i-int(v['6x']))
        elif v['7b']=="-": #fix for npy6r_human
            a['s'] = 'C-term'
        elif i<int(v['7b']):
            a['s'] = 'ECL3'
        elif i<=int(v['7e']):
            a['s'] = 'TM7'
            a['bw'] = '7.'+str(50+i-int(v['7x']))
        elif v['8x']!="-":
            if i<=int(v['8b']):
                a['s'] = 'ICL4'
            elif i<=int(v['8e']):
                a['s'] = 'H8'
                a['bw'] = '8.'+str(50+i-int(v['8x']))
            else:
                a['s'] = 'C-term'
        else:
            a['s'] = 'C-term'

        return a


    def create_annotations(self, data):
        self.logger.info('CREATING ANNOTATIONS')

        data = self.parse_excel(self.annotation_source_file)


        track_performance = []
        track_performance_gn = []

        # print(data)
        for p, v in data["NonXtal_SegEnds_Prot#"].items():
            s = data["Seqs"][p]['Sequence']
            b_and_c = {}
            for entry,gn in data["NonXtal_Bulges_Constr_GPCRdb#"][p].items():
                if len(entry)<3:
                    continue
                if entry[1]=='x' or entry[2]=='x':
                    if gn!='' and gn!="":
                        seg, number = entry.split("x")
                        if seg not in b_and_c:
                            b_and_c[seg] = []
                        b_and_c[seg].append(number)

            # print(p, b_and_c)

            rs = Residue.objects.filter(protein_conformation__protein__entry_name=p, alternative_generic_numbers__scheme__slug='bw').values_list('sequence_number','alternative_generic_numbers__label')
            check = {}
            for r in rs:
                check[r[0]] = r[1]

            rs = Residue.objects.filter(protein_conformation__protein__entry_name=p).values_list('sequence_number','generic_number__label')
            check_gn = {}
            for r in rs:
                check_gn[r[0]] = r[1]

            if not len(check):
                print(p, 'no residues for protein in DB')
                continue


            # try:
                #annotated list
            al = []
            counter = {'No bw in D': 0, 'Different bw': 0, 'No bw in DB': 0, 'bw match': 0, 'Match of no bw': 0, "Different gn": 0}

            for i,aa in enumerate(s):
                i += 1 #fix, since first postion is 1
                #print(i,aa)

                a = self.generate_bw(i,v,aa)

                ##perform bulges / constriction check! 
                ## only do this on bw numbers
                if a['bw']:
                    seg, number = a['bw'].split(".")

                    offset = 0
                    bulge = False
                    if seg in b_and_c:
                        for bc in b_and_c[seg]:
                            if len(bc)>2: #bulge
                                # print(bc[0:2],number)
                                if int(bc[0:2])<50 and int(number)<int(bc[0:2]): #before x50 and before bulge, do smt
                                    offset += 1 #eg if 5x461, then 5.46 becomes 5x461, 5.45 becomes 5x46
                                elif int(bc[0:2])<50 and int(number)==int(bc[0:2]): #before x50 and is bulge, do smt
                                    bulge = True # eg if 5x461, then 5.46 becomes 5x461
                                    gn = str(int(number))+"1"
                                elif int(bc[0:2])>50 and int(number)>int(bc[0:2])+1: #after x50 and after bulge, do smt
                                    offset -= 1 #eg if 2x551, then 2.56 becomes 2x551, 5.57 becomes 5x56
                                elif int(bc[0:2])>50 and int(number)==int(bc[0:2])+1: #after x50 and 1 after bulge, do smt
                                    bulge = True # eg if 2x551, then 2.56 becomes 2x551
                                    gn = str(int(number)-1)+"1"

                            else: #2 numbers, it's a constriction
                                if int(bc[0:2])<50 and int(number)<=int(bc[0:2]): #before x50 and before or equal constrictions, do smt
                                    offset -= 1 #eg if constriction is 7x44, then 7.44 becomes 7x43, 7.43 becomes 7x42
                                if int(bc[0:2])>50 and int(number)>=int(bc[0:2]): #before x50 and before or equal constrictions, do smt
                                    offset += 1 #eg if constriction is 4x57, then 4.57 becomes 4x58, 4.58 becomes 4x59

                        if bulge!=True:
                            gn = str(int(number)+offset)

                        if i in check_gn:
                            if seg+"x"+gn!=check_gn[i]:
                                if i in check:
                                    if a['bw']==check[i]:
                                        counter["Different gn"] += 1

                        # print(p,a['bw']," something in bc",b_and_c[seg])
                
                if i in check:
                    if a['bw']!=check[i]:
                        # print("Change!",i, a['bw']," db:",check[i])
                        if a['bw']:
                            counter['Different bw'] += 1
                        else:
                            counter['No bw in D'] += 1
                    elif a['bw']==check[i]:
                        counter['bw match'] += 1
                elif len(check) and a['bw']:
                    # print("Change!",i, a['bw'],"db: None")
                    counter['No bw in DB'] += 1
                else:
                    counter['Match of no bw'] += 1

                al.append(a)

                    # print(a)
            # except Exception as e:
            #     exc_type, exc_obj, exc_tb = sys.exc_info()
            #     fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
            #     print(exc_type, fname, exc_tb.tb_lineno)

            #     # print('error',p,e)
            #     # pass
            #     # print(p,"pos",i+1,al[-1]['s'],"some kind of error",e)

            percentage = (counter['bw match']+counter['Match of no bw'])*100/len(s)
            percentage_gn = (counter['Different gn'])*100/len(s)
            # print(p, counter,percentage)
            track_performance.append(percentage)
            track_performance_gn.append(percentage_gn)
        # print(track_performance)
        average_performance = sum(track_performance) / float(len(track_performance))
        average_performance_gn = sum(track_performance_gn) / float(len(track_performance_gn))
        print(average_performance,average_performance_gn)
            # break


        self.logger.info('COMPLETED ANNOTATIONS')