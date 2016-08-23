from django.core.management.base import BaseCommand, CommandError
from django.conf import settings
from django.db import connection
from django.utils.text import slugify
from django.db import IntegrityError


from build.management.commands.base_build import Command as BaseBuild
from residue.models import Residue
from residue.functions import *
from protein.models import Protein, ProteinConformation, ProteinSegment, ProteinFamily

import logging
import os
import sys
import re
from datetime import datetime
from collections import OrderedDict
from itertools import islice
from urllib.request import urlopen, quote
import xlrd
import operator
import traceback

class Command(BaseBuild):
    help = 'Reads source data and creates annotations'
    
    logger = logging.getLogger(__name__)

    # source file directory
    annotation_source_file = os.sep.join([settings.DATA_DIR, 'structure_data', 'Structural_Annotation.xlsx'])
    generic_numbers_source_dir = os.sep.join([settings.DATA_DIR, 'residue_data', 'generic_numbers'])

    segments = ProteinSegment.objects.filter(partial=False)
    all_segments = {ps.slug: ps for ps in ProteinSegment.objects.all()}  # all segments dict for faster lookups
    schemes = parse_scheme_tables(generic_numbers_source_dir)


    def handle(self, *args, **options):
        try:
            self.logger.info('CREATING RESIDUES')
            self.data = self.parse_excel(self.annotation_source_file)

            # run the function twice (second run for proteins without reference positions)
            # self.prepare_input(options['proc'], self.data["NonXtal_SegEnds_Prot#"])
            self.prepare_input(4, self.data["NonXtal_SegEnds_Prot#"])

            self.logger.info('COMPLETED CREATING RESIDUES')
        except Exception as msg:
            print(msg)
            self.logger.error(msg)

    def parse_excel(self,path):
        workbook = xlrd.open_workbook(path)
        worksheets = workbook.sheet_names()
        d = {}
        for worksheet_name in worksheets:
            if worksheet_name in d:
                print('Error, worksheet with this name already loaded')
                continue

            d[worksheet_name] = OrderedDict()
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
        a = {'aa':aa, 'pos':i, 's':'', 'numbers':{'bw':''}}
        if i<int(v['1b']):
            a['s'] = 'N-term'
        elif i<=int(v['1e']):
            a['s'] = 'TM1'
            a['numbers']['bw'] = '1.'+str(50+i-int(v['1x']))
        elif v['i1x']!="-" and v['i1e']!="-" and i<int(v['2b']):
            if i<int(v['i1b']):
                a['s'] = 'ICL1'
            elif i<=int(v['i1e']):
                a['s'] = 'ICL1'
                a['numbers']['bw'] = '12.'+str(50+i-int(v['i1x']))
            else:
                a['s'] = 'ICL1'
        elif i<int(v['2b']):
            a['s'] = 'ICL1'
        elif i<=int(v['2e']):
            a['s'] = 'TM2'
            a['numbers']['bw'] = '2.'+str(50+i-int(v['2x']))
        elif v['e1x']!="-" and i<int(v['3b']):
            if i<int(v['e1b']):
                a['s'] = 'ECL1'
            elif i<=int(v['e1e']):
                a['s'] = 'ECL1'
                a['numbers']['bw'] = '23.'+str(50+i-int(v['e1x']))
            else:
                a['s'] = 'ECL1'
        elif i<int(v['3b']):
            a['s'] = 'ECL1'
        elif i<=int(v['3e']):
            a['s'] = 'TM3'
            a['numbers']['bw'] = '3.'+str(50+i-int(v['3x']))
        elif v['i2x']!="-" and i<int(v['4b']):
            if i<int(v['i2b']):
                a['s'] = 'ICL2'
            elif i<=int(v['i2e']):
                a['s'] = 'ICL2'
                a['numbers']['bw'] = '34.'+str(50+i-int(v['i2x']))
            else:
                a['s'] = 'ICL2'
        elif i<int(v['4b']):
            a['s'] = 'ICL2'
        elif i<=int(v['4e']):
            a['s'] = 'TM4'
            a['numbers']['bw'] = '4.'+str(50+i-int(v['4x']))
        elif v['e2x']!="-" and i<int(v['5b']) and v['e2b']!="-":
            if i<int(v['e2b']):
                a['s'] = 'ECL2'
            elif i<=int(v['e2e']):
                a['s'] = 'ECL2'
                a['numbers']['bw'] = '45.'+str(50+i-int(v['e2x']))
            else:
                a['s'] = 'ECL2'
        elif i<int(v['5b']):
            a['s'] = 'ECL2'
        elif i<=int(v['5e']):
            a['s'] = 'TM5'
            a['numbers']['bw'] = '5.'+str(50+i-int(v['5x']))
        elif i<int(v['6b']):
            a['s'] = 'ICL3'
        elif i<=int(v['6e']):
            a['s'] = 'TM6'
            a['numbers']['bw'] = '6.'+str(50+i-int(v['6x']))
        elif v['7b']=="-": #fix for npy6r_human
            a['s'] = 'C-term'
        elif i<int(v['7b']):
            a['s'] = 'ECL3'
        elif i<=int(v['7e']):
            a['s'] = 'TM7'
            a['numbers']['bw'] = '7.'+str(50+i-int(v['7x']))
        elif v['8x']!="-":
            if i<=int(v['8b']):
                a['s'] = 'ICL4'
            elif i<=int(v['8e']):
                a['s'] = 'H8'
                a['numbers']['bw'] = '8.'+str(50+i-int(v['8x']))
            else:
                a['s'] = 'C-term'
        else:
            a['s'] = 'C-term'

        if a['numbers']['bw'] == '':
            a['numbers'].pop('bw', None)

        return a


    def main_func(self, positions, iteration):
        self.logger.info('CREATING ANNOTATIONS')
        # pconfs
        if not positions[1]:
            # proteins = OrderedDict(islice(self.data["NonXtal_SegEnds_Prot#"].items(),positions[0]))
            proteins = list(self.data["NonXtal_SegEnds_Prot#"])[positions[0]:]
        else:
            # proteins = OrderedDict(islice(self.data["NonXtal_SegEnds_Prot#"].items(),positions[0],positions[1]))
            proteins = list(self.data["NonXtal_SegEnds_Prot#"])[positions[0]:positions[1]]

        # print(data)
        counter = 0
        for p in proteins:
            self.logger.info('DOING {}'.format(p))
            v = self.data["NonXtal_SegEnds_Prot#"][p]
            counter += 1
            if p=='':
                continue
            # if counter>20:
            #     break
            s = self.data["Seqs"][p]['Sequence']
            b_and_c = {}
            for entry,gn in self.data["NonXtal_Bulges_Constr_GPCRdb#"][p].items():
                if len(entry)<3:
                    continue
                if entry[1]=='x' or entry[2]=='x':
                    if gn!='' and gn!="":
                        seg, number = entry.split("x")
                        if seg not in b_and_c:
                            b_and_c[seg] = []
                        b_and_c[seg].append(number)

            self.logger.info('Parsed Seq and B&C {}'.format(p))

            if ProteinConformation.objects.filter(protein__entry_name=p).exists():
                pconf = ProteinConformation.objects.filter(protein__entry_name=p).get()

                al = []

                for i,aa in enumerate(s):
                    i += 1 #fix, since first postion is 1
                    res = self.generate_bw(i,v,aa)

                    segment = self.all_segments[res['s']]

                    ##perform bulges / constriction check! 
                    ## only do this on bw numbers
                    if 'bw' in res['numbers']:
                        seg, number = res['numbers']['bw'].split(".")

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

                        res['numbers']['generic_number'] = seg+"x"+gn

                    create_or_update_residue(pconf, segment, self.schemes,res)

                    al.append(res)

        self.logger.info('COMPLETED ANNOTATIONS')