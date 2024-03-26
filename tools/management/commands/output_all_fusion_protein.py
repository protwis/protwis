from django.core.management.base import BaseCommand, CommandError
from django.core.management import call_command
from django.conf import settings
from django.db import connection
from django.db.models import Q
from django.template.loader import render_to_string
from structure.models import Structure
from protein.models import Protein
from residue.models import ResidueGenericNumber, ResidueGenericNumberEquivalent
from construct.functions import *
from common import definitions
from common.selection import SelectionItem
from common.alignment_gpcr import Alignment
import xlsxwriter, xlrd


import logging, json, os



class Command(BaseCommand):

    help = "Output all uniprot mappings"

    logger = logging.getLogger(__name__)


    def handle(self, *args, **options):

        structures = Structure.objects.all().exclude(structure_type__slug__startswith='af-')
        self.fusions = {}
        self.fusions_mutations = {}
        self.fusions_starts = {}
        self.fusions_ends = {}
        f = open('fusions.json', 'w')
        i = 0
        for s in structures:
            pdbname = str(s)
            print(pdbname)
            failed = []
            # try:
            protein = Protein.objects.filter(entry_name=pdbname.lower()).get()
            protein_name = protein.parent.entry_name
            d = fetch_pdb_info(pdbname,protein)
            # except:
            #     print(pdbname,'failed')
            #     failed.append(pdbname)
            # print(d)
            # print(d['construct_sequences'])
            self.find_fusion(d['construct_sequences'],protein_name,pdbname)
            #print(self.fusions)
            i += 1
            # if i>10:
            #     break

        self.save_excel()

    def find_fusion(self,construct_sequences,protein_name,pdbname):
        list_of_none_fusion = ['Not_Observed','Engineered mutation','','Conflict','Insertion','S-arrestin'] # 'Expression tag','Linker',
        list_of_comfirmed_fusion = ['C8TP59','Q0SXH8','Q9V2J8','Soluble cytochrome b562','Endolysin','Rubredoxin','Lysozyme','Flavodoxin']
        #ODD Rubredoxin
        #Q9V2J8 GlgA glycogen synthase  auto_4ZJ8
        #C8TP59  Cytochrome b562
        # Q0SXH8 Cytochrome b(562)
        result = []
        position = None
        for name, segment in construct_sequences.items():
            # print(segment)
            position = segment['where']
            seq = segment['seq']
            mutations = None
            if 'mutations' in segment:
                mutations = segment['mutations']
            wt_ranges = segment['ranges']
            if name in list_of_none_fusion:
                continue
            if name == protein_name:
                continue
            if "Guanine" in name:
                continue
            if "N/A" in name:
                continue

            if name not in self.fusions:
                self.fusions[name] = {}
                self.fusions_mutations[name] = []
                self.fusions_starts[name] = []
                self.fusions_ends[name] = []

            mutation_list = []
            if mutations:
                for mut,v in mutations.items():
                    mutation_list.append([v[0],str(mut),v[1]])
                    if mut not in self.fusions_mutations[name]:
                        self.fusions_mutations[name].append(mut)
            wt_range_str = ""
            for wt_range in wt_ranges:
                wt_range_str += str(wt_range[0])+"-"+str(wt_range[1])
                if wt_range[0] not in self.fusions_starts[name]:
                    self.fusions_starts[name].append(wt_range[0])
                if wt_range[1] not in self.fusions_ends[name]:
                    self.fusions_ends[name].append(wt_range[1])


            self.fusions[name][protein_name+"_"+pdbname] = {'mutations':mutation_list, 'wt_ranges':wt_ranges, 'position':position, 'sequence' : seq}

    def save_excel(self):
        """Convert fusions info to excel file"""

        workbook = xlsxwriter.Workbook('fusions.xlsx')

        for name in self.fusions:

            worksheet_name = name[0:30]
            worksheet_name = worksheet_name.replace("/", " ")

            worksheet = workbook.add_worksheet(worksheet_name)

            headers = ['name','position','seq']
            headers_lookup = {}
            headers_start_end = {}
            header_i = 3
            for i,start in enumerate(sorted(self.fusions_starts[name])):
                headers_start_end[start] = "start_"+str(start)
            print(headers_start_end)
            for i,end in enumerate(sorted(self.fusions_ends[name])):
                if end in headers_start_end:
                    headers_start_end[end] = "both"
                else:
                    headers_start_end[end] = "end_"+str(end)
            print(headers_start_end)

            for end in sorted(headers_start_end):
                if headers_start_end[end]=="both":
                    headers.append("end_"+str(end))
                    headers.append("start_"+str(end))
                    headers_lookup["end_"+str(end)] = header_i
                    header_i +=1
                    headers_lookup["start_"+str(end)] = header_i
                else:
                    headers.append(headers_start_end[end])
                    headers_lookup[headers_start_end[end]] = header_i
                header_i +=1

            print(headers_lookup)

            for i,mut in enumerate(sorted(self.fusions_mutations[name])):
                headers.append("Mut in "+str(mut))
                headers_lookup["mut_"+str(mut)] = header_i
                header_i +=1

            row = 1
            index = {}
            col = 0
            for h in headers:
                worksheet.write(0, col, h)
                index[h] = col
                col += 1
            for xtal_name,v in self.fusions[name].items():
                worksheet.write(row, 0, xtal_name)
                worksheet.write(row, 1, v['position'])
                worksheet.write(row, 2, v['sequence'])
                for mut in v['mutations']:
                    worksheet.write(row, headers_lookup['mut_'+mut[1]], mut[0]+mut[1]+mut[2])
                for r in v['wt_ranges']:
                    worksheet.write(row, headers_lookup['start_'+str(r[0])], "X")
                    worksheet.write(row, headers_lookup['end_'+str(r[1])], "X")
                row += 1

        workbook.close()
