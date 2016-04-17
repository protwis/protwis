from django.core.management.base import BaseCommand, CommandError
from django.conf import settings

from residue.models import Residue
from structure.models import Structure

import os
import yaml
from collections import OrderedDict

class Command(BaseCommand):
    help = 'Reads a structure annotation file exported from Maestro, and adds sequence numbers'

    def add_arguments(self, parser):
        parser.add_argument('-f', '--filename',
            action='store',
            dest='filename',
            help='Path to Uniprot text file')

    segment_dict = {
        'N': 'N-term',
        'i1': 'ICL1',
        'e1': 'ECL1',
        'i2': 'ICL2',
        'e2': 'ECL2',
        'i3': 'ICL3',
        'e3': 'ECL3',
        'TM1': 'TM1',
        'TM2': 'TM2',
        'TM3': 'TM3',
        'TM4': 'TM4',
        'TM5': 'TM5',
        'TM6': 'TM6',
        'TM7': 'TM7',
        'H8': 'H8',
    }
    numbering_prefix_dict = {
        'i1': '12',
        'e1': '23',
        'i2': '34',
        'e2': '45',
        'i3': '56',
        'e3': '67',
    }

    def handle(self, *args, **options):
        if 'filename' in options and options['filename']:
            filename = options['filename']
            if os.path.isfile(filename):
                self.update_structure_annotations(filename)
            else:
                print('Filename not found!')
        else:
            print('No filename specified, aborting')

    def update_structure_annotations(self, filename):

        found = 0
        not_found = 0

        general_cols = ['Title', 'SignProt']
        extra_cols = ['Construct', 'AuxProt', 'PrefChain', 'State']
        segment_cols = [
            ['coord'],
            ['engin'],
            ['SS'],
            ['x50'],
            ['bgr'],
            ['bp', 'bgr'],
            ['bgs'],
            ['bps', 'bgs'],
            ['egr'],
            ['ep', 'egr'],
            ['egs'],
            ['eps', 'egs'],
        ]
        segments = []
        cols = []
        pdbs = OrderedDict()
        with open(filename, 'r') as f:

            for line_num, line in enumerate(f):
                # split the line by tabs
                sline = line.split('\t')
                for col_num, col in enumerate(sline):
                    col = col.strip().strip('"')

                    if line_num == 0:
                        cols.append(col)
                        continue

                    col_title = cols[col_num]
                    segment = False

                    for k, v in self.segment_dict.items():
                        if col_title.startswith(k):
                            segment = v
                            if segment not in segments:
                                segments.append(segment)
                            col_title = col_title[len(k):]
                            break
                    else:
                        if col_title == 'Title':
                            pdb_code = col.split('_')[2]
                            if pdb_code not in pdbs:
                                pdbs[pdb_code] = {}
                                pdbs[pdb_code][col_title] = col

                    if col:
                        if segment and segment not in pdbs[pdb_code]:
                            pdbs[pdb_code][segment] = {}

                        # correct gen. numbers case, e.g. 8X50 -> 8x50
                        col = col.replace('X', 'x')

                        # correct segment name in gen. numbers, e.g. i1x50 -> 12x50
                        if col[:2] in self.numbering_prefix_dict:
                            col = self.numbering_prefix_dict[col[:2]] + col[2:]

                        if segment:
                            pdbs[pdb_code][segment][col_title] = col
                        else:
                            pdbs[pdb_code][col_title] = col

        # update sequence numbers and write line to export file
        export_filename = filename + '.export'
        with open(export_filename, 'w') as f:
            export_line = ''
            for col in general_cols:
                export_line += col + '\t'
            for col in extra_cols:
                export_line += col + '\t'
            for segment in segments:
                for col in segment_cols:
                    export_line += segment + '_' + col[0] + '\t'
            f.write(export_line + '\n')

            for pdb_code, pdb in pdbs.items():
                # fetch structure and protein conformation objects (used for getting sequence numbers later)
                try:
                    structure = Structure.objects.get(pdb_code__index=pdb_code)
                except Structure.DoesNotExist:
                    structure = False
                if structure:
                    pc = structure.protein_conformation

                export_line = ''
                for col in general_cols:
                    export_value = ''
                    if col in pdb:
                        export_value = pdb[col]
                    export_line += export_value + '\t'

                # fetch more columns from db
                try:
                    struct = Structure.objects.get(pdb_code__index=pdb_code);
                    # construct
                    export_line += struct.protein_conformation.protein.entry_name + '\t'
                    # auxiliary protein
                    export_line += ', '.join(struct.stabilizing_agents.values_list('name',
                        flat=True).order_by('id')) + '\t'
                    # preferred chain
                    export_line += struct.preferred_chain + '\t'
                    # state
                    export_line += struct.state.name + '\t'
                except:
                    pass

                for segment in segments:
                    for col in segment_cols:
                        export_value = ''
                        col_name = col[0]

                        # get column value from the original data
                        if col_name in pdb[segment]:
                            export_value = pdb[segment][col_name]
                        # if this is an empty sequence number column, find the sequence number
                        elif len(col) > 1:
                            if (pc and col[1] in pdb[segment] and pdb[segment][col[1]] != '-'
                                and pdb[segment][col[1]] != 'Absent'):
                                res_gn = pdb[segment][col[1]]

                                try:
                                    res = Residue.objects.get(protein_conformation=pc, generic_number__label=res_gn)
                                    export_value = pdbs[pdb_code][segment][col[0]] = res.amino_acid + str(res.sequence_number)
                                    found += 1
                                except Residue.DoesNotExist:
                                    # search for 'adjacent' numbers
                                    sgn = res_gn.split('x')
                                    gen_index = ori_gen_index = int(sgn[1])
                                    gn_found = False
                                    tries = 0

                                    if col_name.startswith('b'):
                                        while not gn_found and tries < 10:
                                            gen_index += 1
                                            search_gn = sgn[0] + 'x' + str(gen_index)
                                            try:
                                                res = Residue.objects.get(protein_conformation=pc,
                                                    generic_number__label=search_gn)
                                                seq_num = res.sequence_number - (gen_index - ori_gen_index)
                                                aa = pc.protein.parent.sequence[seq_num-1]
                                                export_value = pdbs[pdb_code][segment][col[0]] = aa + str(seq_num)
                                                gn_found = True
                                                found += 1
                                            except Residue.DoesNotExist:
                                                tries += 1
                                        else:
                                            if tries == 10:
                                                not_found += 1
                                    if col_name.startswith('e'):
                                        while not gn_found and tries < 10:
                                            gen_index -= 1
                                            search_gn = sgn[0] + 'x' + str(gen_index)
                                            try:
                                                res = Residue.objects.get(protein_conformation=pc,
                                                    generic_number__label=search_gn)
                                                seq_num = res.sequence_number + (ori_gen_index - gen_index)
                                                aa = pc.protein.parent.sequence[seq_num-1]
                                                export_value = pdbs[pdb_code][segment][col[0]] = aa + str(seq_num)
                                                gn_found = True
                                                found += 1
                                            except Residue.DoesNotExist:
                                                tries += 1
                                        else:
                                            if tries == 10:
                                                not_found += 1

                        export_line += export_value + '\t'
                f.write(export_line + '\n')
