from django.core.management.base import BaseCommand, CommandError
from django.conf import settings

from residue.models import Residue
from structure.models import Structure

import os
import yaml
from collections import OrderedDict

class Command(BaseCommand):
    help = 'Reads a structure annotation file exported from Maestro, and exports as multiple yaml files'

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

    def add_arguments(self, parser):
        parser.add_argument('-f', '--filename',
            action='store',
            dest='filename',
            help='Path to Uniprot text file')

    # source file directory
    structure_data_dir = os.sep.join([settings.DATA_DIR, 'structure_data', 'structures'])

    def handle(self, *args, **options):
        # self.update_existing_files()
        if 'filename' in options and options['filename']:
            filename = options['filename']
            if os.path.isfile(filename):
                self.read_structure_annotations(filename)
            else:
                print('Filename not found!')
        else:
            print('No filename specified, aborting')

    def update_existing_files(self):
        filenames = os.listdir(self.structure_data_dir)

        for filename in filenames:
            if filename.endswith('yaml'):
                file_path = os.sep.join([self.structure_data_dir, filename])
                with open(file_path) as f:
                    lines = f.read().splitlines()

                with open(file_path, 'w') as f:
                    for line in lines:
                        if line.startswith('g_protein'):
                            line = line.replace('g_protein', 'signaling_protein')
                        if not line.startswith('representative'):
                            f.write(line + '\n')

    def read_structure_annotations(self, filename):
        cols = []
        pdbs = OrderedDict()
        non_segment_specific_cols = ['SignProt', '!Representative']
        with open(filename, 'r') as f:
            for line_num, line in enumerate(f):
                # split the line by tabs
                sline = line.split('\t')

                # dict to hold exportable row values before they are assigned to a structure
                row_values = {}

                pdb_code = False
                segment = False

                for col_num, col in enumerate(sline):
                    col = col.strip('"')

                    if line_num == 0:
                        cols.append(col)
                        continue

                    col_title = cols[col_num]
                    processed_col_title = col_title.strip('!')

                    if col_title == '!Title':
                        pdb_code = col.split('_')[2]
                        if pdb_code not in pdbs:
                            pdbs[pdb_code] = OrderedDict()
                    elif col_title == 'Seg':
                        segment = col
                        pdbs[pdb_code][segment] = OrderedDict()
                    elif pdb_code and col_title in non_segment_specific_cols:
                        pdbs[pdb_code][processed_col_title] = col
                    elif segment and col and col != '-' and col_title.startswith('!'):
                        pdbs[pdb_code][segment][processed_col_title] = col

        # export a yaml file for each structure
        for pdb_code, pdb in pdbs.items():

            # read existing yaml file and extract data
            old_filename = pdb_code + '.yaml'
            old_file_path = os.sep.join([self.structure_data_dir, old_filename])
            existing_data = {}
            if os.path.isfile(old_file_path):
                with open(old_file_path, 'r') as f:
                    existing_data = yaml.load(f, Loader=yaml.FullLoader)
                sd = OrderedDict()
                cols_to_extract = ['construct', 'pdb', 'representative', 'preferred_chain', 'state',
                    'endogenous_ligand', 'ligand', 'auxiliary_protein', 'signaling_protein']
                for col in cols_to_extract:
                    if col in existing_data:
                        sd[col] = existing_data[col]

            # make structure representative
            if pdb['Representative']:
                sd['representative'] = True

            # Add SignProt data
            if 'SignProt' in pdb and pdb['SignProt']:
                sd['signaling_protein'] = pdb['SignProt']

            # add segment border data
            sd['segments'] = OrderedDict()
            sd['segments_in_structure'] = OrderedDict()
            for segment, segment_data in pdb.items():
                if 'bAlP' in segment_data and 'eAlP' in segment_data:
                    borders = [int(segment_data['bAlP'][1:]), int(segment_data['eAlP'][1:])]
                    sd['segments'][segment] = borders
                borders = []
                if 'bStP' in segment_data:
                    borders.append(int(segment_data['bStP'][1:]))
                else:
                    borders.append(0)
                if 'eStP' in segment_data:
                    borders.append(int(segment_data['eStP'][1:]))
                else:
                    borders.append(0)
                if borders and borders != [0] and borders != [0,0]:
                    sd['segments_in_structure'][segment] = borders
                # anomalies, GPCR specific!
                if isinstance(segment_data, dict):
                    for key, value in segment_data.items():
                        if key.startswith('dist'):
                            if len(value) == 5: # bulge
                                if 'bulges' not in sd:
                                    sd['bulges'] = OrderedDict()
                                if segment not in sd['bulges']:
                                    sd['bulges'][segment] = []
                                sd['bulges'][segment].append(value)
                            if len(value) == 4: # constriction
                                if 'constrictions' not in sd:
                                    sd['constrictions'] = OrderedDict()
                                if segment not in sd['constrictions']:
                                    sd['constrictions'][segment] = []
                                sd['constrictions'][segment].append(value)
                            

            # add distortion data
            sd['distortion'] = OrderedDict()
            for segment, segment_data in pdb.items():
                if 'distortion' in segment_data:
                    sd['distortion'][segment] = segment_data['distortion']

            # yaml filename
            export_filename = pdb_code + '.yaml'
            export_file_path = os.sep.join([self.structure_data_dir, export_filename])

            # export as yaml
            with open(export_file_path, 'w') as f:
                yaml.dump(sd, f)

