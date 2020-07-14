from django.core.management.base import BaseCommand, CommandError
from django.conf import settings

from protein.models import ProteinConformation
from residue.models import Residue

import os
import yaml
from collections import OrderedDict

class Command(BaseCommand):
    help = 'Updates x50 positions in loops'

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
    export_dir_path = os.sep.join([settings.DATA_DIR, 'residue_data', 'reference_positions'])

    def handle(self, *args, **options):
        self.update_loops()

    def update_loops(self):
        # fetch all wild-type proteins
        pcs = ProteinConformation.objects.filter(protein__sequence_type__slug='wt', protein__species_id=1,
            state__slug=settings.DEFAULT_PROTEIN_STATE).select_related('protein')

        for pc in pcs:
            export_data = {}

            # get reference position file to update
            ref_file_path = os.sep.join([self.export_dir_path, pc.protein.entry_name + '.yaml'])
            try:
                with open(ref_file_path) as ref_file:
                    export_data = yaml.load(ref_file, Loader=yaml.FullLoader)
            except:
                print('Error opening file {}'.format(ref_file_path))
                continue

            # order data
            ordered_export_data = OrderedDict(sorted(export_data.items(), key=lambda x: x[0]))

            # ICL1, ECL1, ICL2 (class dependent)
            icl1_ref = False
            ecl1_ref = False
            icl2_ref = False
            rs = Residue.objects.filter(protein_conformation=pc).select_related('generic_number')
            for res_num, r in enumerate(rs):
                if not pc.protein.family.slug.startswith('003'):
                    # ICL1
                    if r.generic_number and r.generic_number.label == '1x50':
                        icl1_ref = r.sequence_number
                    if icl1_ref and r.protein_segment.id == 3 and r.sequence_number == icl1_ref+13:
                        ordered_export_data['12x50'] = r.sequence_number

                if pc.protein.family.slug.startswith('001'):
                    # ECL1
                    try:
                        res_plus_two = rs[res_num+2]
                    except:
                        pass
                    else:
                        if (not ecl1_ref and r.protein_segment.id == 5 and r.amino_acid == 'W'
                            and res_plus_two.protein_segment.id == 5
                            and res_plus_two.amino_acid in AMINO_ACID_GROUPS['hp']):
                            ordered_export_data['23x50'] = ecl1_ref = r.sequence_number
                    # ICL2
                    if r.generic_number and r.generic_number.label == '3x50':
                        icl2_ref = r.sequence_number
                    if icl2_ref and r.protein_segment.id == 7 and r.sequence_number == icl2_ref+7:
                        ordered_export_data['34x50'] = r.sequence_number

            # ECL2 for all classes
            try:
                tm3_cys = Residue.objects.get(protein_conformation=pc, generic_number__label='3x25')
                if tm3_cys.amino_acid == 'C':
                    ecl_2_residues = Residue.objects.filter(protein_conformation=pc, protein_segment_id=9)
                    for residue in ecl_2_residues:
                        if residue.amino_acid == 'C':
                            ordered_export_data['45x50'] = residue.sequence_number
            except Residue.DoesNotExist:
                pass

            # H8 for all classes
            h8_ref = False
            for res_num, r in enumerate(rs):
                if r.generic_number and r.generic_number.label == '7x50':
                    h8_ref = r.sequence_number
                if h8_ref and r.protein_segment.id == 16 and r.sequence_number == h8_ref+10:
                    ordered_export_data['8x50'] = r.sequence_number


            # open export file for writing
            with open(ref_file_path, 'w') as ref_file:
                yaml.dump(ordered_export_data, ref_file, default_flow_style=False)
