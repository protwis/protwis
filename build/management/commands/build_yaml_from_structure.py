from django.core.management.base import BaseCommand, CommandError
from django.conf import settings
from django.db import connection

from residue.models import Residue
from protein.models import Protein, ProteinConformation, ProteinSegment

from optparse import make_option
from datetime import datetime
import logging, os, yaml


class Command(BaseCommand):
    
    def add_arguments(self, parser):
        parser.add_argument('filename', nargs='+', type=str)

    logger = logging.getLogger(__name__)
    structure_build_data_dir = os.sep.join([settings.DATA_DIR, 'structure_data'])

    csv_fields = {
        'id' : 0, 
        'prot_name' : 1, 
        'class' : 2, 
        'pdb_code' : 3, 
        'endogenous_ligand' :4, 
        'resolution' : 5, 
        'xray_ligand' : 6, 
        'ligand_role' : 7, 
        'chain' : 8, 
        'pubmed_id' : 9, 
        'date' : 10, 
        'G_protein' : 11, 
        'stabilizing_agent' : 12, 
        'n-term' : 13, 
        'icl1' : 14, 
        'ecl1' : 15, 
        'icl2' : 16, 
        'ecl2.1' : 17, 
        'ecl2.2' : 18, 
        'icl3' : 19, 
        'ecl3m' : 20, 
        'c-term' : 21,
        }

    help = 'Creates structures from .csv formatted data strored in the directory {}'.format(structure_build_data_dir)

    def handle(self, *args, **options):

        for arg in options['filename']:
            struct_fh = open(os.sep.join([self.structure_build_data_dir, arg]), 'r')
            for line in struct_fh:
                try:
                    data = [x.strip().strip('"') for x in line.split('\t')]
                    Protein.objects.get(entry_name=data[self.csv_fields['prot_name']])
                    print(data[self.csv_fields['prot_name']])
                    self.create_yaml(data)
                except Protein.DoesNotExist as e:
                    continue

    def create_yaml(self, data):

        yaml_pdb_data = {
            'pdb' : data[self.csv_fields['pdb_code']], 
            'resolution' : data[self.csv_fields['resolution']],
            'publication_date' : data[self.csv_fields['date']][:10],
            'pubmed_id' : data[self.csv_fields['pubmed_id']],                       
            }

        if data[self.csv_fields['ligand_role']] == 'Agonist' and data[self.csv_fields['G_protein']] != 'None':
            state = 'Active'
        else:
            state = 'Inactive'

        yaml_struct_annotations = {
            'ligand' : [{'name' : data[self.csv_fields['xray_ligand']],
                'role' : data[self.csv_fields['ligand_role']]}],
            'state' : state,
            'preferred_chain' : data[self.csv_fields['chain']],
            'g_protein' : data[self.csv_fields['G_protein']],
            }

        yaml_other_data = {
            'construct' : data[self.csv_fields['pdb_code']].lower(),
            'segments' : self.get_segments_data(data[self.csv_fields['prot_name']]),
            }
        out_fh = open('{}.yaml'.format(os.sep.join([self.structure_build_data_dir, 'structures', data[self.csv_fields['pdb_code']]])), 'w')
        out_fh.write('# PDB data\n\n')
        yaml.dump(yaml_pdb_data, out_fh, default_flow_style=False)
        out_fh.write('\n# Structure annotations\n\n')
        yaml.dump(yaml_struct_annotations, out_fh)
        out_fh.write('\n# Structure annotations\n\n')
        yaml.dump(yaml_other_data, out_fh)
        out_fh.close()


        yaml_construct = {
            'name' : data[self.csv_fields['pdb_code']].lower(),
            'protein' : data[self.csv_fields['prot_name']],
            }

        construct_fh = open('{}.yaml'.format(os.sep.join([self.structure_build_data_dir, 'constructs', data[self.csv_fields['pdb_code']]])), 'w')
        yaml.dump(yaml_construct, construct_fh, indent=4)
        construct_fh.close()


    def get_segments_data(self, prot_entry_name):
        output = {}
        segments = ProteinSegment.objects.all()
        for segment in segments:
            resi = list(Residue.objects.filter(protein_segment__slug = segment.slug,
                protein_conformation__protein__entry_name = prot_entry_name).order_by('sequence_number'))
            try:
                if resi:
                    output[segment.slug] = [resi[0].sequence_number, resi[-1].sequence_number]
            except Exception as e:
                output[segment.slug] = ['-,-']
        return output
