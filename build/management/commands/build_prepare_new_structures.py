from django.core.management.base import BaseCommand, CommandError
from django.conf import settings
from django.db import connection

from structure.management.commands.list_new_xtals import QueryPDB
from residue.models import Residue
from protein.models import Protein, ProteinConformation, ProteinSegment
from structure.sequence_parser import SequenceParser

from Bio.PDB import parse_pdb_header

import logging, os, urllib, yaml

class Command(BaseCommand):

    logger = logging.getLogger(__name__)

    # source file directory
    structure_data_dir = os.sep.join([settings.DATA_DIR, 'structure_data', 'structures'])
    structure_build_data_dir = os.sep.join([settings.DATA_DIR, 'structure_data'])
    pdb_data_dir = os.sep.join([settings.DATA_DIR, 'structure_data', 'pdbs'])

    def handle(self, *args, **options):

        q = QueryPDB()
        q.list_xtals(verbose=False)

        for record in q.new_structures:
            pdb_code = record[0]
            wt_id = Protein.objects.get(entry_name=record[1]).id
            if not os.path.exists(os.sep.join([self.pdb_data_dir, "{}.pdb".format(pdb_code)])):
                self.download_pdb(pdb_code)
            self.parser = SequenceParser(os.sep.join([self.pdb_data_dir, "{}.pdb".format(pdb_code)]), wt_protein_id=wt_id)
            header = parse_pdb_header(os.sep.join([self.pdb_data_dir, "{}.pdb".format(pdb_code)]))
            self.create_yaml(pdb_code, record[1], header)


    def download_pdb(self, pdb_code):

        url = "https://www.rcsb.org/pdb/files/{}.pdb".format(pdb_code)
        urllib.request.urlretrieve(url, os.sep.join([self.pdb_data_dir, "{}.pdb".format(pdb_code)]))


    def create_yaml(self, pdb_code, prot_name, data):

        yaml_pdb_data = {
            'pdb' : pdb_code,
            'resolution' : data['resolution'],
            'publication_date' : data['release_date'][:10],
            #PDB header contains full citation
            #'pubmed_id' : data['journal_reference'],
            }
        yaml_struct_annotations = {
            'fusion_proteins' : { x : [y[0], y[1]] for x,y in enumerate(self.parser.fusions) }
            }
        yaml_other_data = {
            'construct' : pdb_code.lower(),
            'segments' : self.get_segments_data(prot_name),
            'segments_in_structure' : self.parser.get_segments(),
            }

        out_fh = open('{}.yaml'.format(os.sep.join([self.structure_build_data_dir, 'structures', pdb_code])), 'w')
        out_fh.write('# PDB data\n\n')
        yaml.dump(yaml_pdb_data, out_fh, default_flow_style=False)
        out_fh.write('\n# Structure annotations\n\n')
        yaml.dump(yaml_struct_annotations, out_fh)
        out_fh.write('\n# Structure annotations\n\n')
        yaml.dump(yaml_other_data, out_fh)
        out_fh.close()


        yaml_construct = {
            'name' : pdb_code.lower(),
            'protein' : prot_name,
            }

        construct_fh = open('auto_{}.yaml'.format(os.sep.join([self.structure_build_data_dir, 'constructs', pdb_code])), 'w')
        yaml.dump(yaml_construct, construct_fh, indent=4)
        construct_fh.close()


    def get_segments_data(self, prot_entry_name):
        output = {}
        segments = ProteinSegment.objects.all()
        for segment in segments:
            resi = list(Residue.objects.filter(protein_segment = segment,
                protein_conformation__protein__entry_name = prot_entry_name).order_by('sequence_number'))
            try:
                if resi:
                    output[segment.slug] = [resi[0].sequence_number, resi[-1].sequence_number]
            except Exception as e:
                output[segment.slug] = ['-,-']
        return output
