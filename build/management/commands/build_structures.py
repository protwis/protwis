from django.core.management.base import BaseCommand, CommandError
from django.conf import settings
from django.db import connection

from protein.models import Protein
from residue.models import Residue
from common.models import WebLink, WebResource, Publication
from structure.models import Structure, StructureType

from optparse import make_option
from datetime import datetime
import logging, os


class Command(BaseCommand):

    logger = logging.getLogger(__name__)
    structure_build_data_dir = os.sep.join([settings.DATA_DIR, 'structure_data'])

    csv_fields = {
        'id' : 0, 
        'prot_name' : 1, 
        'class' : 2, 
        'pdb_code' : 3, 
        'endogenous_ligand' :4, 
        'resolution' : 5, 
        'x-ray_ligand' : 6, 
        'ligand_function' : 7, 
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
        # delete any existing protein data
        try:
            self.truncate_structure_tables()
        except Exception as msg:
            # print(msg)
            self.logger.error(msg)
        self.create_structures(args)
    
    def truncate_structure_tables(self):
        cursor = connection.cursor()
        
        tables_to_truncate = [
            #Following the changes in the models - SM
            'structure',                
        ]

        for table in tables_to_truncate:
            cursor.execute("TRUNCATE TABLE {!s} CASCADE".format(table))

    def create_structures(self, args):
        self.logger.info('BUILDING UP STRUCTURE RECORDS')
        for arg in args:
            structures = []
            if os.path.exists(os.sep.join([self.structure_build_data_dir, arg])):
                structures = self.parse_csv_data(os.sep.join([self.structure_build_data_dir, arg]))
                self.logger.info('USING DATA FROM {} FILE'.format(arg))

                for structure in structures:
                    s = Structure() 
                    s.resolution = structure[self.csv_fields['resolution']]
                    s.pdb_publication_date = "{!s}".format(datetime.strptime(structure[self.csv_fields['date']], "%Y-%m-%d %H:%M:%S").date())
                    s.preferred_chain = structure[self.csv_fields['chain']]
                    try:
                        s.protein = Protein.objects.get(entry_name=structure[self.csv_fields['pdb_code']])
                    except Protein.DoesNotExist:
                        print("Failed to save the structure {} Protein {} does not exist.".format(structure[self.csv_fields['pdb_code']], structure[self.csv_fields['prot_name']]))
                        self.logger.error("Failed to save the structure {} Protein {} does not exist.".format(structure[self.csv_fields['pdb_code']], structure[self.csv_fields['prot_name']]))
                        continue
                    #We assume that the proper web resources are defined already
                    try:
                        s.pdb_code = WebLink.objects.get(index=structure[self.csv_fields['pdb_code']])
                    except WebLink.DoesNotExist:
                        code = WebLink(index=structure[self.csv_fields['pdb_code']], web_resource = WebResource.objects.get(slug='pdb'))
                        code.save()
                        s.pdb_code = code
                    #FIXME: Temporary, since at first we are just using crystals
                    try:
                        s.structure_type = StructureType.objects.get(slug='x-ray')
                    except StructureType.DoesNotExist as msg:
                        print(msg)
                        xray = StructureType(slug='x-ray', description='X-Ray Diffraction')
                        xray.save()
                        s.structure_type = xray
                    try:
                        s.publication = Publication.objects.get(web_link__index=structure[self.csv_fields['pubmed_id']])
                    except:
                        p = Publication()
                        try:
                            p.web_link = WebLink.objects.get(index=structure[self.csv_fields['pubmed_id']], web_resource__slug='pubmed')
                        except WebLink.DoesNotExist:
                            code = WebLink(index=structure[self.csv_fields['pubmed_id']], web_resource = WebResource.objects.get(slug='pubmed'))
                            code.save()
                            p.web_link = code
                        p.update_from_pubmed_data(index=structure[self.csv_fields['pubmed_id']])
                        p.save()
                        s.publication = p

                    try:
                        s.save()
                    except Exception as e:
                        print(e)
                        self.logger.error("Failed to save the structure {}".format(structure[self.csv_fields['pdb_code']]))


    def parse_csv_data(self, file_name):
        print('Reading up the structures from file {}'.format(file_name))
        structures_fh = open(file_name, 'r')
        structures = []
        for line in structures_fh:
            structures.append([x.strip().strip('"') for x in line.split(';')])
        print('Done. {:n} structure records parsed'.format(len(structures)))
        return structures
