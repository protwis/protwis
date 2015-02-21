from django.core.management.base import BaseCommand, CommandError
from django.conf import settings
from django.db import connection

from protein.models import Protein
from residue.models import Residue
from publication.models import Publication
from common.models import WebLink, WebResource
from structure.models import Structure, StructureType

from optparse import make_option
import logging, os


class Command(BaseCommand):

    logger = logging.getLogger(__name__)
    structure_build_data_dir = os.sep.join([settings.DATA_DIR, 'structure_data'])

    help = 'Creates structures from .csv formatted data strored in the directory {}'.format(structure_build_data_dir)

    def handle(self, *args, **options):
        

        pass

    def create_structures(self, args):
        self.logger.info('BUILDING UP STRUCTURE RECORDS')
        for arg in args:
            structures = []
            if os.path.exists(os.sep.join([self.generic_numbers_source_dir, arg])):
                structures = self.parse_csv_data(os.sep.join([self.generic_numbers_source_dir, arg]))
                self.logger.info('USING DATA FROM {} FILE'.format(arg))

                for structure in structures:
                    s = Structure() 
                    s.resolution = structure[2]
                    s.pdb_publication_date = structure[6]
                    s.preferred_chain = structure[7]
                    s.save()
                    s.protein = Protein.objects.get(entry_name=structure[0])

                    #We assume that the proper web resources are defined already
                    try:
                        s.pdb_code = WebLink.objects.get(index=structure[1])
                    except WebLink.DoesNotExist as msg:
                        print(msg)
                        code = WebLink(index=structure[1], web_resource = WebResource.objects.get(slug='pdb'))
                    #FIXME: Temporary, since at first we are just using crystals
                    try:
                        s.structure_type = StructureType.objects.get(slug='x-ray')
                    except StructureType.DoesNotExist as msg:
                        print(msg)
                        xray = StructureType(slug='x-ray', name='X-Ray Diffraction')
                        xray.save()
                        s.structure_type = xray
                    try:
                        s.publication = Publication.objects.get(weblink.index=structure[5])
                    except:
                        p = Publication()
                        p.update_from_pubmed_data(index=structure[5])
                        p.save()
                        s.publication = p

                    try:
                        s.save()
                    except Exception as e:
                        print(e)
                        self.logger.error("Failed to save the structure {}".format(structure[1]))
    def parse_csv_data(self, file_name):
        print('Reading up the structures from file {}'.format(file_name))
        structures_fh = open(file_name, 'r')
        structures = []
        for line in structures_fh:
            #The fileds in the csv I use are:
            #prot_name, pdb_code, resolution, endogenous_ligand, crystalized_ligand, pubmed_id, pub_date, preferred_chain
            structures.append([x.strip('"') for x in line.split(',')])
        print('Done. {:n} structure records parsed'.format(len(structures)))
        return structures
