from django.core.management.base import BaseCommand, CommandError
from django.conf import settings
from django.utils.text import slugify
from django.utils.html import strip_tags

from protein.models import (Protein, ProteinConformation, ProteinState, ProteinSequenceType, ProteinSegment,ProteinSource)
from residue.models import Residue
from construct.models import *

from optparse import make_option
from datetime import datetime
import logging, os
import yaml


class Command(BaseCommand):
    help = 'Reads source data and creates protein records for constructs'

    def add_arguments(self, parser):
        parser.add_argument('--filename', action='append', dest='filename',
            help='Filename to import. Can be used multiple times')
        parser.add_argument('--purge', action='store_true', dest='purge', default=False,
            help='Purge existing construct records')

    logger = logging.getLogger(__name__)

    # source file directory
    construct_data_dir = os.sep.join([settings.DATA_DIR, 'structure_data', 'constructs'])

    def handle(self, *args, **options):
        # delete any existing construct data
        if options['purge']:
            try:
                self.purge_constructs()
            except Exception as msg:
                print(msg)
                self.logger.error(msg)

        # import the structure data
        try:
            self.create_constructs(options['filename'])
        except Exception as msg:
            print(msg)
            self.logger.error(msg)

    def purge_constructs(self):
        try:
            pst = ProteinSequenceType.objects.get(slug='mod')
            Protein.objects.filter(sequence_type=pst).delete()
        except ProteinSequenceType.DoesNotExist:
            self.logger.warning('ProteinSequenceType mod not found: nothing to delete.')

    def create_constructs(self, filenames):
        self.logger.info('CREATING CONSTRUCTS')

        # what files should be parsed?
        if not filenames:
            filenames = os.listdir(self.construct_data_dir)

        # parse files
        for source_file in filenames:
            source_file_path = os.sep.join([self.construct_data_dir, source_file])
            if os.path.isfile(source_file_path) and source_file[0] != '.':
                self.logger.info('Reading file {}'.format(source_file_path))
                # read the yaml file
                with open(source_file_path, 'r') as f:
                    sd = yaml.load(f, Loader=yaml.FullLoader)

                    # is a protein specified?
                    if 'protein' not in sd:
                        self.logger.error('Protein not specified for construct, skipping')
                        continue

                    # fetch the parent protein
                    try:
                        ppc = ProteinConformation.objects.select_related('protein__family', 'protein__species',
                            'protein__residue_numbering_scheme').get(protein__entry_name=sd['protein'],
                            state__slug=settings.DEFAULT_PROTEIN_STATE)
                    except ProteinConformation.DoesNotExist:
                        # abort if parent protein is not found
                        self.logger.error('Parent protein {} for construct {} not found, aborting!'.format(
                            sd['protein'], sd['name']))
                        continue

                    if not Protein.objects.filter(name=sd['name']).exists():
                        # create a protein record
                        p = Protein()
                        p.parent = ppc.protein
                        p.family = ppc.protein.family
                        p.species = ppc.protein.species
                        p.residue_numbering_scheme = ppc.protein.residue_numbering_scheme
                        p.sequence_type, created = ProteinSequenceType.objects.get_or_create(slug='mod',
                            defaults={'name': 'Modified'})
                        p.source, created = ProteinSource.objects.get_or_create(name='OTHER')
                        p.entry_name = slugify(strip_tags(sd['name']))
                        p.name = sd['name']
                        p.sequence = ppc.protein.sequence
                        # save protein (construct)
                        try:
                            p.save()
                            self.logger.info('Created construct {} with parent protein {}'.format(p.name,
                                ppc.protein.entry_name))
                        except Exception as e:
                            print(e)
                            self.logger.error('Failed creating construct {} with parent protein {}'.format(p.name,
                                ppc.protein.entry_name))
                            continue
                    else:
                        p = Protein.objects.get(name=sd['name'])

                    if not ProteinConformation.objects.filter(protein=p).exists():
                        # create protein conformation record
                        pc = ProteinConformation()
                        pc.protein = p
                        pc.state = ProteinState.objects.get(slug=settings.DEFAULT_PROTEIN_STATE)
                        try:
                            pc.save()
                            self.logger.info('Created conformation {} of protein {}'.format(pc.state.name, p.name))
                        except:
                            self.logger.error('Failed creating conformation {} of protein {}'.format(pc.state.name,
                                p.entry_name))

        self.logger.info('COMPLETED CREATING CONSTRUCTS')
