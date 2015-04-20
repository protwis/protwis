from django.core.management.base import BaseCommand, CommandError
from django.conf import settings

from protein.models import Protein, ProteinSequenceType
from residue.models import Residue

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
                self.puge_constructs()
            except Exception as msg:
                print(msg)
                self.logger.error(msg)
        
        # import the structure data
        try:
            self.create_constructs(options['filename'])
        except Exception as msg:
            print(msg)
            self.logger.error(msg)
    
    def puge_constructs(self):
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
                    sd = yaml.load(f)

                    # fetch the parent protein
                    try:
                        pp = Protein.objects.get(entry_name=sd['protein'])
                    except Protein.DoesNotExist:
                        # abort if parent protein is not found
                        raise Exception('Parent protein {} for construct {} not found, aborting!'.format(
                            sd['protein'], sd['name']))

                    # create a protein record
                    p = Protein()
                    p.parent = pp
                    p.family = pp.family
                    p.species = pp.species
                    p.residue_numbering_scheme = pp.residue_numbering_scheme
                    p.sequence_type, created = ProteinSequenceType.objects.get_or_create(slug='mod',
                        defaults={'name': 'Modified'})
                    p.name = sd['name']
                    p.sequence = pp.sequence
                    
                    # save construct
                    try:
                        p.save()
                        self.logger.info('Created construct {} with parent protein {}'.format(p.name, pp.name))
                    except:
                        self.logger.error('Failed creating construct {} with parent protein {}'.format(p.name, pp.name))
                        continue

                    # create residue records
                    truncations = []
                    for t in sd['truncations']:
                        truncations += list(range(t[0],t[1]+1))

                    mutations = {}
                    for m in sd['mutations']:
                        res_num = m[1:-1]
                        mutations[res_num] = {
                            'wt_res': m[0],
                            'mut_res': m[-1],
                            'full': m,
                        }

                    prs = Residue.objects.filter(protein=pp).prefetch_related(
                        'protein', 'protein_segment', 'generic_number', 'display_generic_number__scheme',
                        'alternative_generic_number__scheme')
                    for pr in prs:
                        if pr.sequence_number not in truncations:
                            r = Residue()
                            r.protein = p
                            r.protein_segment = pr.protein_segment
                            r.generic_number = pr.generic_number
                            r.display_generic_number = pr.display_generic_number
                            r.sequence_number = pr.sequence_number

                            # amino acid, check for mutations
                            if r.sequence_number in mutations:
                                if mutations[r.sequence_number]['wt_res'] == pr.amino_acid:
                                    r.amino_acid = mutations[r.sequence_number]['mut_res']
                                else:
                                    self.logger.error('Mutation {} in construct {} does not match wild-type sequence' \
                                    + ' of {}'.format(mutations[r.sequence_number]['full'], p.name, pp.entry_name))
                            else:
                                r.amino_acid = pr.amino_acid

                            # save residue before populating M2M relations
                            r.save()

                            # alternative generic numbers
                            agns = pr.alternative_generic_number.all()
                            for agn in agns:
                                r.alternative_generic_number.add(agn)

        self.logger.info('COMPLETED CREATING CONSTRUCTS')