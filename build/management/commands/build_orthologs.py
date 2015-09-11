from django.core.management.base import BaseCommand, CommandError
from django.conf import settings
from django.db import connection
from django.db.models import Q

from build.management.commands.build_human_proteins import Command as BuildHumanProteins
from residue.functions import *
from protein.models import Protein, Gene

import logging
import os
import yaml
from optparse import make_option


class Command(BuildHumanProteins):
    help = 'Reads uniprot text files and creates protein entries of orthologs of human proteins'

    def add_arguments(self, parser):
        parser.add_argument('--only-constructs', action='store_true', dest='only_constructs',
            help='Only import orthologs for which there are annotated constructs. Useful for building a small ' \
            + 'database with all structures')
        parser.add_argument('--purge', action='store_true', dest='purge', default=False,
            help='Purge existing orthologs records')

    ref_position_source_dir = os.sep.join([settings.DATA_DIR, 'residue_data', 'reference_positions'])
    auto_ref_position_source_dir = os.sep.join([settings.DATA_DIR, 'residue_data', 'auto_reference_positions'])
    construct_data_dir = os.sep.join([settings.DATA_DIR, 'structure_data', 'constructs'])
    uniprot_url = 'http://www.uniprot.org/uniprot/?query={}&columns=id&format=tab'

    def handle(self, *args, **options):
        if options['purge']:
            try:
                self.purge_orthologs()
            except:
                self.logger.error('Could not purge orthologs')

        if options['only_constructs']:
            only_constructs = True
        else:
            only_constructs = False
        
        try:
            self.create_orthologs(only_constructs)
        except Exception as msg:
            print(msg)
            self.logger.error(msg)

    def purge_orthologs(self):
        Protein.objects.filter(~Q(species__id=1)).delete()

    def create_orthologs(self, only_constructs):
        self.logger.info('CREATING ORTHOLOGS')

        construct_entry_names = []
        if only_constructs:
            # go through constructs and finding their entry names for lookup
            self.logger.info('Getting construct accession codes')
            filenames = os.listdir(self.construct_data_dir)
            for source_file in filenames:
                source_file_path = os.sep.join([self.construct_data_dir, source_file])
                self.logger.info('Getting protein name from construct file {}'.format(source_file))
                split_filename = source_file.split(".")
                construct_name = split_filename[0]
                extension = split_filename[1]
                if extension != 'yaml':
                    continue

                # read the yaml file
                with open(source_file_path, 'r') as f:
                    sd = yaml.load(f)

                # check whether protein is specified
                if 'protein' not in sd:
                    continue

                # append entry_name to lookup list
                construct_entry_names.append(sd['protein'])

        # parse files
        filenames = os.listdir(self.local_uniprot_dir)
        for source_file in filenames:
            source_file_name = os.sep.join([self.local_uniprot_dir, source_file])
            self.logger.info('Processing accession ' + source_file)
            split_filename = source_file.split(".")
            accession = split_filename[0]
            extension = split_filename[1]
            if extension != 'txt':
                continue

            up = self.parse_uniprot_file(accession)

            # skip human proteins
            if 'species_latin_name' in up and up['species_latin_name'] == 'Homo sapiens':
                continue

            # should proteins that are not constructs be skipped?
            if only_constructs and up['entry_name'] not in construct_entry_names:
                continue

            # is there already an entry for this protein?
            try:
                p = Protein.objects.get(entry_name=up['entry_name'])
                continue
            except Protein.DoesNotExist:
                p = None
                # get human ortholog using gene name
                for gene in up['genes']:
                    try:
                        g = Gene.objects.get(name__iexact=gene, species__id=1, position=0)
                        ps = g.proteins.all().order_by('id')
                        p = ps[0]
                        self.logger.info("Human ortholog found: {}".format(p.entry_name))
                        break
                    except Gene.DoesNotExist:
                        self.logger.info("No gene found for {}".format(gene))
                        continue
                # if gene name not found, try using entry name
                if not p:
                    split_entry_name = up['entry_name'].split('_')

                    # add _ to the split entry name to avoid e.g. gp1 matching gp139
                    entry_name_query = split_entry_name[0] + '_'
                    try:
                        p = Protein.objects.get(entry_name__startswith=entry_name_query, species__id=1)
                        self.logger.info("Human ortholog found: {}".format(p.entry_name))
                    except Protein.DoesNotExist:
                        self.logger.info("No match found for {}".format(entry_name_query))

            # skip if no ortholog is found FIXME use a profile to find a good template
            if not p:
                continue

            # check whether reference positions exist for this protein, and find them if they do not
            ref_position_file_path = os.sep.join([self.ref_position_source_dir, up['entry_name'] + '.yaml'])
            auto_ref_position_file_path = os.sep.join([self.auto_ref_position_source_dir, up['entry_name'] + '.yaml'])
            if not os.path.isfile(ref_position_file_path):
                # look for the file in the automatically generated reference file dir
                if not os.path.isfile(auto_ref_position_file_path):
                    # get reference positions of human ortholog
                    template_ref_position_file_path = os.sep.join([self.ref_position_source_dir,
                        p.entry_name + '.yaml'])
                    ref_positions = align_protein_to_reference(up, template_ref_position_file_path, p)

                    # write reference positions to a file
                    with open(auto_ref_position_file_path, "w") as auto_ref_position_file:
                        yaml.dump(ref_positions, auto_ref_position_file, default_flow_style=False)

            # create a database entry for the protein
            self.create_protein(p.name, p.family, p.sequence_type, p.residue_numbering_scheme, accession, up)

        self.logger.info('COMPLETED ORTHOLOGS')