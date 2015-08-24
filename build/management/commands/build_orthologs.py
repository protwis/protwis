from django.core.management.base import BaseCommand, CommandError
from django.conf import settings
from django.db import connection

from build.management.commands.build_human_proteins import Command as BuildHumanProteins
from residue.functions import *
from protein.models import Protein, Gene

import logging
import os
import yaml


class Command(BuildHumanProteins):
    help = 'Reads uniprot text files and creates protein entries of orthologs of human proteins'

    ref_position_source_dir = os.sep.join([settings.DATA_DIR, 'residue_data', 'reference_positions'])
    auto_ref_position_source_dir = os.sep.join([settings.DATA_DIR, 'residue_data', 'auto_reference_positions'])

    def handle(self, *args, **options):
        # create proteins
        try:
            self.create_orthologs()
        except Exception as msg:
            print(msg)
            self.logger.error(msg)

    def create_orthologs(self):
        self.logger.info('CREATING ORTHOLOGS')

        # parse files
        filenames = os.listdir(self.local_uniprot_dir)
        for source_file in filenames:
            source_file_name = os.sep.join([self.local_uniprot_dir, source_file])
            self.logger.info('Processing accession ' + source_file)
            accession = source_file.split(".")[0]
            up = self.parse_uniprot_file(accession)

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
            auto_ref_position_file_path = os.sep.join([self.ref_position_source_dir, up['entry_name'] + '.yaml'])
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