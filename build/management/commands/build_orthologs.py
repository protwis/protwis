from django.core.management.base import BaseCommand, CommandError
from django.conf import settings
from django.db import connection

from build.management.commands.build_proteins import Command as BuildProteins
from residue.functions import *
from protein.models import Protein, Gene

import logging
import os
from Bio import AlignIO
from Bio.Align.Applications import ClustalOmegaCommandline
import yaml


class Command(BuildProteins):
    help = 'Reads uniprot text files and creates protein entries of orthologs of human proteins'

    ref_position_source_dir = os.sep.join([settings.DATA_DIR, 'residue_data', 'reference_positions'])

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
                # get human ortholog
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

            # skip if no ortholog is found FIXME use a profile to find a good template
            if not p:
                continue

            # check whether reference positions exist for this protein, and find them if they do not
            ref_position_file_path = os.sep.join([self.ref_position_source_dir, up['entry_name'] + '.yaml'])
            if not os.path.isfile(ref_position_file_path):
                # get reference positions of human ortholog
                template_ref_position_file_path = os.sep.join([self.ref_position_source_dir, p.entry_name + '.yaml'])
                if not os.path.isfile(template_ref_position_file_path):
                    self.logger.error("File {} not found, skipping!".format(template_ref_position_file_path))
                    continue
                template_ref_positions = load_reference_positions(template_ref_position_file_path)


                # write sequences to files
                seq_filename = "/tmp/" + accession + ".fa"
                with open(seq_filename, 'w') as seq_file:
                    seq_file.write("> ref\n")
                    seq_file.write(p.sequence + "\n")
                    seq_file.write("> seq\n")
                    seq_file.write(up['sequence'] + "\n")

                try:
                    ali_filename = "/tmp/out.fa"
                    acmd = ClustalOmegaCommandline(infile=seq_filename, outfile=ali_filename, force=True)
                    stdout, stderr = acmd()
                    a = AlignIO.read(ali_filename, "fasta")
                    self.logger.info("{} aligned to {}".format(up['entry_name'], p.entry_name))
                except:
                    self.logger.error('Alignment failed for {}'.format(up['entry_name']))
                    continue

                # find reference positions
                ref_positions = {}
                ref_positions_in_ali = {}
                for position_generic_number, rp in template_ref_positions.items():
                    gaps = 0
                    for i, r in enumerate(a[0].seq, 1):
                        if r == "-":
                            gaps += 1
                        if i-gaps == rp:
                            ref_positions_in_ali[position_generic_number] = i
                for position_generic_number, rp in ref_positions_in_ali.items():
                    gaps = 0
                    for i, r in enumerate(a[1].seq, 1):
                        if r == "-":
                            gaps += 1
                        if i == rp:
                            ref_positions[position_generic_number] = i - gaps

                # write reference positions to a file
                with open(ref_position_file_path, "w") as ref_position_file:
                    yaml.dump(ref_positions, ref_position_file, default_flow_style=False)

            # create a database entry for the protein
            self.create_protein(p.name, p.family, p.residue_numbering_scheme, accession, up)

        self.logger.info('COMPLETED ORTHOLOGS')