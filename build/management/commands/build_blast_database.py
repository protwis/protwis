from django.core.management.base import BaseCommand, CommandError
from django.conf import settings
from django.db import connection

from protein.models import Protein, ProteinSegment
from residue.models import Residue
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

import logging, sys, os, tempfile
from subprocess import Popen, PIPE

class Command(BaseCommand):
    help = "Generates Blast database from sequences deposited within protwis dataset. The sequences are indexed by " \
        + "Protein.id to speed up searching and make querying easier."

    logger = logging.getLogger(__name__)

    # The command puts blastdb in the first directory from common static directories
    blast_database_dir = os.sep.join([settings.STATICFILES_DIRS[0], "blast"])
    tmp_file_path = os.sep.join([tempfile.gettempdir(), "temp.fa"])


    def handle(self, *args, **options):
        # BLAST database all sequences
        proteins = Protein.objects.filter(sequence_type__slug="wt")
        self.build_database(proteins, "protwis_blastdb")

        # BLAST database all GPCR sequences
        proteins = Protein.objects.filter(sequence_type__slug="wt", family__slug__startswith="00")
        self.build_database(proteins, "protwis_gpcr_blastdb")

        # BLAST database only human sequences
        proteins = Protein.objects.filter(sequence_type__slug="wt", species__common_name="Human")
        self.build_database(proteins, "protwis_human_blastdb")

        # BLAST database only human GPCR sequences
        proteins = Protein.objects.filter(sequence_type__slug="wt", species__common_name="Human", family__slug__startswith="00")
        self.build_database(proteins, "protwis_human_gpcr_blastdb")

        # BLAST database only human GPCR sequences for 7TM bundle
        proteins = Protein.objects.filter(sequence_type__slug="wt", species__common_name="Human", family__slug__startswith="00")

        c = 0
        for protein in proteins:
            # Process sequences to obtain bundle only
            segment_start = ProteinSegment.objects.get(slug="TM1")
            segment_end = ProteinSegment.objects.get(slug="TM7")
            residues = Residue.objects.filter(protein_conformation__protein=protein, protein_segment_id__gte=segment_start.pk, protein_segment_id__lte=segment_end.pk) \
                        .order_by("id").values_list("amino_acid", flat=True)
            protein.sequence = "".join(residues)
            if protein.sequence=="":
                c+=1

        if len(proteins)!=c:
            self.build_database(proteins, "protwis_human_bundle_blastdb")

        self.build_database(os.sep.join([settings.DATA_DIR, 'g_protein_data', 'g_protein_chimeras.fasta']), 'g_protein_chimeras')

    def build_database(self, proteins, blast_db_dir):
        # All sequences
        print("BUILDING BLAST DATABASE", blast_db_dir)
        self.logger.info("BUILDING BLAST DATABASE" + blast_db_dir)
        db_output_path = os.sep.join([self.blast_database_dir, blast_db_dir])

        # fetch sequences
        sequences = []
        if type(proteins)==type(''):
            sequences = SeqIO.parse(open(proteins), 'fasta')
        else:
            for protein in proteins:
                sequences.append(SeqRecord(Seq(protein.sequence), id=str(protein.id),
                    description=protein.entry_name))

        # Write sequences to file
        try:
            if os.path.exists(self.tmp_file_path):
                os.unlink(self.tmp_file_path)
            SeqIO.write(sequences, self.tmp_file_path, "fasta")
            self.logger.info("Saving sequences into {}".format(self.tmp_file_path))
        except Exception as e:
            self.logger.error("Saving the sequences failed")

        # Create BLAST database
        self.logger.info("Running makeblastdb")
        try:
            # No need to unlink the previously existing database - makeblastdb overwrites it
            makeblastdb = Popen("makeblastdb -in {} -dbtype prot -title protwis_blastdb -out {} -parse_seqids".format(
                self.tmp_file_path, db_output_path), universal_newlines=True, stdout=PIPE, shell=True, stderr=PIPE)
            out, err = makeblastdb.communicate()
            if len(err) != 0:
                self.logger.error(err)
                print('ERROR:', blast_db_dir)
        except Exception as e:
            self.logger.error("Makeblastdb failed")

        # remove tmp sequence file
        if os.path.exists(self.tmp_file_path):
            os.unlink(self.tmp_file_path)

        self.logger.info("COMPLETED BUILDING BLAST DATABASE" + blast_db_dir)
