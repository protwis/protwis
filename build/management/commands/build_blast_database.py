from django.core.management.base import BaseCommand, CommandError
from django.conf import settings
from django.db import connection

from protein.models import Protein
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

import logging, sys, os, tempfile
from subprocess import Popen, PIPE

class Command(BaseCommand):
    help = 'Generates Blast database from sequences deposited within protwis dataset. The sequences are indexed by ' \
        + 'Protein.id to speed up searching and make querying easier.'

    logger = logging.getLogger(__name__)
    
    # The command puts blastdb in the first directory from common static directories
    blast_database_dir = os.sep.join([settings.STATICFILES_DIRS[0], 'blast'])
    db_file_path = os.sep.join([blast_database_dir, 'protwis_blastdb'])
    human_db_file_path = os.sep.join([blast_database_dir, 'protwis_human_blastdb'])
    tmp_file_path = os.sep.join([tempfile.gettempdir(), 'temp.fa'])

    def handle(self, *args, **options):
        # All sequences
        self.logger.info('BUILDING BLAST DATABASE')

        sequences = []
        self.logger.info('Building blast database in {}'.format(self.blast_database_dir))
        
        # fetch proteins
        proteins = Protein.objects.filter(sequence_type__slug='wt')
        for protein in proteins:
            sequences.append(SeqRecord(Seq(protein.sequence), id=str(protein.id),
                description=protein.entry_name))
        try:
            if os.path.exists(self.tmp_file_path):
                os.unlink(self.tmp_file_path)
            SeqIO.write(sequences, self.tmp_file_path, 'fasta')
            self.logger.info('Saving sequences into {}'.format(self.tmp_file_path))
        except Exception as e:
            self.logger.error('Saving the sequences failed')
        
        self.logger.info('Running makeblastdb')
        try:
            # No need to unlink the previously existing database - makeblastdb overwrites it
            makeblastdb = Popen("makeblastdb -in {} -dbtype prot -title protwis_blastdb -out {} -parse_seqids".format(
                self.tmp_file_path, self.db_file_path), universal_newlines=True, stdout=PIPE, shell=True, stderr=PIPE)
            out, err = makeblastdb.communicate()
            if len(err) != 0:
                self.logger.error(err)
        except Exception as e:
            self.logger.error('Makeblastdb failed')

        
        # remove tmp sequence file
        if os.path.exists(self.tmp_file_path):
            os.unlink(self.tmp_file_path)

        self.logger.info('COMPLETED BUILDING BLAST DATABASE')

        # Human sequences only

        self.logger.info('BUILDING BLAST DATABASE WITH HUMAN SEQUENCES')

        sequences = []
        self.logger.info('Building human blast database in {}'.format(self.blast_database_dir))
        
        # fetch proteins
        proteins = Protein.objects.filter(sequence_type__slug='wt', species__common_name='Human')
        for protein in proteins:
            sequences.append(SeqRecord(Seq(protein.sequence), id=str(protein.id),
                description=protein.entry_name))
        try:
            if os.path.exists(self.tmp_file_path):
                os.unlink(self.tmp_file_path)
            SeqIO.write(sequences, self.tmp_file_path, 'fasta')
            self.logger.info('Saving sequences into {}'.format(self.tmp_file_path))
        except Exception as e:
            self.logger.error('Saving the sequences failed')
        
        self.logger.info('Running makeblastdb')
        try:
            # No need to unlink the previously existing database - makeblastdb overwrites it
            makeblastdb = Popen("makeblastdb -in {} -dbtype prot -title protwis_blastdb -out {} -parse_seqids".format(
                self.tmp_file_path, self.human_db_file_path), universal_newlines=True, stdout=PIPE, shell=True, stderr=PIPE)
            out, err = makeblastdb.communicate()
            if len(err) != 0:
                self.logger.error(err)
        except Exception as e:
            self.logger.error('Makeblastdb failed')

        
        # remove tmp sequence file
        if os.path.exists(self.tmp_file_path):
            os.unlink(self.tmp_file_path)

        self.logger.info('COMPLETED BUILDING BLAST DATABASE WITH HUMAN SEQUENCES')
