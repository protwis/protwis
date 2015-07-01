from django.core.management.base import BaseCommand, CommandError
from django.conf import settings
from django.db import connection

from protein.models import Protein
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq

import logging,sys,os
from subprocess import Popen,PIPE

class Command(BaseCommand):

    logger = logging.getLogger(__name__)
    
    #The command puts blastdb in the first directory from common static directories
    blast_database_dir = os.sep.join([settings.STATICFILES_DIRS[0], 'blast'])
    help = 'Generates Blast database from sequences deposited within protwis dataset. The sequences are indexed by Protein.Id to speed up searching and make querying easier.'

    def handle(self, *args, **options):

        sequences = []
        print(self.blast_database_dir)
        print('Fetching sequence data and saving it into temp file')
        self.logger.info('Attempting to fetch protein objects...')
        try:
            proteins = Protein.objects.all()
        except Exception as e:
            print(e)
            self.logger.error('Something went wrong while fetching proteins!')
        for protein in proteins:
            sequences.append(SeqRecord(Seq(protein.sequence, IUPAC.protein), id=str(protein.id), description=protein.entry_name))
        self.logger.info('Saving sequences into temp file')
        try:
            if os.path.exists(os.sep.join([self.blast_database_dir, 'temp.fa'])):
                os.unlink(os.sep.join([self.blast_database_dir, 'temp.fa']))
            SeqIO.write(sequences, os.sep.join([self.blast_database_dir, 'temp.fa']), 'fasta')
        except Exception as e:
            print(e)
            self.logger.error('Failed saving the sequences')
        
        print('Running makeblastdb')
        self.logger.info('Issuing makeblastdb')
        try:
            #No need to unlink the previously existing database - makeblastdb overwrites it
            makeblastdb = Popen("makeblastdb -in {} -dbtype prot -title protwis_blastdb -out {} -parse_seqids".format(os.sep.join([self.blast_database_dir, 'temp.fa']), os.sep.join([self.blast_database_dir, 'protwis_blastdb'])), universal_newlines=True, stdout=PIPE, shell=True, stderr=PIPE)
        except Exception as e:
            print(e)
            self.logger.error('Failed to launch makeblastdb')
        out, err = makeblastdb.communicate()
        if len(err) != 0:
            print(err)
        if os.path.exists(os.sep.join([self.blast_database_dir, 'temp.fa'])):
            os.unlink(os.sep.join([self.blast_database_dir, 'temp.fa']))
        print('Success')
