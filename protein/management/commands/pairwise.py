from build.management.commands.base_build import Command as BaseBuild
from django.conf import settings
from protein.models import Protein, ProteinSegment, ProteinConformation, ProteinState
from structure.models import Structure, Rotamer
from structure.functions import BlastSearch
from Bio.Blast import NCBIXML, NCBIWWW
from Bio import pairwise2

import subprocess, shlex, os


class Command(BaseBuild):  
	help = 'Blastp search custom dbs'

	def add_arguments(self, parser):
		super(Command, self).add_arguments(parser=parser)
		parser.add_argument('-p1', help='Entry name or sequence of protein 1', default=False, type=str)
		parser.add_argument('-p2', help='Entry name or sequence of protein 2', default=False, type=str)
	
	def handle(self, *args, **options):
		p1 = Protein.objects.filter(entry_name=options['p1'])
		if len(p1)==0:
			seq1 = options['p1']
		else:
			seq1 = p1[0].sequence
		p2 = Protein.objects.filter(entry_name=options['p2'])
		if len(p2)==0:
			seq2 = options['p2']
		else:
			seq2 = p2[0].sequence

		pw2 = pairwise2.align.localms(seq1, seq2, 3, -4, -3, -1)
		ref_seq, temp_seq = str(pw2[0][0]), str(pw2[0][1])
		for r,t in zip(ref_seq,temp_seq):
			print(r,t)