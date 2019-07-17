from django.core.management.base import BaseCommand, CommandError
from build.management.commands.base_build import Command as BaseBuild
from django.conf import settings
from django.db import connection
from django.db import IntegrityError

from protein.models import (Protein, ProteinGProtein,ProteinGProteinPair, ProteinConformation, ProteinState, ProteinFamily, ProteinAlias,
        ProteinSequenceType, Species, Gene, ProteinSource, ProteinSegment)

from residue.models import (ResidueNumberingScheme, ResidueGenericNumber, Residue, ResidueGenericNumberEquivalent)

from signprot.models import SignprotComplex

import re
import os
from Bio import pairwise2
from collections import OrderedDict
import logging
import shlex, subprocess
from io import StringIO
from Bio.PDB import PDBParser,PPBuilder
from Bio import pairwise2
import yaml
import pprint


AA = {'ALA':'A', 'ARG':'R', 'ASN':'N', 'ASP':'D',
     'CYS':'C', 'GLN':'Q', 'GLU':'E', 'GLY':'G',
     'HIS':'H', 'ILE':'I', 'LEU':'L', 'LYS':'K',
     'MET':'M', 'PHE':'F', 'PRO':'P', 'SER':'S',
     'THR':'T', 'TRP':'W', 'TYR':'Y', 'VAL':'V', 
     'YCM':'C', 'CSD':'C', 'TYS':'Y', 'SEP':'S'} #non-standard AAs


class Command(BaseBuild):
	wt_annotation_file = os.sep.join([settings.DATA_DIR, 'structure_data', 'annotation', 'ECD_wt.yaml'])

	def add_arguments(self, parser):
		parser.add_argument('--filename', action='append', dest='filename',
							help='Filename to import. Can be used multiple times')
		parser.add_argument('--wt', default=False, type=str, help='Add wild type protein sequence to data')
		parser.add_argument('--xtal', default=False, type=str, help='Add xtal to data')
		parser.add_argument('--build_datafile', default=False, action='store_true', help='Build PDB_UNIPROT_ENSEMBLE_ALL file')
		parser.add_argument('--purge', default=False, action='store_true', help='Purge G protein structures from database')

	def handle(self, *args, **options):
		self.options = options
		self.build_scheme()
		self.build_segments()
		self.build_residues()

	def build_scheme(self):
		self.scheme, created = ResidueNumberingScheme.objects.get_or_create(slug='ecd', short_name='ECD', name='Class B GPCR Extracellular domain', parent=None)

	def build_segments(self):
		self.segments = ['H1', 'h1b1', 'B1', 'b1b2', 'B2', 'b2b3', 'B3', 'b3h2', 'H2', 'h2b4', 'B4', 'b4h3', 'H3', 'h3tm1']
		self.segments = ProteinSegment.objects.filter(name__startswith='ECD')
		self.segments_dict = OrderedDict()
		for s in self.segments:
			if s[0]=='H':
				category = 'helix'
			elif s[0]=='B':
				category = 'sheet'
			else:
				category = 'loop'
			seg = ProteinSegment.objects.get(slug=s, name='ECD '+s, proteinfamily='GPCR')
			self.segments_dict[s] = seg

	def build_ECD_gn(self, gn, segment):
		ecd_gn, created = ResidueGenericNumber.objects.get_or_create(label=gn, protein_segment=segment, scheme=self.scheme)
		return ecd_gn

	def build_residues(self):
		with open(self.wt_annotation_file, 'r') as f:
			wt_annotation = yaml.load(f)
		# wt_annotation = OrderedDict([('calcr_human', wt_annotation['calcr_human'])])
		
		for entry_name, val in wt_annotation.items():
			protein = Protein.objects.get(entry_name=entry_name)
			prot_conf = ProteinConformation.objects.get(protein=protein)
			residues = Residue.objects.filter(protein_conformation=prot_conf, protein_segment__slug__in=['N-term', 'TM1'])
			annotation_dict = OrderedDict()
			print(entry_name, val)
			for i, seg in enumerate(self.segments):
				if seg[0] in ['H','B']:
					segment = OrderedDict()
					x = 50
					for i in range(int(val[seg+'x50']), int(val[seg+'b'])-1, -1):
						segment[i] = seg+'x'+str(x)
						x-=1
					x = 50
					for i in range(int(val[seg+'x50']), int(val[seg+'e'])+1):
						segment[i] = seg+'x'+str(x)
						x+=1
					key_order = sorted(list(segment.keys()))
					ordered_segment = OrderedDict([(k, segment[k]) for k in key_order])
					annotation_dict[seg] = ordered_segment
				else:
					segment = OrderedDict()
					if seg=='h3tm1':
						H3e = int(val['H3e'])
						TM1b = residues.filter(protein_segment__slug='TM1')[0].sequence_number
						for i in range(H3e+1, TM1b):
							segment[i] = None
					else:
						for j in range(int(val[self.segments[i-1]+'e'])+1, int(val[self.segments[i+1]+'b'])):
							segment[j] = None
					annotation_dict[seg] = segment
			

			# pprint.pprint(annotation_dict)
			# print(protein, prot_conf)
			# print(residues)
			for seg, val in annotation_dict.items():
				for num, gn in val.items():
					res = residues.get(sequence_number=num)
					if gn:
						gn = self.build_ECD_gn(gn, self.segments_dict[seg])
						res.display_generic_number = gn
						res.generic_number = gn
					res.protein_segment = self.segments_dict[seg]
					res.save()
					# print(res,gn)
