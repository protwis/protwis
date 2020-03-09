from django.core.management.base import BaseCommand, CommandError
from build.management.commands.base_build import Command as BaseBuild
from django.conf import settings
from django.db import connection
from django.db import IntegrityError

from protein.models import (Protein, ProteinGProtein,ProteinGProteinPair, ProteinConformation, ProteinState, ProteinFamily, ProteinAlias,
        ProteinSequenceType, Species, Gene, ProteinSource, ProteinSegment, ProteinAnomaly)

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
	anomalies_file = os.sep.join([settings.DATA_DIR, 'structure_data', 'annotation', 'ECD_anomalies.yaml'])

	def add_arguments(self, parser):
		pass

	def handle(self, *args, **options):
		self.options = options
		self.build_scheme()
		self.get_segments()
		self.build_ECD_anomalies()
		self.build_residues()

	def build_scheme(self):
		self.scheme, created = ResidueNumberingScheme.objects.get_or_create(slug='ecd', short_name='ECD', name='Class B GPCR Extracellular domain', parent=None)

	def get_segments(self):
		self.segments = ProteinSegment.objects.filter(name__startswith='ECD')
		self.segments_dict = OrderedDict()
		for s in self.segments:
			self.segments_dict[s.slug] = s

	def build_ECD_anomalies(self):
		with open(self.anomalies_file, 'r') as fa:
			self.anomalies = yaml.load(fa)

	def build_ECD_gn(self, gn, segment):
		ecd_gn, created = ResidueGenericNumber.objects.get_or_create(label=gn, protein_segment=segment, scheme=self.scheme)
		return ecd_gn

	def build_residues(self):
		with open(self.wt_annotation_file, 'r') as f:
			wt_annotation = yaml.load(f, Loader=yaml.FullLoader)
		# wt_annotation = OrderedDict([('sctr_human', wt_annotation['sctr_human'])])
		with open(self.B1_annotation_file, 'r') as fB1:
			B1_annotation = yaml.load(fB1, Loader=yaml.FullLoader)
		
		for entry_name, val in wt_annotation.items():
			protein = Protein.objects.get(entry_name=entry_name)
			prot_conf = ProteinConformation.objects.get(protein=protein)
			residues = Residue.objects.filter(protein_conformation=prot_conf, protein_segment__slug__in=['N-term', 'TM1'])
			annotation_dict = OrderedDict()
			anomalies_to_add = []
			for anom, has_anom in self.anomalies[entry_name].items():
				if has_anom!='-':
					anomalies_to_add.append(anom)
			for i, seg_obj in enumerate(self.segments):
				seg = seg_obj.slug
				add_anom = False
				anom_shift = False
				for a in anomalies_to_add:
					if a.startswith(seg):
						add_anom = True
				if seg_obj.fully_aligned:
					segment = OrderedDict()
					x = 50
					if val[seg+'x50']=='-':
						annotation_dict[seg] = OrderedDict()
						continue
					for i in range(int(val[seg+'x50']), int(val[seg+'b'])-1, -1):
						segment[i] = seg+'x'+str(x)
						x-=1
					x = 50
					for i in range(int(val[seg+'x50']), int(val[seg+'e'])+1):
						if anom_shift:
							i+=1
							if i==int(val[seg+'e'])+1:
								continue
						if add_anom and int(a.split('x')[1][:2])==x-1:
							segment[i] = seg+'x'+str(x-1)+'1'
							anom_shift = True
							segment[i+1] = seg+'x'+str(x)
							x+=1
						else:
							segment[i] = seg+'x'+str(x)
							x+=1
					key_order = sorted(list(segment.keys()))
					ordered_segment = OrderedDict([(k, segment[k]) for k in key_order])
					annotation_dict[seg] = ordered_segment
				else:
					segment = OrderedDict()
					if seg=='h4tm1':
						last_helix_end = int(val[self.segments.filter(category='helix').reverse()[0].slug+'e'])
						TM1b = residues.filter(protein_segment__slug='TM1')[0].sequence_number
						for i in range(last_helix_end+1, TM1b):
							segment[i] = None
					else:
						prev_i = i
						while val[self.segments[prev_i-1].slug+'e']=='-':
							prev_i-=1
							if prev_i==0:
								break
						for j in range(int(val[self.segments[prev_i-1].slug+'e'])+1, int(val[self.segments[i+1].slug+'b'])):
							segment[j] = None
					annotation_dict[seg] = segment
			for seg, val in annotation_dict.items():
				for num, gn in val.items():
					res = residues.get(sequence_number=num)
					if gn:
						gn = self.build_ECD_gn(gn, self.segments_dict[seg])
						res.display_generic_number = gn
						res.generic_number = gn
					res.protein_segment = self.segments_dict[seg]
					res.save()
