from build.management.commands.base_build import Command as BaseBuild
from django.db.models import Q
from django.conf import settings

from protein.models import Protein, ProteinConformation, ProteinAnomaly, ProteinState, ProteinSegment
from residue.models import Residue, ResidueGenericNumberEquivalent
from residue.functions import dgn, ggn
from structure.models import *
from structure.functions import HSExposureCB, PdbStateIdentifier
from common.alignment import AlignedReferenceTemplate, GProteinAlignment
from common.definitions import *
from common.models import WebLink
from signprot.models import SignprotComplex
import structure.structural_superposition as sp
import structure.assign_generic_numbers_gpcr as as_gn
import structure.homology_models_tests as tests

import Bio.PDB as PDB
from modeller import *
from modeller.automodel import *
from collections import OrderedDict
import os
import subprocess
import shlex
import logging
import pprint
from io import StringIO
import sys
import re
import zipfile
import shutil
import math
from copy import deepcopy
from datetime import datetime, date
import yaml
import traceback
import math
import numpy as np


startTime = datetime.now()


class Command(BaseBuild):
	help = 'Test scripts'

	def add_arguments(self, parser):
		parser.add_argument('-s', help='Run StructureAnomalyRecognition on specific crystal structure', default=False)
		parser.add_argument('-r', help='Specify residue number range. E.g. 10-25', default=False)
		parser.add_argument('-a', help='Run on all structures', default=False, action='store_true')
		parser.add_argument('--verbose', help='Verbose', default=False, action='store_true')

	def handle(self, *args, **options):
		if options['a']:
			structures = Structure.objects.all().exclude(structure_type__slug__startswith='af-').filter(protein_conformation__protein__family__parent__parent__parent__slug='001')
			for s in structures:
				sar = StructureAnomalyRecognition(s)
				sar.run_recog()
		else:
			sar = StructureAnomalyRecognition(options['s'], options['r'], options['verbose'])
			sar.run_recog()

class StructureAnomalyRecognition(object):
	def __init__(self, xtal=False, num_range=False, verbose=False):
		if xtal:
			self.verbose = verbose
			try:
				xtal.pdb_code
				self.structure = xtal
			except:
				self.structure = Structure.objects.get(pdb_code__index=xtal.upper())
			self.parent_prot_conf = ProteinConformation.objects.get(protein=self.structure.protein_conformation.protein.parent)
			io = StringIO(self.structure.pdb_data.pdb)
			self.pdb_struct = PDB.PDBParser(QUIET=True).get_structure(self.structure.pdb_code.index, io)[0]
			self.range = []
			if num_range:
				self.range = [[int(i) for i in num_range.split('-')]]
			else:
				for t in ProteinSegment.objects.filter(proteinfamily='GPCR',category='helix'):
					resis = Residue.objects.filter(protein_conformation__protein=self.structure.protein_conformation.protein.parent, protein_segment=t)
					if len(resis)==0:
						continue
					self.range.append([resis[0].sequence_number, resis.reverse()[0].sequence_number])

	def run_recog(self):
		chain = self.pdb_struct[self.structure.preferred_chain[0]]
		constrictions, bulges, values = OrderedDict(),OrderedDict(),OrderedDict()
		for r in chain:
			skip = True
			for ran in self.range:
				if ran[0]-1<r.get_id()[1]<ran[1]:
					skip=False
			if skip:
				continue
			try:
				ca = r['CA'].get_coord()
				ca2 = chain[r.get_id()[1]+2]['CA'].get_coord()
				ca3 = chain[r.get_id()[1]+3]['CA'].get_coord()
				ca5 = chain[r.get_id()[1]+5]['CA'].get_coord()
				b0 = -1.0*(ca2-ca)
				b1 = ca3-ca2
				b2 = ca5-ca3
				b0xb1 = np.cross(b0,b1)
				b1xb2 = np.cross(b2,b1)
				b0xb1_x_b1xb2 = np.cross(b0xb1, b1xb2)
				y = np.dot(b0xb1_x_b1xb2, b1)*(1.0/np.linalg.norm(b1))
				x = np.dot(b0xb1, b1xb2)
				if self.verbose:
					print(chain[r.get_id()[1]+2].get_id()[1],'-',chain[r.get_id()[1]+3].get_id()[1],np.degrees(np.arctan2(y, x)))
				values[r.get_id()[1]] = np.degrees(np.arctan2(y, x))
			except:
				pass
		for num, val in values.items():
			if abs(val)>150:
				count = 1
				for i in range(1,4):
					try:
						if abs(values[num+i])>150:
							count+=1
						else:
							raise Exception
					except:
						break
				if count==3:
					constrictions[num+2] = [val]
			if abs(val)<100:
				count = 1
				for i in range(1,3):
					try:
						if abs(values[num+i])<100:
							count+=1
						else:
							raise Exception
					except:
						break
				if count==3:
					bulges[num+2] = [val]

		found_c, missed_c, found_b, missed_b = OrderedDict(),OrderedDict(),OrderedDict(),OrderedDict()
		remove_c = []
		db_constrictions = self.structure.protein_anomalies.filter(anomaly_type__slug='constriction')
		db_constrictions_dict = OrderedDict()
		for c in db_constrictions:
			gn = c.generic_number.label
			prev_gn = gn[:-1]+str(int(gn[-1])-1)
			prev_resi = Residue.objects.get(protein_conformation=self.parent_prot_conf,
										    display_generic_number__label=dgn(prev_gn, self.parent_prot_conf))
			db_constrictions_dict[gn] = prev_resi.sequence_number
			for ca2 in constrictions:
				if ca2-2<=prev_resi.sequence_number<=ca2+3:
					if gn not in found_c:
						found_c[gn] = prev_resi.sequence_number
						remove_c.append(ca2)
					else:
						remove_c.append(ca2)
			if gn not in found_c:
				missed_c[gn] = ''
		for r in remove_c:
			del constrictions[r]
		for ca2 in constrictions:
			ca2_found = False
			for key, value in found_c.items():
				if ca2-2<=value<=ca2+3:
					ca2_found = True
					break
			if not ca2_found:
				missed_c[ca2] = ca2

		remove_b = []
		db_bulges = self.structure.protein_anomalies.filter(anomaly_type__slug='bulge')
		db_bulges_dict = OrderedDict()
		for b in db_bulges:
			gn = b.generic_number.label
			try:
				resi = Residue.objects.get(protein_conformation=self.parent_prot_conf,
										   display_generic_number__label=dgn(b.generic_number.label, self.parent_prot_conf))
			except ResidueGenericNumberEquivalent.DoesNotExist:
				print('Warning: {} ResidueGenericNumberEquivalent object missing from db'.format(gn))
				continue

			db_bulges_dict[gn] = resi.sequence_number
			for ca2 in bulges:
				if ca2-2<=resi.sequence_number<=ca2+2:
					if gn not in found_b:
						found_b[gn] = resi.sequence_number
						remove_b.append(ca2)
					else:
						remove_b.append(ca2)
			if gn not in found_b:
				missed_b[gn] = ''
		for r in remove_b:
			del bulges[r]
		for ca2 in bulges:
			ca2_found = False
			for key, value in found_b.items():
				if ca2-2<=value<=ca2+3:
					ca2_found = True
					break
			if not ca2_found:
				missed_b[ca2] = ca2

		print('#################')
		print(self.structure)
		print('DB constrictions: {}'.format(db_constrictions_dict))
		print('DB bulges: {}'.format(db_bulges_dict))
		print('Found constrictions: {}'.format(found_c))
		print('Missed constrictions: {}'.format(missed_c))
		print('Found bulges: {}'.format(found_b))
		print('Missed bulges: {}'.format(missed_b))
		for c in missed_c:
			if type(c)==type(0):
				for i in range(c-3,c+3):
					try:
						print(i,values[i])
					except:
						pass
		for b in missed_b:
			if type(b)==type(0):
				for i in range(b-3,b+3):
					try:
						print(i,values[i])
					except:
						pass
