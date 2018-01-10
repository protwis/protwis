from build.management.commands.base_build import Command as BaseBuild
from django.db.models import Q
from django.conf import settings

from protein.models import Protein, ProteinConformation, ProteinAnomaly, ProteinState, ProteinSegment
from residue.models import Residue
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


startTime = datetime.now()


class Command(BaseBuild):  
	help = 'Test scripts'
	
	def add_arguments(self, parser):
		parser.add_argument('--gns', help="Specifiy generic numbers involved in calculation", default=False, nargs='+')


	def handle(self, *args, **options):
		t = TestStateIdentifier(options['gns'])
		t.run()

class TestStateIdentifierSets(object):
	def __init(self):
		pass

	def run(self):
		tm2 = ['2x39','2x40','2x41','2x42']
		tm6 = ['6x35','6x36','6x37','6x38']
		tm3 = ['3x47','3x46','3x45','3x44']
		tm7 = ['7x53','7x52','7x51','7x50']
		best = 1000
		best_params = []
		for t2 in tm2:
			for t6 in tm6:
				for t3 in tm3:
					for t7 in tm7:
						t = TestStateIdentifier([t2, t6, t3, t7])
						t.run()
						if t.mismatch<best:
							best = t.mismatch
							best_params = [t2, t6, t3, t7]
						print(t2,t6,t3,t7, t.mismatch, t.match, t.exceptions)


class TestStateIdentifier(object):
	def __init__(self, gns):
		self.tm2_gn, self.tm6_gn, self.tm3_gn, self.tm7_gn = '2x39', '6x35', '3x47', '7x53'
		if gns:
			for value in gns:
				if value.startswith('2'):
					self.tm2_gn = value
				elif value.startswith('6'):
					self.tm6_gn = value
				elif value.startswith('3'):
					self.tm3_gn = value
				elif value.startswith('7'):
					self.tm7_gn = value

	def run(self):
		strs = Structure.objects.filter(refined=False)
		self.match, self.mismatch, self.exceptions = 0,0,0
		for s in strs:
			try:
				r_s = Structure.objects.get(pdb_code__index=s.pdb_code.index+'_refined')
				psi = PdbStateIdentifier(r_s, self.tm2_gn, self.tm6_gn, self.tm3_gn, self.tm7_gn)
				psi.run()
				psis = PdbStateIdentifier(s, self.tm2_gn, self.tm6_gn, self.tm3_gn, self.tm7_gn)
				psis.run()
				if psi.state!=psis.state:
					print(s, psis.state.slug, psis.activation_value, psi.state.slug, psi.activation_value)
					self.mismatch+=1
				else:
					self.match+=1
			except:
				# print('Exception: ', s)
				self.exceptions+=1
		hommods = StructureModel.objects.all()
		for h in hommods:
			try:
				psih = PdbStateIdentifier(h, self.tm2_gn, self.tm6_gn, self.tm3_gn, self.tm7_gn)
				psih.run()
				psiss = PdbStateIdentifier(h.main_template, self.tm2_gn, self.tm6_gn, self.tm3_gn, self.tm7_gn)
				psiss.run()
				if psih.state!=psiss.state:
					# print(h, psiss.state.slug, psiss.activation_value, psih.state.slug, psih.activation_value)
					self.mismatch+=1
				else:
					self.match+=1
			except:
				# print('Exception hommod:', h)
				self.exceptions+=1
		# print('match:', self.match, 'mismatch:', self.mismatch, 'exceptions:', self.exceptions)
		return 0