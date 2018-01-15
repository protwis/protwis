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
		super(Command, self).add_arguments(parser=parser)
		parser.add_argument('--gns', help="Specifiy generic numbers involved in calculation for TestStateIdentifier", default=False, nargs='+')
		parser.add_argument('--only_xtals', help="Only run TestStateIdentifier on xtals", default=False, action='store_true')
		parser.add_argument('--cutoffs', help="Set inactive and intermediate cutoffs", default=False, nargs='+')


	def handle(self, *args, **options):
		# strs = Structure.objects.filter(refined=False).exclude(protein_conformation__protein__parent__family__parent__parent__parent__slug__in=['002','005'])
		# for s in strs:
		# 	print(s,',',s.state.slug)
		# return 0
		# t = TestStateIdentifierSets(options['only_xtals'])
		# t.run()
		# with open('./structure/cutoff_test_02.csv','w') as fw:
		# with open('./structure/cutoff_test_04.csv','r') as f:
		# 	lines = f.readlines()
		# 	big_l = []
		# 	for line in lines:
		# 		li = line.split(',')
		# 		li[4] = float(li[4])
		# 		li[-1] = float(li[-1])
		# 		big_l.append(li)
		# 	# big_l = sorted(big_l, key=lambda x: (x[4]))
		# 	c = 0
		# 	for l in big_l:
		# 		c+=1
		# 		if c<104:
		# 			print("['{}','{}','{}','{}',{},{}],".format(l[0],l[1],l[2],l[3],l[7],l[8]))
		# self.data = [['2x39','6x35','3x47','7x53',-0.5,8],
		# 		['2x39','6x35','3x44','7x53',-2,6],
		# 		['2x39','6x38','3x47','7x52',-0.5,6],
		# 		['2x39','6x38','3x44','7x52',-2,4],
		# 		['2x39','6x35','3x47','7x53',1,8],
		# 		['2x39','6x37','3x46','7x53',-0.5,6],
		# 		['2x39','6x35','3x45','7x53',-2,6],
		# 		['2x39','6x37','3x46','7x53',-2,6],
		# 		['2x39','6x38','3x47','7x52',1,6],
		# 		['2x39','6x38','3x46','7x52',1,6],
		# 		['2x39','6x37','3x47','7x53',-2,4],
		# 		['2x39','6x35','3x46','7x53',1,8],
		# 		['2x40','6x35','3x47','7x53',1,8],
		# 		['2x39','6x38','3x47','7x53',-0.5,8],
		# 		['2x42','6x35','3x44','7x53',-0.5,8],
		# 		['2x39','6x38','3x45','7x53',-2,6],
		# 		['2x42','6x35','3x45','7x53',1,8],
		# 		['2x42','6x35','3x44','7x53',1,8],
		# 		['2x39','6x35','3x47','7x52',-0.5,6],
		# 		['2x39','6x35','3x45','7x53',-0.5,6],
		# 		['2x39','6x35','3x44','7x53',-0.5,6],
		# 		['2x39','6x36','3x46','7x53',-0.5,6],
		# 		['2x40','6x35','3x45','7x53',-0.5,6],
		# 		['2x40','6x37','3x46','7x53',-0.5,6],
		# 		['2x41','6x37','3x45','7x53',-2,6],
		# 		['2x41','6x37','3x44','7x53',-2,6],
		# 		['2x39','6x35','3x45','7x52',-2,4],
		# 		['2x39','6x35','3x44','7x52',-2,4],
		# 		['2x42','6x37','3x44','7x53',-2,4],
		# 		['2x41','6x37','3x46','7x53',1,8],
		# 		['2x39','6x35','3x45','7x53',1,6],
		# 		['2x39','6x36','3x46','7x53',1,6],
		# 		['2x41','6x38','3x45','7x52',1,6],
		# 		['2x42','6x38','3x46','7x52',1,6],
		# 		['2x39','6x38','3x45','7x53',-0.5,6],
		# 		['2x40','6x35','3x44','7x53',-0.5,6],
		# 		['2x41','6x37','3x45','7x53',-0.5,6],
		# 		['2x41','6x37','3x44','7x53',-0.5,6],
		# 		['2x39','6x35','3x45','7x52',-0.5,4],
		# 		['2x40','6x38','3x45','7x52',-0.5,4],
		# 		['2x42','6x38','3x44','7x52',-0.5,4],
		# 		['2x42','6x38','3x44','7x53',-2,6],
		# 		['2x39','6x35','3x47','7x50',-2,4],
		# 		['2x39','6x35','3x46','7x51',-2,4],
		# 		['2x39','6x38','3x45','7x52',-2,4],
		# 		['2x41','6x38','3x45','7x51',-2,4],
		# 		['2x39','6x38','3x47','7x53',1,8],
		# 		['2x41','6x37','3x47','7x53',1,8],
		# 		['2x42','6x38','3x47','7x53',1,8],
		# 		['2x39','6x35','3x44','7x53',1,6],
		# 		['2x39','6x37','3x46','7x53',1,6],
		# 		['2x40','6x35','3x45','7x53',1,6],
		# 		['2x40','6x35','3x44','7x53',1,6],
		# 		['2x40','6x37','3x46','7x53',1,6],
		# 		['2x39','6x35','3x46','7x50',-0.5,6],
		# 		['2x39','6x38','3x44','7x53',-0.5,6],
		# 		['2x41','6x35','3x44','7x51',-0.5,6],
		# 		['2x42','6x37','3x46','7x53',-0.5,6],
		# 		['2x42','6x38','3x45','7x53',-0.5,6],
		# 		['2x39','6x38','3x44','7x52',-0.5,4],
		# 		['2x39','6x38','3x44','7x53',-2,6],
		# 		['2x39','6x36','3x44','7x53',-2,4],
		# 		['2x39','6x35','3x47','7x52',1,6],
		# 		['2x39','6x35','3x46','7x52',1,6],
		# 		['2x40','6x38','3x47','7x52',1,6],
		# 		['2x41','6x38','3x46','7x51',1,6],
		# 		['2x42','6x35','3x45','7x53',1,6],
		# 		['2x40','6x35','3x44','7x53',-0.5,8],
		# 		['2x41','6x37','3x47','7x53',-0.5,8],
		# 		['2x42','6x38','3x47','7x53',-0.5,8],
		# 		['2x39','6x35','3x47','7x53',-0.5,6],
		# 		['2x39','6x38','3x46','7x52',-0.5,6],
		# 		['2x40','6x38','3x44','7x52',-0.5,6],
		# 		['2x42','6x35','3x44','7x52',-0.5,6],
		# 		['2x42','6x38','3x44','7x53',-0.5,6],
		# 		['2x39','6x35','3x46','7x51',-0.5,4],
		# 		['2x39','6x35','3x44','7x52',-0.5,4],
		# 		['2x39','6x37','3x47','7x53',-0.5,4],
		# 		['2x39','6x35','3x44','7x53',-2,4],
		# 		['2x39','6x35','3x44','7x50',-2,4],
		# 		['2x39','6x38','3x46','7x51',-2,4],
		# 		['2x40','6x37','3x44','7x53',-2,4],
		# 		['2x42','6x35','3x44','7x51',-2,4],
		# 		['2x42','6x37','3x45','7x53',-2,4],
		# 		['2x42','6x38','3x44','7x52',-2,4],
		# 		['2x42','6x38','3x44','7x50',-2,4],
		# 		['2x40','6x35','3x44','7x53',1,8],
		# 		['2x40','6x38','3x47','7x53',1,8],
		# 		['2x41','6x38','3x44','7x52',1,8],
		# 		['2x39','6x35','3x46','7x50',1,6],
		# 		['2x41','6x37','3x46','7x52',1,6],
		# 		['2x41','6x37','3x45','7x53',1,6],
		# 		['2x42','6x37','3x46','7x53',1,6],
		# 		['2x42','6x38','3x45','7x53',1,6],
		# 		['2x39','6x35','3x44','7x52',1,4],
		# 		['2x42','6x35','3x44','7x50',-0.5,6],
		# 		['2x39','6x35','3x47','7x52',-0.5,4],
		# 		['2x39','6x37','3x46','7x53',-0.5,4],
		# 		['2x42','6x37','3x44','7x53',-0.5,4],
		# 		['2x39','6x35','3x45','7x53',-2,4],
		# 		['2x39','6x37','3x46','7x53',-2,4],
		# 		['2x42','6x38','3x46','7x51',-2,4],
		# 		['2x42','6x38','3x45','7x52',-2,4]]

		# self.data = [['2x39','6x35','3x44','7x53',-2,6.5],
		# 			['2x39','6x35','3x47','7x53',-0.5,8.0],
		# 			['2x39','6x35','3x47','7x53',0,8.0],
		# 			['2x39','6x35','3x44','7x53',-2,6.0],
		# 			['2x39','6x35','3x47','7x53',-0.5,8.5],
		# 			['2x39','6x35','3x47','7x53',0,8.5],
		# 			['2x39','6x35','3x47','7x53',-0.5,7.5],
		# 			['2x39','6x35','3x47','7x53',0,7.5],
		# 			['2x39','6x35','3x44','7x53',-2,5.5],
		# 			['2x39','6x38','3x47','7x52',-1,6.0],
		# 			['2x39','6x38','3x47','7x52',-0.5,6.0],
		# 			['2x39','6x38','3x47','7x52',0,6.0],
		# 			['2x39','6x38','3x47','7x52',0.5,6.0],
		# 			['2x39','6x35','3x44','7x52',-2.5,5.0],
		# 			['2x39','6x35','3x44','7x52',-2,5.0],
		# 			['2x39','6x35','3x44','7x53',-1.5,6.5],
		# 			['2x39','6x35','3x45','7x53',-2.5,6.0],
		# 			['2x39','6x35','3x45','7x53',-1.5,6.0],
		# 			['2x39','6x35','3x46','7x53',1,8.5],
		# 			['2x39','6x35','3x46','7x53',1,9.0],
		# 			['2x39','6x35','3x47','7x53',0.5,8.0],
		# 			['2x39','6x37','3x46','7x53',-1.5,5.5],
		# 			['2x39','6x37','3x46','7x53',-1.5,6.0],
		# 			['2x39','6x38','3x44','7x52',-2,4.0],
		# 			['2x39','6x38','3x44','7x52',-2,5.0],
		# 			['2x40','6x35','3x44','7x53',0,7.5],
		# 			['2x42','6x35','3x44','7x53',0,8.5],
		# 			['2x39','6x35','3x44','7x52',-1.5,5.0],
		# 			['2x39','6x35','3x44','7x53',-2,7.0],
		# 			['2x39','6x35','3x44','7x53',-1.5,6.0],
		# 			['2x39','6x35','3x45','7x53',-2.5,5.5],
		# 			['2x39','6x35','3x45','7x53',-2.5,6.5],
		# 			['2x39','6x35','3x45','7x53',-2,6.0],
		# 			['2x39','6x35','3x45','7x53',-1.5,5.5],
		# 			['2x39','6x35','3x45','7x53',-1.5,6.5],
		# 			['2x39','6x35','3x45','7x53',-1,6.0],
		# 			['2x39','6x35','3x46','7x53',1.5,8.5],
		# 			['2x39','6x35','3x46','7x53',1.5,9.0],
		# 			['2x39','6x35','3x47','7x53',0.5,8.5],
		# 			['2x39','6x35','3x47','7x53',1,8.0],
		# 			['2x39','6x36','3x44','7x53',-3,3.0],
		# 			['2x39','6x37','3x46','7x53',-2,5.5],
		# 			['2x39','6x37','3x46','7x53',-2,6.0],
		# 			['2x39','6x37','3x46','7x53',-1,5.5],
		# 			['2x39','6x37','3x46','7x53',-1,6.0],
		# 			['2x39','6x37','3x46','7x53',-0.5,5.5],
		# 			['2x39','6x37','3x46','7x53',-0.5,6.0],
		# 			['2x39','6x37','3x47','7x53',-3,4.0],
		# 			['2x39','6x37','3x47','7x53',-3,4.5],
		# 			['2x39','6x37','3x47','7x53',-2.5,4.0],
		# 			['2x39','6x37','3x47','7x53',-2.5,4.5],
		# 			['2x39','6x38','3x44','7x52',-2,4.5],
		# 			['2x39','6x38','3x44','7x52',-1.5,4.0],
		# 			['2x39','6x38','3x44','7x52',-1.5,5.0],
		# 			['2x39','6x38','3x44','7x52',-1,4.0],
		# 			['2x39','6x38','3x44','7x52',-1,5.0],
		# 			['2x39','6x38','3x45','7x52',-2.5,3.0],
		# 			['2x39','6x38','3x45','7x52',-2.5,3.5],
		# 			['2x39','6x38','3x46','7x52',0.5,6.0],
		# 			['2x39','6x38','3x46','7x52',0.5,6.5],
		# 			['2x39','6x38','3x47','7x52',-1,5.5],
		# 			['2x39','6x38','3x47','7x52',-0.5,5.5],
		# 			['2x39','6x38','3x47','7x52',0,5.5],
		# 			['2x39','6x38','3x47','7x52',0.5,5.5],
		# 			['2x40','6x35','3x45','7x53',-0.5,7.0],
		# 			['2x40','6x38','3x44','7x52',0,5.5],
		# 			['2x40','6x38','3x47','7x53',1.5,8.5],
		# 			['2x41','6x37','3x44','7x53',-2,6.5],
		# 			['2x41','6x38','3x44','7x52',1,7.5],
		# 			['2x41','6x38','3x44','7x52',1.5,7.5],
		# 			['2x41','6x38','3x44','7x52',2,7.5],
		# 			['2x42','6x35','3x44','7x53',-0.5,8.5],
		# 			['2x42','6x35','3x44','7x53',0.5,8.5],
		# 			['2x42','6x38','3x47','7x53',0,8.0]]
		# self.only_xtals = options['only_xtals']
		# self.processors = options['proc']
		# self.prepare_input(options['proc'], self.data)
		t = TestStateIdentifier(options['gns'], options['only_xtals'], float(options['cutoffs'][0]), float(options['cutoffs'][1]))
		t.run()
		print(datetime.now()-startTime)

	def main_func(self, positions, iteration, count, lock):
		processor_id = round(self.processors*positions[0]/len(self.data))+1
		i = 0
		while count.value<len(self.data):
			i += 1
			with lock:
				d = self.data[count.value]
				count.value +=1
			t = TestStateIdentifier([d[0],d[1],d[2],d[3]],self.only_xtals,d[4],d[5])
			t.run()


class TestStateIdentifierSets(object):
	def __init__(self, only_xtals=False):
		self.only_xtals = only_xtals

	def run(self):
		tm2 = ['2x39','2x40','2x41','2x42']
		tm6 = ['6x35','6x36','6x37','6x38']
		tm3 = ['3x47','3x46','3x45','3x44']
		tm7 = ['7x53','7x52','7x51','7x50']
		inact_cutoffs = [1, -0.5, -2]
		inter_cutoffs = [8, 6, 4]
		best = 1000
		best_params = []
		counter = 0
		with open('./structure/cutoff_test_01.csv', 'w') as f:
			for iac in inact_cutoffs:
				for inc in inter_cutoffs:
					for t2 in tm2:
						for t6 in tm6:
							for t3 in tm3:
								for t7 in tm7:
									counter+=1
									t = TestStateIdentifier([t2, t6, t3, t7], self.only_xtals, iac, inc)
									t.run()
									if t.mismatch<best:
										best = t.mismatch
										best_params = [t2, t6, t3, t7, iac, inc]
									print(counter, t2,t6,t3,t7, t.mismatch, t.match, t.exceptions, iac, inc)
									f.write('{},{},{},{},{},{},{},{},{}\n'.format(t2, t6, t3, t7, t.mismatch, t.match, t.exceptions, iac, inc))
		print(best_params, best)


class TestStateIdentifierBestSets(TestStateIdentifierSets):
	def run(self, data):
		cutoff_finetune = [-1,-0.5,0,0.5,1]
		for plus_iac in cutoff_finetune:
			for plus_inc in cutoff_finetune:
				t = TestStateIdentifier([data[0],data[1],data[2],data[3]], self.only_xtals, data[4]+plus_iac, data[5]+plus_inc)
				t.run()
				print('{}-{}-{}-{},{},{},{},{},{}'.format(data[0],data[1],data[2],data[3],t.mismatch,t.match,t.exceptions,data[4]+plus_iac,data[5]+plus_inc))


class TestStateIdentifier(object):
	def __init__(self, gns, only_xtals=False, inact_cutoff=-1, inter_cutoff=8):
		self.tm2_gn, self.tm6_gn, self.tm3_gn, self.tm7_gn = '2x39', '6x35', '3x44', '7x53'
		self.inact_cutoff = inact_cutoff
		self.inter_cutoff = inter_cutoff
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
		self.only_xtals = only_xtals

	def run(self):
		strs = Structure.objects.filter(refined=False).exclude(protein_conformation__protein__parent__family__parent__parent__parent__slug__in=['002','005'])
		self.match, self.mismatch, self.exceptions = 0,0,0
		for s in strs:
			try:
				if self.only_xtals:
					psis = PdbStateIdentifier(s, self.tm2_gn, self.tm6_gn, self.tm3_gn, self.tm7_gn, self.inact_cutoff, self.inter_cutoff)
					psis.run()
					if psis.state!=s.state:
						print(s, s.state, s.distance, psis.state, psis.activation_value, 'mismatch')
						self.mismatch+=1
					else:
						# print(s, s.state, s.distance, psis.activation_value)
						self.match+=1
				else:
					r_s = Structure.objects.get(pdb_code__index=s.pdb_code.index+'_refined')
					psi = PdbStateIdentifier(r_s, self.tm2_gn, self.tm6_gn, self.tm3_gn, self.tm7_gn, self.inact_cutoff, self.inter_cutoff)
					psi.run()
					psis = PdbStateIdentifier(s, self.tm2_gn, self.tm6_gn, self.tm3_gn, self.tm7_gn, self.inact_cutoff, self.inter_cutoff)
					psis.run()
					if psi.state!=psis.state:
						print(s, psis.state.slug, psis.activation_value, psi.state.slug, psi.activation_value)
						self.mismatch+=1
					else:
						self.match+=1
			except:
				# print('Exception: ', s)
				self.exceptions+=1
		if not self.only_xtals:
			hommods = StructureModel.objects.all().exclude(protein__family__parent__parent__parent__slug__in=['002','003','005'])
			for h in hommods:
				try:
					psih = PdbStateIdentifier(h, self.tm2_gn, self.tm6_gn, self.tm3_gn, self.tm7_gn, self.inact_cutoff, self.inter_cutoff)
					psih.run()
					psiss = PdbStateIdentifier(h.main_template, self.tm2_gn, self.tm6_gn, self.tm3_gn, self.tm7_gn, self.inact_cutoff, self.inter_cutoff)
					psiss.run()
					if psih.state!=psiss.state:
						print(h, psiss.state.slug, psiss.activation_value, psih.state.slug, psih.activation_value)
						self.mismatch+=1
					else:
						self.match+=1
				except:
					# print('Exception hommod:', h)
					self.exceptions+=1
		print('match:', self.match, 'mismatch:', self.mismatch, 'exceptions:', self.exceptions)
		print('{}-{}-{}-{},{},{},{},{},{}'.format(self.tm2_gn, self.tm6_gn, self.tm3_gn, self.tm7_gn, self.mismatch, self.match, self.exceptions, self.inact_cutoff, self.inter_cutoff))
		return 0