from build.management.commands.base_build import Command as BaseBuild
from build.management.commands.build_homology_models_zip import Command as UploadModel
from django.db.models import Q
from django.conf import settings

from protein.models import Protein, ProteinConformation, ProteinAnomaly, ProteinState, ProteinSegment
from residue.models import Residue
from residue.functions import dgn, ggn
from structure.models import *
from structure.functions import HSExposureCB, PdbStateIdentifier, update_template_source, StructureSeqNumOverwrite
from common.alignment import AlignedReferenceTemplate, GProteinAlignment
from common.definitions import *
from common.models import WebLink
from signprot.models import SignprotComplex
import structure.structural_superposition as sp
import structure.assign_generic_numbers_gpcr as as_gn
import structure.homology_models_tests as tests
from structure.signprot_modeling import SignprotModeling 
from structure.homology_modeling_functions import GPCRDBParsingPDB, ImportHomologyModel, Remodeling

import Bio.PDB as PDB
from modeller import *
from modeller.automodel import *
from collections import OrderedDict
import os
import shlex
import logging
import pprint
from io import StringIO, BytesIO
import sys
import re
import zipfile
import shutil
import math
from copy import deepcopy
from datetime import datetime, date
import yaml
import traceback
import subprocess

startTime = datetime.now()

class Command(BaseBuild):  
	help = 'Build automated chimeric GPCR homology models'
	models_with_knots = []

	def add_arguments(self, parser):
		super(Command, self).add_arguments(parser=parser)
		parser.add_argument('--path', help='Path for zipped model files to check', default=False, type=str)

	def handle(self, *args, **options):
		self.remove_models_with_knots()
		return 0
		if options['path']:
			self.path = options['path']
		else:
			self.path = './structure/complex_models_zip/'
		self.c1, self.c2 = 0,0
		self.file_list = [i for i in os.listdir(self.path) if i.endswith('zip')]
		# self.file_list = ['ClassA_mc5r_human-gnai2_human_6D9H_2019-03-22_GPCRdb.zip']

		self.processors = options['proc']
		self.prepare_input(options['proc'], self.file_list)
		
		# print(self.models_with_knots)
		# print(len(self.file_list), len(self.models_with_knots), str(len(self.models_with_knots)/len(self.file_list)*100)+'%')
		datetime.now() - startTime

	def main_func(self, positions, iteration, count, lock):
		processor_id = round(self.processors*positions[0]/len(self.file_list))+1
		i = 0
		while count.value<len(self.file_list):
			i += 1
			with lock:
				zipname = self.file_list[count.value]
				
				count.value +=1
			modelname = zipname.split('.')[0]
			zip_ref = zipfile.ZipFile(self.path+zipname, 'r')
			zip_ref.extractall(self.path+modelname+'/')
			zip_ref.close()
			p = PDB.PDBParser()
			classname, receptor, species_signprot, species, main_temp, date, placeholder = modelname.split('_')
			receptor_obj = Protein.objects.get(entry_name=receptor+'_'+species_signprot.split('-')[0])
			signprot_obj = Protein.objects.get(entry_name=species_signprot.split('-')[1]+'_'+species)
			model = p.get_structure('model', self.path+modelname+'/'+modelname+'.pdb')
			hse = HSExposureCB(model, radius=9.5, check_chain_breaks=True, check_knots=True, receptor=receptor_obj, signprot=signprot_obj, restrict_to_chain=['A','R'])
			# print(zipname)
			# print(hse.remodel_resis)
			if len(hse.remodel_resis)>0:
				print(zipname)
				# self.models_with_knots.append(zipname)
				# print(self.models_with_knots)
			shutil.rmtree(self.path+modelname)

	def remove_models_with_knots(self):
		with open('./complex_knots_HE_hdhe.txt', 'r') as f:
			zips = f.readlines()
			for z in zips:
				z = z[:-1]
				try:
					shutil.move('./structure/complex_models_zip/'+z, './structure/models_with_knots/'+z)
				except:
					pass
