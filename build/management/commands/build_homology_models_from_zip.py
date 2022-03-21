from build.management.commands.base_build import Command as BaseBuild

from protein.models import Protein, ProteinConformation, ProteinState
from residue.models import Residue
from structure.models import *

import zipfile
import os
import shutil

class Command(BaseBuild):  
	help = 'Build GPCRdb homology models from archived zip file'
    
	def add_arguments(self, parser):
		parser.add_argument('-f', help='''Path to zipfile to be parsed''', 
    						default=True, type=str)

	def handle(self, *args, **options):
		path_to_zip = options['f']
		zip_ref = zipfile.ZipFile(path_to_zip, 'r')
		path_to_files = '/'.join(path_to_zip.split('/')[:-1])
		zip_ref.extractall(path_to_files)
		zip_ref.close()
		all_files = os.listdir(path_to_files+'/homology_models')
		for f in all_files:
			if 'post'  in f:
				continue
			if f.endswith('.pdb'):
				p = UploadHomologyModelsToDB(path_to_files+'/homology_models/'+f)
				p.parse_model()
				p.parse_template_stats()
				p.parse_seqsim()
		shutil.rmtree(path_to_files+'/homology_models')


class UploadHomologyModelsToDB():
	''' Parse and upload to local GPCRdb all 3 types of homology model files (.pdb, .templates.csv, .template_similarities.csv) from .zip file.
	'''
	def __init__(self, path_to_model_pdb):
		self.path_to_model_pdb = path_to_model_pdb
		pdb_filename = self.path_to_model_pdb.split('/')[-1]
		self.model_class, self.receptor, self.species, self.state, self.main_template, fileformat = pdb_filename.split('_')
		self.homology_model = None

	def parse_model(self):
		''' Parse .pdb model file
		'''
		with open(self.path_to_model_pdb, 'r') as f:
			lines = f.readlines()
			version = lines[0][-11:-1]
			pdb = '\n'.join(lines)
		self.protein_obj = Protein.objects.get(entry_name=self.receptor+'_'+self.species)
		main_temp_obj = Structure.objects.get(pdb_code__index=self.main_template)
		state_obj = ProteinState.objects.get(name=self.state)
		try:
			StructureModel.objects.filter(protein=self.protein_obj, state=state_obj).delete()
		except:
			pass
		self.homology_model, created = StructureModel.objects.update_or_create(protein=self.protein_obj, state=state_obj, main_template=main_temp_obj,
																			   pdb=pdb, version=version)

	def parse_template_stats(self):
		''' Parse .templates.csv model file
		'''
		try:
			StructureModelStatsRotamer.objects.filter(homology_model=self.homology_model).delete()
		except:
			pass
		with open(self.path_to_model_pdb[:-3]+'templates.csv', 'r') as f:
			lines = f.readlines()
			for line in lines[1:]:
				seg, seqnum, gn, rec, back, rot = line.split(',')
				rot = rot[:-1]
				try:
					b_temp = Structure.objects.get(pdb_code__index=back)
				except:
					b_temp = None
				try:
					r_temp = Structure.objects.get(pdb_code__index=rot)
				except:
					r_temp = None
				res = Residue.objects.get(protein_conformation__protein=self.protein_obj, sequence_number=int(seqnum))
				rota, created = StructureModelStatsRotamer.objects.update_or_create(homology_model=self.homology_model, backbone_template=b_temp, rotamer_template=r_temp,
																					residue=res)

	def parse_seqsim(self):
		''' Parse .template_similarities.csv model file
		'''
		try:
			StructureModelSeqSim.objects.filter(homology_model=self.homology_model).delete()
		except:
			pass
		with open(self.path_to_model_pdb[:-3]+'template_similarities.csv') as f:
			lines = f.readlines()
			for line in lines[1:]:
				temp, sim, resolution, repres, state = line.split(',')
				temp_obj = Structure.objects.get(pdb_code__index=temp)
				sim_obj, created = StructureModelSeqSim.objects.update_or_create(homology_model=self.homology_model, template=temp_obj, similarity=int(sim))
