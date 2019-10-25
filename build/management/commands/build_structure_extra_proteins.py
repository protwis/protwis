from build.management.commands.base_build import Command as BaseBuild
from django.conf import settings

from signprot.models import SignprotComplex
from structure.models import StructureExtraProteins, Structure
from protein.models import Protein, ProteinConformation
from residue.models import Residue

from interaction.models import StructureLigandInteraction

import os
import math
import yaml

class Command(BaseBuild):
	help = "Build StructureExtraProteins data"

	g_prot_dict = {'GNAS2':'Gs', 'GNAL':'Golf', 
				   'GNAI1':'Gi1', 'GNAI2':'Gi2', 'GNAI3':'Gi3', 'GNAT1':'Gt1', 'GNAT2':'Gt2', 'GNAT3':'Gt3', 'GNAZ':'Gz', 'GNAO':'Go',
				   'GNAQ':'Gq', 'GNA11':'G11', 'GNA14':'G14', 'GNA15':'G15',
				   'GNA12':'G12', 'GNA13':'G13'}
	arrestin_dict = {'arrs_mouse':'S-arrestin', 'arrs_bovin':'S-arrestin'}

	with open(os.sep.join([settings.DATA_DIR, 'structure_data','extra_protein_notes.yaml']), 'r') as note_file:
		notes = yaml.load(note_file)

	def add_arguments(self, parser):
		super(Command, self).add_arguments(parser=parser)
		parser.add_argument('--purge', help='Purge existing entries in GPCRdb', default=False, action='store_true')

	def handle(self, *args, **options):
		if options['purge']:
			StructureExtraProteins.objects.all().delete()
		self.build_g_protein_alpha_subunits()
		self.build_from_notes()		

	def build_g_protein_alpha_subunits(self):
		sc = SignprotComplex.objects.all()
		for s in sc:
			sep = StructureExtraProteins()
			sep.wt_protein = s.protein
			sep.structure = s.structure
			sep.protein_conformation = ProteinConformation.objects.get(protein__entry_name=s.structure.pdb_code.index.lower()+'_a')
			sep.display_name = self.g_prot_dict[s.protein.family.name]
			if s.structure.pdb_code.index in self.notes:
				sep.note = self.notes[s.structure.pdb_code.index]['note']
			else:
				sep.note = None
			sep.chain = s.alpha
			sep.category = 'G alpha'
			wt_resis = Residue.objects.filter(protein_conformation__protein=s.protein)
			struct_resis = Residue.objects.filter(protein_conformation=sep.protein_conformation)
			sep.wt_coverage = round(len(struct_resis)/len(wt_resis)*100)
			sep.save()
			s.structure.extra_proteins.add(sep)
			s.structure.save()

	def build_from_notes(self):
		for struct, vals in self.notes.items():
			if 'category' in vals:
				sep = StructureExtraProteins()
				if vals['category']=='G alpha':
					try:
						wt_protein = Protein.objects.get(entry_name=vals['prot'].lower()+'_bovin')
					except Protein.DoesNotExist:
						wt_protein = Protein.objects.get(entry_name=vals['prot'].lower()+'_human')
					sep.display_name = self.g_prot_dict[vals['prot']]
				elif vals['category']=='Arrestin':
					wt_protein = Protein.objects.get(entry_name=vals['prot'].lower())
					sep.display_name = self.arrestin_dict[vals['prot']]

				sep.wt_protein = wt_protein
				sep.structure = Structure.objects.get(pdb_code__index=struct)
				sep.protein_conformation = None
				sep.note = vals['note']
				sep.chain = vals['chain']
				sep.category = vals['category']
				wt_resis = Residue.objects.filter(protein_conformation__protein=wt_protein)
				sep.wt_coverage = round(vals['length']/len(wt_resis)*100)

				sep.save()
				sep.structure.extra_proteins.add(sep)
				sep.structure.save()

			