from build.management.commands.base_build import Command as BaseBuild

from protein.models import Protein
from residue.models import Residue
from structure.functions import PdbStateIdentifier
from structure.models import *

import Bio.PDB as PDB
from collections import OrderedDict


class Command(BaseBuild):
	def add_arguments(self, parser):
		super(Command, self).add_arguments(parser=parser)
		parser.add_argument('-s', help="PDB code of GPCR structures", default=False, type=str)
		parser.add_argument('--state', help="Activation state in case of homology model", default=False, type=str)
		parser.add_argument('--gns', help="Specifiy generic numbers involved in calculation", default=False, nargs='+')
		parser.add_argument('--cutoffs', help="Specify inactive-intermediate and intermediate-active cutoffs", default=False, nargs='+')

	def handle(self, *args, **options):
		try:
			s = Structure.objects.get(pdb_code__index=options['s'])
		except:
			s = StructureModel.objects.get(protein__entry_name=options['s'], state__slug=options['state'])
		if options['gns']:
			if s.protein_conformation.protein.family.slug.startswith('002') or s.protein_conformation.protein.family.slug.startswith('003'):
				tm2_gn, tm6_gn, tm3_gn, tm7_gn = '2x41', '6x33', '3x44', '7x51'
			else:	
				tm2_gn, tm6_gn, tm3_gn, tm7_gn = '2x41', '6x38', '3x44', '7x52'
			for value in options['gns']:
				if value.startswith('2'):
					tm2_gn = value
				elif value.startswith('6'):
					tm6_gn = value
				elif value.startswith('3'):
					tm3_gn = value
				elif value.startswith('7'):
					tm7_gn = value
			if options['cutoffs']:
				psi = PdbStateIdentifier(s, tm2_gn=tm2_gn, tm6_gn=tm6_gn, tm3_gn=tm3_gn, tm7_gn=tm7_gn, inactive_cutoff=float(options['cutoffs'][0]), intermediate_cutoff=float(options['cutoffs'][1]), )
			else:
				psi = PdbStateIdentifier(s, tm2_gn=tm2_gn, tm6_gn=tm6_gn, tm3_gn=tm3_gn, tm7_gn=tm7_gn)
		else:
			psi = PdbStateIdentifier(s)
		psi.run()
		print(options['s'], psi.activation_value, psi.state)
