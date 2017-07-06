from build.management.commands.base_build import Command as BaseBuild

from protein.models import Protein
from residue.models import Residue
from structure.functions import PdbChainSelector
from structure.models import *

import Bio.PDB as PDB
from collections import OrderedDict


class Command(BaseBuild):
	def add_arguments(self, parser):
		super(Command, self).add_arguments(parser=parser)
		parser.add_argument('-s', help="PDB code of GPCR structures", default=False, type=str)
		parser.add_argument('-r', help="UniProt common name of GPCR", default=False, type=str)

	def handle(self, *args, **options):
		print(options['s'],options['r'])
		pcs = PdbChainSelector(options['s'], options['r'])
		pcs.run_dssp()
		preferred_chain = pcs.select_chain()
		a1 = pcs.get_seqnums_by_secondary_structure("A", "H")
		a2 = pcs.get_seqnums_by_secondary_structure("A", "G")
		a3 = pcs.get_seqnums_by_secondary_structure("A", "I")
		b1 = pcs.get_seqnums_by_secondary_structure("B", "H")
		b2 = pcs.get_seqnums_by_secondary_structure("B", "G")
		b3 = pcs.get_seqnums_by_secondary_structure("B", "I")
		a = a1+a2+a3
		b = b1+b2+b3
		print("not in b")
		for i in a:
			if i not in b:
				print(i)
		print("not in a")
		for j in b:
			if j not in a:
				print(j)
		print(options['s'], preferred_chain)
