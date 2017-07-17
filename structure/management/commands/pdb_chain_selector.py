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
		print(options['s'], preferred_chain)
