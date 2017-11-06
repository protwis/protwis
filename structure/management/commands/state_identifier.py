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

	def handle(self, *args, **options):
		psi = PdbStateIdentifier(options['s'])
		psi.run()
		print(options['s'], psi.activation_value, psi.state)
