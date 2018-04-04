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
		parser.add_argument('--ss', help="""Output number of number of residues with given secondary structure. 
							H: alpha helix, B: isolated beta-bridge, E: strand, G: 3-10 helix, I: pi-helix, T: turn, S: bend, -: Other""",
							default=False, type=str)

	def handle(self, *args, **options):
		print(options['s'],options['r'])
		pcs = PdbChainSelector(options['s'], options['r'])
		pcs.run_dssp()
		preferred_chain = pcs.select_chain()
		print(options['s'], preferred_chain)
		if options['ss']:
			ss = pcs.get_seqnums_by_secondary_structure(preferred_chain, options['ss'])
			print('Secondary structure: {}'.format(options['ss']), ss)
