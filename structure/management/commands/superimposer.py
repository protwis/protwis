from build.management.commands.base_build import Command as BaseBuild

from protein.models import Protein
from residue.models import Residue
from structure.structural_superposition import RotamerSuperpose
from structure.models import *

import Bio.PDB as PDB
from collections import OrderedDict


class Command(BaseBuild):
	def add_arguments(self, parser):
		super(Command, self).add_arguments(parser=parser)
		parser.add_argument('-f', help="Files", nargs='+')
		parser.add_argument('-c', help="Chain IDs", nargs='+')
		parser.add_argument('--r1', help="Reference residue numbers", default=False, nargs='+')
		parser.add_argument('--r2', help="Target residue numbers", default=False, nargs='+')


	def handle(self, *args, **options):
		if len(options['f'])!=2:
			raise AssertionError('Error: Please give two files as arguments for -f, one reference, one target')
		if len(options['c'])!=2:
			raise AssertionError('Error: Please give two chain IDs as arguments for -c, one reference, one target')
		p1 = PDB.PDBParser()
		ref = p1.get_structure('ref',options['f'][0])
		p2 = PDB.PDBParser()
		tar = p2.get_structure('tar',options['f'][1])
		r1_list, a1_list, r2_list, a2_list = [], [], [], []
		for i in options['r1']:
			if '-' in i:
				i_s, i_e = i.split('-')
				r1_list+=range(int(i_s),int(i_e)+1)
			else:
				r1_list.append(int(i))
		for j in options['r2']:
			if '-' in j:
				j_s, j_e = j.split('-')
				r2_list+=range(int(j_s),int(j_e)+1)
			else:
				r2_list.append(int(j))
		for r1 in r1_list:
			for atom1 in ref[0][options['c'][0]][r1]:
				a1_list.append(atom1)
		for r2 in r2_list:
			for atom2 in tar[0][options['c'][1]][r2]:
				a2_list.append(atom2)
		rs = RotamerSuperpose(a1_list, a2_list)
		new_atoms = rs.run()
		for n in new_atoms:
			tar[0][options['c'][1]][n.get_parent().get_full_id()[-1][1]][n.name].coord = n.coord
		io = PDB.PDBIO()
		io.set_structure(tar)
		f_name, extension = options['f'][1].split('.')
		io.save('{}_superimposed.{}'.format(f_name, extension))