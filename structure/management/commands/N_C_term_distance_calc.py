from Bio.PDB import PDBParser,PDBList

import os, sys

AA = {'ALA':'A', 'ARG':'R', 'ASN':'N', 'ASP':'D',
     'CYS':'C', 'GLN':'Q', 'GLU':'E', 'GLY':'G',
     'HIS':'H', 'ILE':'I', 'LEU':'L', 'LYS':'K',
     'MET':'M', 'PHE':'F', 'PRO':'P', 'SER':'S',
     'THR':'T', 'TRP':'W', 'TYR':'Y', 'VAL':'V', 
     'YCM':'C', 'CSD':'C', 'TYS':'Y', 'SEP':'S'} #non-standard AAs

class NCTermDistanceCalc():
	def __init__(self):
		args = sys.argv
		if len(args)!=3:
			raise AssertionError("Error: Use two arguments, first as input, second as output file")
		else:
			self.input_file = args[1]
			self.output_file = args[2]
			self.pdb_ids = []

	def run(self):
		self.parse_input_file()
		output = []
		for i in self.pdb_ids:
			dist = self.download_parse_calc_pdb(i)
			output.append([i,dist])
		self.prepare_output(output)

	def parse_input_file(self):
		with open(self.input_file, 'r') as f:
			lines = f.readlines()
		for i in lines:
			self.pdb_ids.append(i.replace('\n',''))

	def download_parse_calc_pdb(self, pdb_id):
		pdbl = PDBList()
		pdbl.retrieve_pdb_file(pdb_id, pdir='./', file_format="pdb")
		pdb = PDBParser()
		struct = pdb.get_structure('pdb', './pdb{}.ent'.format(pdb_id).lower())[0]
		for c in struct:
			AAs = []
			for r in c:
				if r.get_resname() in AA:
					try:
						r['CA']
						AAs.append(r)
					except:
						pass
			dist = self.calculate_distance(AAs[0], AAs[-1])
		os.remove('./pdb{}.ent'.format(pdb_id).lower())
		return dist

	def calculate_distance(self, residue1, residue2):
		v1 = residue1['CA'].get_vector()
		v2 = residue2['CA'].get_vector()
		return (v1-v2).norm()

	def prepare_output(self, output_list):
		with open(self.output_file, 'w') as f:
			for i in output_list:
				f.write('{},{}\n'.format(i[0],i[1]))


nctdc = NCTermDistanceCalc()
nctdc.run()


