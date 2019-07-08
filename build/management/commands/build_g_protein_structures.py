from django.core.management.base import BaseCommand, CommandError
from build.management.commands.base_build import Command as BaseBuild
from django.conf import settings
from django.db import connection
from django.db import IntegrityError

from protein.models import (Protein, ProteinGProtein,ProteinGProteinPair, ProteinConformation, ProteinState, ProteinFamily, ProteinAlias,
        ProteinSequenceType, Species, Gene, ProteinSource, ProteinSegment)

from residue.models import (ResidueNumberingScheme, ResidueGenericNumber, Residue, ResidueGenericNumberEquivalent)

from signprot.models import SignprotComplex

import re
from Bio import pairwise2
from collections import OrderedDict
import logging
import shlex, subprocess
from io import StringIO
from Bio.PDB import PDBParser,PPBuilder
from Bio import pairwise2
import pprint


AA = {'ALA':'A', 'ARG':'R', 'ASN':'N', 'ASP':'D',
     'CYS':'C', 'GLN':'Q', 'GLU':'E', 'GLY':'G',
     'HIS':'H', 'ILE':'I', 'LEU':'L', 'LYS':'K',
     'MET':'M', 'PHE':'F', 'PRO':'P', 'SER':'S',
     'THR':'T', 'TRP':'W', 'TYR':'Y', 'VAL':'V', 
     'YCM':'C', 'CSD':'C', 'TYS':'Y', 'SEP':'S'} #non-standard AAs


class Command(BaseBuild):
	def add_arguments(self, parser):
		parser.add_argument('--purge', default=False, action='store_true', help='Purge G protein structures from database')

	def handle(self, *args, **options):
		self.options = options
		if self.options['purge']:
			Residue.objects.filter(protein_conformation__protein__entry_name__endswith='_a', protein_conformation__protein__family__parent__parent__name='Alpha').delete()
			ProteinConformation.objects.filter(protein__entry_name__endswith='_a', protein__family__parent__parent__name='Alpha').delete()
			Protein.objects.filter(entry_name__endswith='_a', family__parent__parent__name='Alpha').delete()

		# Building protein and protconf objects for g protein structure in complex
		scs = SignprotComplex.objects.all()
		for sc in scs:
			self.logger.info('Protein, ProteinConformation and Residue build for alpha subunit of {} is building'.format(sc))
			try:
				# Alpha subunit
				try:
					alpha_protein = Protein.objects.get(entry_name=sc.structure.pdb_code.index.lower()+'_a')
				except:
					alpha_protein = Protein()
					alpha_protein.entry_name = sc.structure.pdb_code.index.lower()+'_a'
					alpha_protein.accession = None
					alpha_protein.name = sc.structure.pdb_code.index.lower()+'_a'
					alpha_protein.sequence = sc.protein.sequence
					alpha_protein.family = sc.protein.family
					alpha_protein.parent = sc.protein
					alpha_protein.residue_numbering_scheme = sc.protein.residue_numbering_scheme
					alpha_protein.sequence_type = ProteinSequenceType.objects.get(slug='mod')
					alpha_protein.source = ProteinSource.objects.get(name='OTHER')
					alpha_protein.species = sc.protein.species
					alpha_protein.save()
				try:
					alpha_protconf = ProteinConformation.objects.get(protein__entry_name=sc.structure.pdb_code.index.lower()+'_a')
				except:
					alpha_protconf = ProteinConformation()
					alpha_protconf.protein = alpha_protein
					alpha_protconf.state = ProteinState.objects.get(slug='active')
					alpha_protconf.save()
				pdbp = PDBParser(PERMISSIVE=True, QUIET=True)
				s = pdbp.get_structure('struct', StringIO(sc.structure.pdb_data.pdb))
				chain = s[0][sc.alpha]
				nums = []
				for res in chain:
					try:
						res['CA']
						nums.append(res.get_id()[1])
					except:
						pass
				
				resis = Residue.objects.filter(protein_conformation__protein=sc.protein)
				num_i = 0
				temp_seq2 = ''
				pdb_num_dict = OrderedDict()
				# Create first alignment based on sequence numbers
				for n in nums:
					if sc.structure.pdb_code.index=='6OIJ' and n<30:
						nr = n+6
					else:
						nr = n
					pdb_num_dict[n] = [chain[n], resis.get(sequence_number=nr)]
				# Find mismatches
				mismatches = []
				for n, res in pdb_num_dict.items():
					if AA[res[0].get_resname()]!=res[1].amino_acid:
						mismatches.append(res)

				pdb_lines = sc.structure.pdb_data.pdb.split('\n')
				seqadv = []
				for l in pdb_lines:
					if l.startswith('SEQADV'):
						seqadv.append(l)
				mutations, shifted_mutations = OrderedDict(), OrderedDict()
				# Search for annotated engineered mutations in pdb SEQADV
				for s in seqadv:
					line_search = re.search('SEQADV\s{1}[A-Z\s\d]{4}\s{1}([A-Z]{3})\s{1}([A-Z]{1})\s+(\d+)[\s\S\d]{5}([\s\S\d]{12})([A-Z]{3})\s+(\d+)(\s\S+)',s)
					if line_search!=None:
						if line_search.group(2)==sc.alpha:
							if line_search.group(4).strip()==sc.protein.accession:
								if line_search.group(3)==line_search.group(6):
									mutations[int(line_search.group(3))] = [line_search.group(1), line_search.group(5)]
								else:
									shifted_mutations[int(line_search.group(3))] = [line_search.group(1), line_search.group(5), int(line_search.group(6))]
							else:
								# Exception for 6G79
								if line_search.group(3)!=line_search.group(6) and 'CONFLICT' in line_search.group(7):
									mutations[int(line_search.group(3))] = [line_search.group(1), line_search.group(5)]
								# Exception for 5G53
								if line_search.group(4).strip()!=sc.protein.accession:
									mutations[int(line_search.group(3))] = [line_search.group(1), line_search.group(5)]
				remaining_mismatches = []

				# Check and clear mismatches that are registered in pdb SEQADV as engineered mutation
				for m in mismatches:
					num = m[0].get_id()[1]
					if num in mutations:
						if m[0].get_resname()!=mutations[num][0] and m[1].amino_acid!=AA[mutations[num][1]]:
							remaining_mismatches.append(m)
					elif num in shifted_mutations:
						remaining_mismatches.append(m)
					else:
						remaining_mismatches.append(m)

				### sanity check
				# print(mutations)
				# print(shifted_mutations)
				# print(mismatches)
				# print(remaining_mismatches)
				# pprint.pprint(pdb_num_dict)

				# Mismatches remained possibly to seqnumber shift, making pairwise alignment to try and fix alignment
				if len(remaining_mismatches)>0 and sc.structure.pdb_code.index!='6OIJ':
					ppb = PPBuilder()
					seq = ''
					for pp in ppb.build_peptides(chain, aa_only=False):
						seq += str(pp.get_sequence())
					pw2 = pairwise2.align.localms(sc.protein.sequence, seq, 2, -1, -.5, -.1)
					ref_seq, temp_seq = str(pw2[0][0]), str(pw2[0][1])
					wt_pdb_dict = OrderedDict()
					pdb_wt_dict = OrderedDict()
					j, k = 0, 0
					for i, ref, temp in zip(range(0,len(ref_seq)), ref_seq, temp_seq):
						if ref!='-' and temp!='-':
							wt_pdb_dict[resis[j]] = pdb_num_dict[nums[k]]
							pdb_wt_dict[pdb_num_dict[nums[k]][0]] = resis[j]
							j+=1
							k+=1
						elif ref=='-':
							wt_pdb_dict[i] = pdb_num_dict[nums[k]]
							pdb_wt_dict[pdb_num_dict[nums[k]][0]] = i
							k+=1
						elif temp=='-':
							wt_pdb_dict[resis[j]] = i
							pdb_wt_dict[i] = resis[j]
							j+=1
					for i, r in enumerate(remaining_mismatches):
						# Adjust for shifted residue when residue is a match
						if r[0].get_id()[1]-remaining_mismatches[i-1][0].get_id()[1]>1:
							pdb_num_dict[r[0].get_id()[1]-1][1] = pdb_wt_dict[chain[r[0].get_id()[1]-1]]
						# Adjust for shifted residue when residue is mutated and it's logged in SEQADV
						if r[0].get_id()[1] in shifted_mutations:
							pdb_num_dict[r[0].get_id()[1]][1] = resis.get(sequence_number=shifted_mutations[r[0].get_id()[1]][2])
						# Adjust for shift
						else:
							pdb_num_dict[r[0].get_id()[1]][1] = pdb_wt_dict[r[0]]

				bulked_residues = []
				for key, val in pdb_num_dict.items():
					# print(key, val) # sanity check
					res_obj = Residue()
					res_obj.sequence_number = val[0].get_id()[1]
					res_obj.amino_acid = AA[val[0].get_resname()]
					res_obj.display_generic_number = val[1].display_generic_number
					res_obj.generic_number = val[1].generic_number
					res_obj.protein_conformation = alpha_protconf
					res_obj.protein_segment = val[1].protein_segment
					bulked_residues.append(res_obj)
				Residue.objects.bulk_create(bulked_residues)
				self.logger.info('Protein, ProteinConformation and Residue build for alpha subunit of {} is finished'.format(sc))
			except Exception as msg:
				print('Protein, ProteinConformation and Residue build for alpha subunit of {} has failed'.format(sc))
				print(msg)
				self.logger.info('Protein, ProteinConformation and Residue build for alpha subunit of {} has failed'.format(sc))
