from Bio import SeqIO
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from collections import OrderedDict
import os
import yaml
import pprint
import math
import matplotlib.pyplot as plt
import numpy as np
from datetime import datetime
import requests


class NewSeqToAlignment():
	def __init__(self, path_to_master_file, path_to_master_file2=None):
		self.path_to_master_file = path_to_master_file
		self.aligned_seq_len = 0
		self.aligned_seqs, self.seqs, self.master_uniprots = self.parse_input(path_to_master_file)
		self.query_seqs = OrderedDict()
		if path_to_master_file2:
			self.aligned_seqs2, self.seqs2, self.master_uniprots2 = self.parse_input(path_to_master_file2)

	# def run_add(self, path_to_file, path_to_query):
	# 	# self.parse_input(path_to_file, path_to_query)
	# 	for q_id, q_seq in self.query_seqs.items():
	# 		ref_id, ref_score, ref_alignment = self.find_closest(q_id, q_seq)
	# 		print(ref_id, ref_score, q_id)
	# 		self.align_to_ref(ref_id, q_id, q_seq)
	# 		break

	def add_to_master(self, entries_to_add):
		with open(self.path_to_master_file, 'a') as f:
			for i, j in entries_to_add.items():
				f.write('>{}\n{}\n'.format(i,j))
			

	def run_edit(self, path_to_file, ref_id, cutoff=0.5):
		list_to_edit = []
		query_align, query_no_gaps = self.parse_input(path_to_file)
		for key, seq in query_align.items():
			if ref_id in key:
				ref_key = key
				ref_seq = seq
				break
		refseq_len = len(ref_seq.replace('-',''))
		ref_gaps = len(ref_seq)-refseq_len
		c=0
		for i, j in query_align.items():
			if len(j.replace('-',''))-ref_gaps<refseq_len*cutoff:
				if '/' in i:
					i = i.split('/')[0]
				list_to_edit.append(i)
				c+=1
		return list_to_edit

	def parse_input(self, path_to_file, check_seqlen=False):
		aligned_seqs, seqs = OrderedDict(), OrderedDict()
		for i, val in enumerate(SeqIO.parse(path_to_file, "fasta")):
			if check_seqlen:
				if i>0 and len(val.seq)!=self.aligned_seq_len:
					print('Warning: aligned sequence length difference: {} and {}'.format(self.aligned_seq_len, len(val.seq)))
			self.aligned_seq_len = len(val.seq)
			aligned_seqs[val.id] = str(val.seq)
		for i, j in aligned_seqs.items():
			seqs[i] = j.replace('-','')
		uniprots = self.get_uniprots(seqs)
		return aligned_seqs, seqs, uniprots

	def get_uniprots(self, seqs):
		uniprot_list = {}
		for i, j in seqs.items():
			if '|' in i:
				split1 = i.split('|')
				if len(split1[0])<6 and len(split1[1])>=6:
					uni_key = split1[1]
					uniprot_list[uni_key] = i
			elif '.' in i:
				split2 = i.split('.')
				if len(split2[0])>=6:
					uni_key = split2[0]
					uniprot_list[uni_key] = i
			else:
				uniprot_list[i] = i
		if len(seqs)!=len(uniprot_list):
			raise AssertionError('Error: Parsing issues with UniProt accessions in master file')
		return uniprot_list
		
	def find_closest(self, query_id, query_seq):
		best_score = 0
		best_id = None
		best_alignment = None
		for ref_id, ref_seq in self.seqs.items():
			pw = pairwise2.align.localms(ref_seq, query_seq, 3, 1, -3, -.1)
			score = pw[0][2]
			if score>best_score:
				best_score = score
				best_id = ref_id
				best_alignment = pw
		return best_id, best_score, best_alignment
		
	def remove_from_master(self, keys_to_remove):
		out = OrderedDict()
		for i, j in self.aligned_seqs.items():
			if i not in keys_to_remove:
				out[i] = j
		return out

	def write_to_file(self, out_file, dict):
		with open(out_file, 'w') as f:
			for i, j in dict.items():
				f.write('>{}\n{}\n'.format(i,j))

	def align_to_ref(self, ref_seq, query_seq, ident_score=4, sim_score=2, gap_open=-2, gap_ext=-.5, verbose=False):
		pw = pairwise2.align.localms(ref_seq, query_seq, ident_score, sim_score, gap_open, gap_ext)
		score = pw[0][2]
		if verbose:
			print(format_alignment(*pw[0]))
			print(score)
			print(self.aligned_seq_len, len(pw[0][1]))
		return score

	def check_content_in_master(self, file_to_compare):
		new_set_with_gaps, new_set = self.parse_input(file_to_compare)
		missing_from_master = OrderedDict()
		c=0
		for i, j in new_set.items():
			orig_key = i
			if '.' in i:
				i = i.split('.')[0]
			if i not in self.master_uniprots:
				missing_from_master[orig_key] = j
				c+=1
		print('Found {} entries missing from master'.format(c))
		return missing_from_master

	def filter(self, input_f):
		if type(input_f)==type('') and os.path.exists(input_f):
			f_aligned, f_seqs, f_uniprots = self.parse_input(input_f)
		else:
			f_uniprots = input_f
		print(f_uniprots)
		o1, o2 = OrderedDict(), OrderedDict()
		for i, j in f_uniprots.items():
			if i in self.master_uniprots:
				o1[i] = self.seqs[self.master_uniprots[i]]
			elif i in self.master_uniprots2:
				o2[i] = self.seqs2[self.master_uniprots2[i]]
		return o1, o2

	def ref_sim_matrix(self, refs=[], files=[], out_file=None):
		sim_matrix = {'References': refs}
		ref_seqs = []
		parsed_files = []
		for f in files:
			if f==self.path_to_master_file:
				for r in refs:
					if r in self.master_uniprots:
						if self.master_uniprots[r] in self.seqs:
							ref_seqs.append(self.seqs[self.master_uniprots[r]])
				parsed_files.append([self.aligned_seqs, self.seqs, self.master_uniprots])
			else:
				parsed_aligned, parsed_seqs, parsed_uniprots = self.parse_input(f)
				for r in refs:
					if r in parsed_uniprots:
						if parsed_uniprots[r] in parsed_seqs:
							ref_seqs.append(parsed_seqs[parsed_uniprots[r]])
				parsed_files.append([parsed_aligned, parsed_seqs, parsed_uniprots])
		if len(refs)==len(ref_seqs):
			print('{} Reference sequences found'.format(len(ref_seqs)))
		
		for f in parsed_files:
			for key, seq in f[1].items():
				sim_matrix[key] = []
				for i, ref in enumerate(ref_seqs):
					if refs[i] in f[2] and f[2][refs[i]]==key:
						continue
					score = self.align_to_ref(ref, seq)
					sim_matrix[key].append(score)
				mean = sum(sim_matrix[key])/len(sim_matrix[key])
				variance = sum([((x - mean) ** 2) for x in sim_matrix[key]]) / len(sim_matrix[key]) 
				res = variance ** 0.5
				sim_matrix[key].append(mean)
				sim_matrix[key].append(res)
		if out_file:
			with open(out_file, 'w') as f:
				yaml.dump(sim_matrix, f, indent=4, default_flow_style=False)
			# with open(out_file, 'r') as f:
			# 	yf = yaml.load(f, Loader=yaml.FullLoader)
			# pprint.pprint(yf)
		# print(sim_matrix)
		return sim_matrix

	def sim_matrix_two_sets(self, out_file, one_file=False):
		startTime = datetime.now()
		print(startTime, 'Running...')
		length = len(self.seqs)
		progress_percentages = [1, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90]
		c=1
		if one_file:
			with open(out_file, 'w') as f:
				for i, j in self.seqs.items():
					if round(c/length*100) in progress_percentages:
						print(datetime.now() - startTime, 'Progress: {}%'.format(round(c/length*100)))
						if len(progress_percentages)>1:
							progress_percentages = progress_percentages[1:]
						else:
							progress_percentages = []
						for k, l in self.seqs2.items():
							score = self.align_to_ref(j, l)
							f.write('{},{},{}\n'.format(i, k, score))
					c+=1
		else:
			for i, j in self.seqs.items():
				if round(c/length*100) in progress_percentages:
					print(datetime.now() - startTime, 'Progress: {}%'.format(round(c/length*100)))
					if len(progress_percentages)>1:
						progress_percentages = progress_percentages[1:]
					else:
						progress_percentages = []
				with open('./classd/sim_matrix/{}_to_ste3.csv'.format(i), 'w') as f:
					for k, l in self.seqs2.items():
						score = self.align_to_ref(j, l)
						f.write('{},{},{}\n'.format(i, k, score))
				c+=1
						

	def process_sim_matrix(self, sim_matrix, score_cutoff=None, sd_cutoff=None):
		if os.path.exists(sim_matrix):
			with open(sim_matrix, 'r') as f:
				sim_matrix = yaml.load(f, Loader=yaml.FullLoader)
		data = []
		for i, j in sim_matrix.items():
			if i=='References':
				continue
			else:
				data.append([i, round(j[-2]), round(j[-1])]+j[:-2])
		sorted_data = sorted(data, key=lambda x: (-x[1], x[2]))
		filtered_data = OrderedDict()
		for i in sorted_data:
			# print(i[0],i[1],i[2],i[3])
			if score_cutoff and sd_cutoff:
				if i[1]>=score_cutoff and i[2]<=sd_cutoff:
					try:
						filtered_data[i[0]] = self.seqs[i[0]]
					except:
						filtered_data[i[0]] = self.seqs2[i[0]]
			else:
				try:
					filtered_data[i[0]] = self.seqs[i[0]]
				except:
					filtered_data[i[0]] = self.seqs2[i[0]]
		return filtered_data

	def remove_duplicates(self, in_file):
		out = {}
		with open(in_file, 'r') as f:
			lines = f.readlines()
			for l in lines:
				if l not in out:
					out[l.replace('\n','')] = ''
		return out

	def check_on_uniprot(self, uniprots):
		for i in uniprots:
			x = requests.get('https://www.uniprot.org/uniprot/{}.xml'.format(i))
			lines = x.text.split('\n')
			try:
				if 'Swiss' in lines[2]:
					print(i, lines[2])
			except:
				print(i, 'Error')

			

		



###STE2
# nsta = NewSeqToAlignment('./classd/Uniprot_STE2_IPR000366_yeast.fasta')

### Add to master
# nsta.run_add('./classd/Focused_STE2.fa','./classd/STE2_to_add.fa')

### Remove based on occupancy cutoff
# list_to_edit = nsta.run_edit('./classd/ste2_occupancy.fasta', 'D6VTK4', cutoff=0.5)
# print('Removed {} entries'.format(len(list_to_edit)))
# print(list_to_edit)
# out = nsta.remove_from_master(list_to_edit)
# nsta.write_to_file('./classd/ste2_occupancy_50_cutoff.fasta', out)


###STE3
# nsta = NewSeqToAlignment('./classd/Uniprot_STE3_IPR001499_yeast.fasta')

### Compare and add to master
# missing = nsta.check_content_in_master('./classd/A0A5C3N0L3_blast.txt')
# nsta.add_to_master(missing)

### Remove based on occupany cutoff
# list_to_edit = nsta.run_edit('./classd/ste3_occupancy_all.fasta', '', cutoff=0.5)
# print('Removed {} entries'.format(len(list_to_edit)))
# print(list_to_edit)
# out = nsta.remove_from_master(list_to_edit)
# nsta.write_to_file('./classd/ste3_occupancy_50_cutoff.fasta', out)

### Combining
# nsta = NewSeqToAlignment('./classd/Uniprot_STE2_IPR000366_yeast.fasta','./classd/Uniprot_STE3_IPR001499_yeast.fasta')
# nsta.ref_sim_matrix(['D6VTK4','P06783'],['./classd/Uniprot_STE2_IPR000366_yeast.fasta','./classd/Uniprot_STE3_IPR001499_yeast.fasta'], './classd/sim_matrix.yaml')
# filtered = nsta.process_sim_matrix('./classd/sim_matrix.yaml', 1000, 50)
# nsta.write_to_file('./classd/filtered_1000_50.fasta',filtered)
# nsta.sim_matrix_two_sets('./classd/all_to_all.csv')
# o = nsta.filter('./classd/combined_smallest_full_tree.fasta')
# nsta.write_to_file('./classd/ste2_filtered.fasta', o[0])
# nsta.write_to_file('./classd/ste3_filtered.fasta', o[1])

### Combining filtered
# nsta = NewSeqToAlignment('./classd/ste2_filtered.fasta', './classd/ste3_filtered.fasta')
# nsta.sim_matrix_two_sets('./classd/all_to_all.csv', one_file=True)
# ste2 = nsta.remove_duplicates('./classd/ste2_filtered_v1.txt')
# ste3 = nsta.remove_duplicates('./classd/ste3_filtered_v1.txt')
# ste2_filter = nsta.filter(ste2)[0]
# ste3_filter = nsta.filter(ste3)[1]
# nsta.write_to_file('./classd/ste2+ste3_filtered.fasta', ste2_filter)
# nsta.write_to_file('./classd/ste2+ste3_filtered2.fasta', ste3_filter)

# nsta = NewSeqToAlignment('./classd/ste2+ste3_109and108_aligned_v2.fasta')
# nsta.write_to_file('./classd/ste2+ste3_109and108_v2.fasta', nsta.seqs)

### Check on UniProt
nsta = NewSeqToAlignment('./classd/ste2+ste3_109and108_v2.fasta')
nsta.check_on_uniprot(nsta.master_uniprots)