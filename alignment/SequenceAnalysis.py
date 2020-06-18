from collections import OrderedDict
from Bio import SeqIO
import subprocess, shlex, os


class SequenceAnalysis():
	def __init__(self, alignment_file):
		self.entries = OrderedDict()
		for record in SeqIO.parse(alignment_file, 'fasta'):
			self.entries[record.id] = str(record.seq)

	def calculate_conservation_for_reference(self, reference):
		set_length = len(self.entries)
		for key, seq in self.entries.items():
			if reference in key:
				ref_key = key
				break
		out = OrderedDict()
		for i in range(1, len(self.entries[ref_key])+1):
			out[i] = [self.entries[ref_key][i-1], 0]
		for key, seq in self.entries.items():
			for i, s in enumerate(seq):
				if s==out[i+1][0]:
					out[i+1][1]+=1
		for i, val in out.items():
			print(i, val, round(val[1]/set_length*100))

sa = SequenceAnalysis('./classd/STE2_395_aligned_final.fasta')
sa.calculate_conservation_for_reference('STE2_YEAST')