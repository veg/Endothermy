####################################################################
### THIS IS A PYTHON SCRIPT THAT TAKES IN A FASTA FILE AND TRIMS 
### EACH SEQUENCE IF NECESSARY SO THEY ARE IN MULTIPLES OF 3
####################################################################

import sys
import re

input_fasta = sys.argv[1]
output_fasta = sys.argv[2]

with open(input_fasta, 'r') as in_f:
	fasta_lines = in_f.readlines()

species_ind = []

for ind, line in enumerate(fasta_lines):
	if line.startswith('>'):
		species_ind.append(ind)

fasta_dict = {}

for ind, value in enumerate(species_ind):
	seq_list = []
	seq_string = ""
	if value == species_ind[-1]:
		start_ind = value + 1
		seq_list = fasta_lines[start_ind:]
	else:
		start_ind = value + 1 
		end_ind = species_ind[ind + 1]
		seq_list = fasta_lines[start_ind:end_ind]
	seq_string = ''.join(seq_list)
	seq_string = re.sub(r'\s+', '', seq_string)
	fasta_dict[fasta_lines[value]] = seq_string

#print(fasta_dict)
with open(output_fasta, 'w') as out_f:
	for key, value in fasta_dict.items():
		if len(value) % 3 == 0: 
			seq_to_write = value + '\n'
		elif len(value) % 3 == 1:
			seq_to_write = value[:-1] + '\n'
		elif len(value) % 3 == 2:
			seq_to_write = value[:-2] + '\n'
		out_f.write(key)
		out_f.write(seq_to_write)
