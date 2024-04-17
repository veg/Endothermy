import sys
import re

input_fasta = sys.argv[1]
output_fasta = sys.argv[2]

with open(input_fasta, 'r') as in_f:
	fasta_lines = in_f.readlines()

cleaned_fasta = fasta_lines

for ind, line in enumerate(fasta_lines):
	cleaned_fasta[ind] = line.replace('-', '_')

with open(output_fasta, 'w') as out_f:
	out_f.writelines(cleaned_fasta)
