import os 
import sys

output_ext = [
	'.unlagined.nodup.fasta.fa_codons.ID.RD.SA.fasta.FILTER.fas',
	'.unlagined.nodup.fasta.fa_codons.ID.RD.SA.fasta.filtered.BUSTED-PH.json',
	'.unlagined.nodup.fasta.fa_codons.ID.RD.SA.fasta.BUSTED-PH.json',
	'.unlagined.nodup.fasta.fa_codons.ID.RD.SA.fasta.BUSTED.json',
	'.unlagined.nodup.fasta.fa_codons.ID.RD.SA.fasta.BUSTED-E.json'
]
output_complete_dir = '/data/shares/veg/ray_finned/Complete/data/sc_1_busteds'
input_ext = '.unlagined.nodup.fasta.fa_codons.ID.RD.SA.fasta.raxml.bestTree'
input_dir = '/home/agselberg/ray_finned/trees_unlab'


in_files_full = [x for x in os.listdir(input_dir) if x.endswith(input_ext)]
in_files_split = [x.split(input_ext)[0] for x in in_files_full]

out_file_incomplete = []
for ind, o_ext in enumerate(output_ext):

	out_file_complete = [x for x in os.listdir(output_complete_dir) if x.endswith(o_ext)]
	out_complete_split = [x.split(o_ext)[0] for x in out_file_complete]
	out_to_concat = [item for item in in_files_split if item not in out_complete_split]
	#print(o_ext)
	#print(len(out_to_concat))
	out_file_incomplete = list(set(out_file_incomplete + out_to_concat))

#print(len(out_file_incomplete))
for ind, item in enumerate(out_file_incomplete):
	os.system('cp ' + os.path.join(input_dir, item + input_ext) + ' /home/agselberg/ray_finned/results-e_sc1/trees')
	os.system('cp ' + os.path.join('/home/agselberg/ray_finned/cleaned_fastas', item + '.unlagined.nodup.fasta.fa_codons.ID.RD.SA.fasta') + ' /home/agselberg/ray_finned/results-e_sc1/fastas')
	
