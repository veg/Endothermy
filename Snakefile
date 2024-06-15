"""
@Description: Fish endothermy examination
@Authors: Alexander G Lucaci, Avery Selberg
@Version: 2024.6
"""

# =============================================================================
# Imports
# =============================================================================

import os
#os.environ['OPENBLAS_NUM_THREADS'] = '1'
import sys
import json
import csv
from pathlib import Path
import glob
from Bio import SeqIO

# =============================================================================
# Declares
# =============================================================================

with open("config.MACSE.json", "r") as input_sc:
  config = json.load(input_sc)
#end with
with open("cluster.json", "r") as input_c:
  cluster = json.load(input_c)
#end with

# =============================================================================
# User Settings
# =============================================================================

BASEDIR = os.getcwd()
print("# We are operating out of base directory:", BASEDIR)
# Which project are we analyzing?
LABEL = config["Label"]
DATA_DIRECTORY = os.path.join(BASEDIR, "data", config["DataDirectory"])
print("# Examining data from:", DATA_DIRECTORY)
EXTENSION = "fasta"
PARTITION_LIST = config["Partition"].split(",")
print("# Testing Partition(s):", PARTITION_LIST)
CANDIDATE_GENES = [os.path.basename(x) for x in glob.glob(os.path.join(DATA_DIRECTORY, "*." + EXTENSION))]
## modified candidate genes, takes in tree files from outdir#
#print("# We will process", len(CANDIDATE_GENES), "candidate gene files")

# Check for Multiples of 3 in Codon alignments
def CheckMSA(msa):
    with open(msa, "r") as fh:
        for n, _record in enumerate(SeqIO.parse(fh, "fasta")):
            _id   = _record.id
            _desc = _record.description
            _seq  = _record.seq
            break
        #end for
    #end with
    if len(str(_seq)) % 3 == 0:
        return True
    else:
        return False
    #end if
#end method

# Subset
#CANDIDATE_GENES = CANDIDATE_GENES[0:20]
#print("# We will process", len(CANDIDATE_GENES), "candidate gene files")
#sys.exit(1)

# =============================================================================
# Software Settings
# =============================================================================

HYPHY    = "/home/agselberg/hyphy-develop/hyphy LIBPATH=/home/agselberg/hyphy-develop/res"
HYPHYMPI = "HYPHYMPI"
HYPHY_ANALYSES = os.path.join(BASEDIR, "hyphy-analyses")
STRIKE_AMBIGS_BF = os.path.join(BASEDIR, "scripts", "strike-ambigs.bf")
OUTDIR = os.path.join(BASEDIR, "results", LABEL + "_MACSE")
#CANDIDATE_GENES = glob.glob(os.path.join(OUTDIR,'*.RD.SA.codons.cln.trim67.fa.raxml.bestTree'))
#print("# We will process", len(CANDIDATE_GENES), "candidate gene files")
#sys.exit(1)

OUTDIR_BUSTED= os.path.join(os.path.join(BASEDIR, "results", LABEL + "_MACSE_ENV"))

# Hyphy-analyses
REMOVE_DUPS_BF = os.path.join(HYPHY_ANALYSES, "remove-duplicates", "remove-duplicates.bf")
BUSTED_PH_BF = os.path.join(HYPHY_ANALYSES, "BUSTED-PH", "BUSTED-PH.bf")
LABEL_TREE_BF = os.path.join(HYPHY_ANALYSES, "LabelTrees", "label-tree.bf")
FILTER_OUTLIERS_BF = os.path.join(HYPHY_ANALYSES, "find-outliers", "find-outliers-slac.bf")

# Other scripts
LABEL_BG = os.path.join(BASEDIR, "scripts", "tree-remaining-bg.py")

# Create output directories
Path(os.path.join(BASEDIR, "results")).mkdir(parents=True, exist_ok=True)
Path(OUTDIR).mkdir(parents=True, exist_ok=True)

# Settings, these can be passed in or set in a config.json type file
PPN = cluster["__default__"]["ppn"] 

# =============================================================================
# rule all
# =============================================================================

rule all:
    input:
        expand(os.path.join(OUTDIR, "{GENE}.fa"), GENE=CANDIDATE_GENES),
        expand(os.path.join(OUTDIR, "{GENE}.codons.trim.fa"), GENE=CANDIDATE_GENES),
        expand(os.path.join(OUTDIR, "{GENE}.codons.fa"), GENE=CANDIDATE_GENES),
        expand(os.path.join(OUTDIR, "{GENE}.codons.cln.fa"), GENE=CANDIDATE_GENES),
        expand(os.path.join(OUTDIR, "{GENE}.SA.codons.cln.fa"), GENE=CANDIDATE_GENES),
        expand(os.path.join(OUTDIR, "{GENE}.RD.SA.codons.cln.fa"), GENE=CANDIDATE_GENES),
        expand(os.path.join(OUTDIR, "{GENE}.RD.SA.codons.cln.trim67.fa"), GENE=CANDIDATE_GENES),
        expand(os.path.join(OUTDIR, "{GENE}.RD.SA.codons.cln.trim67.fa.raxml.bestTree"), GENE=CANDIDATE_GENES),
        expand(os.path.join(OUTDIR, "{GENE}.RD.SA.codons.cln.trim67.fa.raxml.bestTree.labeled_fgOnly_{P}.nwk"), GENE=CANDIDATE_GENES, P=PARTITION_LIST),
        expand(os.path.join(OUTDIR, "{GENE}.RD.SA.codons.cln.trim67.fa.raxml.bestTree.labeled_fg-nuOnly_{P}.nwk"), GENE=CANDIDATE_GENES, P=PARTITION_LIST),
        expand(os.path.join(OUTDIR, "{GENE}.RD.SA.codons.cln.trim67.fa.raxml.bestTree.labeled_{P}.nwk"), GENE=CANDIDATE_GENES, P=PARTITION_LIST),
        expand(os.path.join(OUTDIR_BUSTED, "{GENE}.RD.SA.codons.cln.trim67.fa.BUSTED.json"), GENE=CANDIDATE_GENES),
        #expand(os.path.join(OUTDIR_BUSTED, "{GENE}.RD.SA.codons.cln.trim67.fa.BUSTED_{P}.json"), GENE=CANDIDATE_GENES, P=PARTITION_LIST),
        #expand(os.path.join(OUTDIR_BUSTED, "{GENE}.RD.SA.codons.cln.trim67.fa.BUSTED_2rates.json"), GENE=CANDIDATE_GENES),
        expand(os.path.join(OUTDIR_BUSTED, "{GENE}.RD.SA.codons.cln.trim67.fa.BUSTED-E.json"), GENE=CANDIDATE_GENES),
        #expand(os.path.join(OUTDIR_BUSTED, "{GENE}.RD.SA.codons.cln.trim67.fa.BUSTED-E_2rates.json"), GENE=CANDIDATE_GENES),
        expand(os.path.join(OUTDIR_BUSTED, "{GENE}.RD.SA.codons.cln.trim67.fa.FILTER.json"), GENE=CANDIDATE_GENES),
        expand(os.path.join(OUTDIR_BUSTED, "{GENE}.RD.SA.codons.cln.trim67.fa.FILTER.fasta"), GENE=CANDIDATE_GENES),
        expand(os.path.join(OUTDIR_BUSTED, "{GENE}.RD.SA.codons.cln.trim67.fa.FILTERED.labeled_fgOnly_{P}.nwk"), GENE=CANDIDATE_GENES, P=PARTITION_LIST),
        expand(os.path.join(OUTDIR_BUSTED, "{GENE}.RD.SA.codons.cln.trim67.fa.FILTERED.labeled_fg-nuOnly_{P}.nwk"), GENE=CANDIDATE_GENES, P=PARTITION_LIST),
        expand(os.path.join(OUTDIR_BUSTED, "{GENE}.RD.SA.codons.cln.trim67.fa.FILTERED.labeled_{P}.nwk"), GENE=CANDIDATE_GENES, P=PARTITION_LIST),
        expand(os.path.join(OUTDIR_BUSTED, "{GENE}.RD.SA.codons.cln.trim67.fa.{P}.filtered.BUSTED-PH.json"), GENE=CANDIDATE_GENES, P=PARTITION_LIST),
        #expand(os.path.join(OUTDIR_BUSTED, "{GENE}.RD.SA.codons.cln.trim67.fa.{P}.BUSTED-PH.json"), GENE=CANDIDATE_GENES, P=PARTITION_LIST),
        #expand(), GENE=CANDIDATE_GENES, P=PARTITION_LIST),
        #expand(), GENE=CANDIDATE_GENES),

# =============================================================================
# Individual rules
# =============================================================================

rule clean:
    input:
        input = os.path.join(DATA_DIRECTORY, "{GENE}")
    output:
        output = os.path.join(OUTDIR, "{GENE}.fa")
    shell:
       "python scripts/clean-fasta.py {input.input} {output.output}"
# end rule

"""
for sequences in E*unaligned.fasta;
do
    if [  ! -e ${sequences%.*.*}.macse.txt  ];
    then
    echo macse alignment started > ${sequences%.*.*}.macse.txt;
    java -jar $macsepath -prog trimNonHomologousFragments -seq $sequences -min_homology_to_keep_seq 0.4 -out_NT ${sequences%.*}.NT_trimNonHomologousFragments.fasta -out_AA  ${sequences%.*}.AA_trimNonHomologousFragments.fasta;
    java -jar $macsepath -prog alignSequences -seq ${sequences%.*}.NT_trimNonHomologousFragments.fasta -out_NT ${sequences%.*.*}.NT_aligned.fasta -out_AA ${sequences%.*.*}.AA_aligned.fasta;
# replace '!' insertion character with 'N', since it seems to interfere with other software
    sed -i 's/!/N/g' ${sequences%.*.*}.NT_aligned.fasta;
    fi;
done
"""

rule macse_trim:
    input:
        input = rules.clean.output.output
    output:
        codons = os.path.join(OUTDIR, "{GENE}.codons.trim.fa"),
        aa     = os.path.join(OUTDIR, "{GENE}.aa.trim.fa")
    shell:
        "macse -prog trimNonHomologousFragments -seq {input.input} -min_homology_to_keep_seq 0.4 -out_NT {output.codons} -out_AA {output.aa}"
# end rule

rule macse:
    input:
        input = rules.macse_trim.output.codons
    output:
        codons = os.path.join(OUTDIR, "{GENE}.codons.fa"),
        aa     = os.path.join(OUTDIR, "{GENE}.aa.fa")
    shell:
        """
        macse -prog alignSequences -seq {input.input} -out_NT {output.codons} -out_AA {output.aa} -max_refine_iter 3 -local_realign_init 0.3 -local_realign_dec 0.2
        sed -i 's/!/N/g' {output.codons}
        """
# end rule

# Get rid of internal stop codons
# hyphy cln Universal /Users/sweaver/Programming/hyphy/HIV_RT.nex "No/No" /Users/sweaver/Programming/hyphy/HIV_RT_cleaned.nex
rule cln:
   input:
       input = rules.macse.output.codons
   output:
       output = os.path.join(OUTDIR, "{GENE}.codons.cln.fa")
   shell:
       "{HYPHY} CLN Universal {input.input} 'No/No' {output.output}"
# end rule

rule strike_ambigs_msa:
   input:
       input_msa = rules.cln.output.output
   output:
       output = os.path.join(OUTDIR, "{GENE}.SA.codons.cln.fa")
   shell:
      "{HYPHY} {STRIKE_AMBIGS_BF} --alignment {input.input_msa} --output {output.output}"
# end rule

rule remove_duplicates_msa:
   input:
       input_msa = rules.strike_ambigs_msa.output.output
   output:
       output = os.path.join(OUTDIR, "{GENE}.RD.SA.codons.cln.fa")
   shell:
      "{HYPHY} {REMOVE_DUPS_BF} --msa {input.input_msa} --output {output.output} ENV='DATA_FILE_PRINT_FORMAT=9'"
# end rule

rule trimal:
   input:
       input = rules.remove_duplicates_msa.output.output
   output:
       output = os.path.join(OUTDIR, "{GENE}.RD.SA.codons.cln.trim67.fa")
   shell:
       "/home/agselberg/trimal/source/trimal -in {input.input} -out {output.output} -gt 0.67"
# end rule

rule raxml_ng_msa:
    params:
        THREADS = 1 # note: this will not run if too many threads
    input:
        input = rules.trimal.output.output
    output:
        output = os.path.join(OUTDIR, "{GENE}.RD.SA.codons.cln.trim67.fa.raxml.bestTree")
    shell:
        "raxml-ng --model GTR+G --msa {input.input} --threads {params.THREADS} --seed 2 --redo"
# end rule 

#rule prepare_bg_msa:
#    input:
#        foreground_file = os.path.join(BASEDIR, "data", "Partitions", "{P}_FG.txt"),
#        nuisance_file = os.path.join(BASEDIR, "data", "Partitions", "{P}_Nuisance.txt"),
#        all_species = os.path.join(BASEDIR, "data", "Partitions", "All_Species.txt")
#    output:
#        background_file = os.path.join(BASEDIR, "data", "Partitions", "{P}_BG.txt")
#    run:
#        # --- Read in data
#        with open(input.all_species, "r") as f:
#            all_species = [line for line in f]
#        #end with
#        with open(input.nuisance_file, "r") as f:
#            nu_species = [line for line in f]
#        #end with
#        with open(input.foreground_file) as f:
#            fg_species = [line for line in f]
#        #end with
#        # --- Create BG partition
#        species_no_nu = [item for item in all_species if item not in nu_species]
#        bg_species = [item for item in species_no_nu if item not in fg_species]
#        # --- Write to file
#        with open(output.background_file, "w") as f:
#            for line in bg_species:
#                f.write(line)
#            #end for
#        #end with
# end rule

rule label_tree_fg:
    input:
        tree = rules.raxml_ng_msa.output.output,
        partition_file = os.path.join(BASEDIR, "data", "Partitions", "{P}_FG.txt"),   
    output:
        output = os.path.join(OUTDIR, "{GENE}.RD.SA.codons.cln.trim67.fa.raxml.bestTree.labeled_fgOnly_{P}.nwk")
    shell:
        "{HYPHY} {LABEL_TREE_BF} --tree {input.tree} --list {input.partition_file} --output {output.output} --internal-nodes 'Parsimony' --label FOREGROUND" 
# end rule

rule label_tree_nu:
    input:
        tree = rules.label_tree_fg.output.output,
        partition_file = os.path.join(BASEDIR, "data", "Partitions", "{P}_Nuisance.txt")
    output:
        output = os.path.join(OUTDIR, "{GENE}.RD.SA.codons.cln.trim67.fa.raxml.bestTree.labeled_fg-nuOnly_{P}.nwk")
    shell:
        "{HYPHY} {LABEL_TREE_BF} --tree {input.tree} --list {input.partition_file} --output {output.output} --internal-nodes 'Parsimony' --label NUISANCE" 
# end rule

rule label_tree_bg:
    input:
        tree = rules.label_tree_nu.output.output,
    output:
        output = os.path.join(OUTDIR, "{GENE}.RD.SA.codons.cln.trim67.fa.raxml.bestTree.labeled_{P}.nwk")
    shell:
        "python {LABEL_BG} {input.tree} {output.output}"
# end rule

#---------------------------------------------------------------------------- 
# Run BUSTED WITHOUT ERROR SINK COMPONENT
#---------------------------------------------------------------------------- 

## RUN BUSTED ON ALL BRANCHES##
rule busted:
    input:
        seq = rules.trimal.output.output,
        tree = rules.raxml_ng_msa.output.output
    output:
        json = os.path.join(OUTDIR_BUSTED, "{GENE}.RD.SA.codons.cln.trim67.fa.BUSTED.json"), 
        fits = os.path.join(OUTDIR_BUSTED, "{GENE}.RD.SA.codons.cln.trim67.fa.BUSTED-fit.lf"),
        int_fits = os.path.join(OUTDIR_BUSTED, "{GENE}.RD.SA.codons.cln.trim67.fa.BUSTED-intfit.json")
    shell:
        """
        {HYPHY} BUSTED --alignment {input.seq} --tree {input.tree} --srv Yes --starting-points 10 --error-sink No --save-fit {output.fits} --intermediate-fits {output.int_fits} --output {output.json} ENV='TOLERATE_NUMERICAL_ERRORS=1'
        """
# end rule 

rule busted_fg:
    input:
        seq = rules.trimal.output.output,
        tree = rules.label_tree_bg.output.output
    output:
        json = os.path.join(OUTDIR_BUSTED, "{GENE}.RD.SA.codons.cln.trim67.fa.BUSTED_{P}.json"), 
        fits = os.path.join(OUTDIR_BUSTED, "{GENE}.RD.SA.codons.cln.trim67.fa.BUSTED_{P}-fit.lf")
    shell:
        """
        {HYPHY} BUSTED --alignment {input.seq} --tree {input.tree} --srv Yes --starting-points 10 --error-sink No --branches FOREGROUND --save-fit {output.fits} --output {output.json} ENV='TOLERATE_NUMERICAL_ERRORS=1'
        """
###end rule 

rule busted_2_rates:
    input:
        seq = rules.trimal.output.output,
        tree = rules.raxml_ng_msa.output.output
    output:
        json = os.path.join(OUTDIR_BUSTED, "{GENE}.RD.SA.codons.cln.trim67.fa.BUSTED_2rates.json"), 
        fits = os.path.join(OUTDIR_BUSTED, "{GENE}.RD.SA.codons.cln.trim67.fa.BUSTED-fit_2rates.lf")
    shell:
        """
        {HYPHY} BUSTED --alignment {input.seq} --tree {input.tree} --srv Yes --starting-points 10 --error-sink No --save-fit {output.fits} --rates 2 --syn-rates 2 --output {output.json} ENV='TOLERATE_NUMERICAL_ERRORS=1'
        """

#---------------------------------------------------------------------------- 
# Run BUSTED WITH ERROR SINK COMPONENT
#---------------------------------------------------------------------------- 

## RUN BUSTED-E ON ALL BRANCHES##
rule busted_e:
    input:
        seq = rules.trimal.output.output,
        tree = rules.raxml_ng_msa.output.output,
        int_fits = rules.busted.output.int_fits
    output:
        json = os.path.join(OUTDIR_BUSTED, "{GENE}.RD.SA.codons.cln.trim67.fa.BUSTED-E.json"), 
        fits = os.path.join(OUTDIR_BUSTED, "{GENE}.RD.SA.codons.cln.trim67.fa.BUSTED-E-fit.lf")
    shell:
        """
        {HYPHY} BUSTED --alignment {input.seq} --tree {input.tree} --srv Yes --starting-points 10 --error-sink Yes --save-fit {output.fits} --intermediate-fits {input.int_fits} --output {output.json} ENV='TOLERATE_NUMERICAL_ERRORS=1;'
        """
###end rule 

## RUN BUSTED-E ON ALL BRANCHES##
rule busted_e_2_rates:
    input:
        seq = rules.trimal.output.output,
        tree = rules.raxml_ng_msa.output.output
    output:
        json = os.path.join(OUTDIR_BUSTED, "{GENE}.RD.SA.codons.cln.trim67.fa.BUSTED-E_2rates.json"), 
        fits = os.path.join(OUTDIR_BUSTED, "{GENE}.RD.SA.codons.cln.trim67.fa.BUSTED-E-fit-2rates.lf")
    shell:
        """
        {HYPHY} BUSTED --alignment {input.seq} --tree {input.tree} --srv Yes --starting-points 10 --error-sink Yes --save-fit {output.fits} --rates 2 --syn-rates 2 --output {output.json} ENV='TOLERATE_NUMERICAL_ERRORS=1;'
        """
###end rule 

rule filter:
    input:
        e_json = rules.busted_e.output.json
    output:
        json = os.path.join(OUTDIR_BUSTED, "{GENE}.RD.SA.codons.cln.trim67.fa.FILTER.json"), 
        seq = os.path.join(OUTDIR_BUSTED, "{GENE}.RD.SA.codons.cln.trim67.fa.FILTER.fas")
    shell:
        """
        {HYPHY} error-filter {input.e_json} --output {output.seq} --output-json {output.json};
        """
# end rule 

rule get_filtered_tree:
    input:
        seq = rules.filter.output.seq
    output:
        seq = os.path.join(OUTDIR_BUSTED, "{GENE}.RD.SA.codons.cln.trim67.fa.FILTER.fasta"),
        tree = os.path.join(OUTDIR_BUSTED, "{GENE}.RD.SA.codons.cln.trim67.fa.FILTER.nwk"),
    shell:
        """
        head -n -1 {input.seq} > {output.seq};
        tail -n 1 {input.seq} > {output.tree}
        """

rule label_filtered_fg:
    input:
        tree = rules.get_filtered_tree.output.tree,
        partition_file = os.path.join(BASEDIR, "data", "Partitions", "{P}_FG.txt"),   
    output:
        output = os.path.join(OUTDIR_BUSTED, "{GENE}.RD.SA.codons.cln.trim67.fa.FILTERED.labeled_fgOnly_{P}.nwk")
    shell:
        "{HYPHY} {LABEL_TREE_BF} --tree {input.tree} --list {input.partition_file} --output {output.output} --internal-nodes 'Parsimony' --label FOREGROUND" 
#end rule

rule label_filtered_nu:
    input:
        tree = rules.label_filtered_fg.output.output,
        partition_file = os.path.join(BASEDIR, "data", "Partitions", "{P}_Nuisance.txt")
    output:
        output = os.path.join(OUTDIR_BUSTED, "{GENE}.RD.SA.codons.cln.trim67.fa.FILTERED.labeled_fg-nuOnly_{P}.nwk")
    shell:
        "{HYPHY} {LABEL_TREE_BF} --tree {input.tree} --list {input.partition_file} --output {output.output} --internal-nodes 'Parsimony' --label NUISANCE" 
#end rule

rule label_filtered_bg:
    input:
        tree = rules.label_filtered_nu.output.output,
    output:
        output = os.path.join(OUTDIR_BUSTED, "{GENE}.RD.SA.codons.cln.trim67.fa.FILTERED.labeled_{P}.nwk")
    shell:
        "python {LABEL_BG} {input.tree} {output.output}"
#end rule

#---------------------------------------------------------------------------- 
# Run BUSTED-PH
#---------------------------------------------------------------------------- 

## RUN BUSTED-PH AFTER FILTERING##
rule busted_ph_filtered:
    input:
        seq = rules.get_filtered_tree.output.seq,
        tree = rules.label_filtered_bg.output.output
    output:
        json = os.path.join(OUTDIR_BUSTED, "{GENE}.RD.SA.codons.cln.trim67.fa.{P}.filtered.BUSTED-PH.json"), 
    shell:
        """
        {HYPHY} {BUSTED_PH_BF} --alignment {input.seq} --tree {input.tree} --srv Yes --starting-points 10 --branches FOREGROUND --output {output.json} --comparison BACKGROUND ENV='TOLERATE_NUMERICAL_ERRORS=1'
        """
#end rule 

## RUN BUSTED-PH WITHOUT FILTERING##
rule busted_ph_unfiltered:
    input:
        seq = rules.trimal.output.output,
        tree = rules.label_tree_bg.output.output
    output:
        json = os.path.join(OUTDIR_BUSTED, "{GENE}.RD.SA.codons.cln.trim67.fa.{P}.BUSTED-PH.json"), 
    shell:
        """
        {HYPHY} {BUSTED_PH_BF} --alignment {input.seq} --tree {input.tree} --srv Yes --starting-points 10 --branches FOREGROUND --comparison BACKGROUND --output {output.json} ENV='TOLERATE_NUMERICAL_ERRORS=1'
        """
#end rule 

#---------------------------------------------------------------------------- 
# End of file
#---------------------------------------------------------------------------- 
