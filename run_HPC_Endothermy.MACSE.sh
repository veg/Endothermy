#!/bin/bash

clear

echo ""

# Download HyPhy Standalone analyses
git clone https://github.com/veg/hyphy-analyses.git

# Set up the pipeline failure expectations.
#set -euo pipefail

echo "--- Initialized --- "

# Uncomment if you want to generate an analysis DAG file.
#snakemake --forceall --dag | dot -Tpdf > ENDOTHERMY_DAG.pdf

echo "Creating $LOGDIR directory"
LOGDIR="logs_macse"
mkdir -p $LOGDIR

echo "Executing HPC Snakemake command"

# Execute the Snakemake command

snakemake \
      -s Snakefile_env \
      --cluster-config cluster.json \
      --cluster "qsub -V -l nodes={cluster.nodes}:ppn={cluster.ppn} -q {cluster.name} -l walltime={cluster.walltime} -e $LOGDIR -o $LOGDIR" \
      --jobs 100 all \
      --rerun-incomplete \
      --keep-going \
      --reason \
      --rerun-incomplete \
      --latency-wait 60 


# End Snakemake command

exit 0

# End of file.
