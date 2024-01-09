#!/bin/bash

clear

echo "Version: v0.1 --- "
echo "2023, COMPARE MODELS TO BUSTED-PH"
echo ""

# Set up the pipeline failure expectations.
set -euo pipefail

echo "Initialized --- "

# Uncomment if you want to generate an analysis DAG file.
#snakemake --forceall --dag | dot -Tpdf > RASCL_DAG.pdf

echo "Creating 'logs' directory"
mkdir -p logs

echo "Executing HPC Snakemake command"

# Execute the Snakemake command
snakemake \
      -s Snakefile_filter \
      --cluster-config cluster.json \
      --cluster "qsub -V -l nodes={cluster.nodes}:ppn={cluster.ppn} -q {cluster.name} -l walltime={cluster.walltime} -e logs -o logs" \
      --jobs 100 all \
      --rerun-incomplete \
      --keep-going \
      --reason \

# End Snakemake command

exit 0

# End of file.
