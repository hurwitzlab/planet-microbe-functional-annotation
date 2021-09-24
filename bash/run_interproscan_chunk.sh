#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=5gb
#SBATCH --time=2:00:00

set -u
cd $SLURM_SUBMIT_DIR
${IPS} -appl Pfam -i ${CURR_FILE} -b ${IPS_OUT} -goterms -iprlookup -dra -cpu 1
echo "Done"
