#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --mem=20gb
#SBATCH --time=60:00:00



## PBS -l select=1:ncpus=4:mem=24gb
## PBS -l walltime=12:00:00

set -u

#cd /xdisk/bhurwitz/mig2020/rsgrps/bhurwitz/mattmiller899/my_databases/interproscan/lookup_service*
cd $SLURM_SUBMIT_DIR/tools/lookup_service*
java -Xmx8000m -jar server-*.war --headless --port 1234
