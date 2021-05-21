#!/bin/bash
#PBS -l select=1:ncpus=4:mem=24gb
#PBS -l walltime=12:00:00

### SBATCH --nodes=1
### SBATCH --ntasks=8
### SBATCH --mem=40gb
### SBATCH --time=96:00:00

set -u

cd /xdisk/bhurwitz/mig2020/rsgrps/bhurwitz/mattmiller899/planet_microbe/planet-microbe-functional-annotation
module load singularity
singularity exec /xdisk/bhurwitz/mig2020/rsgrps/bhurwitz/mattmiller899/planet_microbe/planet-microbe-functional-annotation/singularity/bowtie.simg /bowtie2/bowtie2-2.4.2/bowtie2 -x $INDEX -U $IN_FILES -S $SAM -p 8 -f --un-gz $OUT
