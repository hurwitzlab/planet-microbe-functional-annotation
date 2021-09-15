#!/bin/bash

set -u
#module load singularity

BOWTIE_IDX="human+phiX"
cd $SLURM_SUBMIT_DIR/data/bowtie_index
INPUT="$1"
OUTPUT_DIR="$2"
THREADS=4
NO_EXT=${INPUT%%.*} 
OUTPUT=${OUTPUT_DIR}/$(basename $NO_EXT).fastq.gz
echo "$OUTPUT"
singularity exec /groups/bhurwitz/planet-microbe-functional-annotation/singularity/bowtie.simg /bowtie2/bowtie2-2.4.2/bowtie2 -x $BOWTIE_IDX -U $INPUT --un-gz $OUTPUT -p $THREADS
