#!/bin/bash

set -u
#module load singularity

BOWTIE_IDX="human+phiX"
cd /xdisk/bhurwitz/mig2020/rsgrps/bhurwitz/mattmiller899/planet_microbe/planet-microbe-functional-annotation/data/bowtie_index
INPUT="$1"
OUTPUT_DIR="$2"
THREADS=4
NO_EXT=${INPUT%%.*} 
OUTPUT=${OUTPUT_DIR}/$(basename $NO_EXT).fastq.gz
echo "$OUTPUT"
singularity exec /xdisk/bhurwitz/mig2020/rsgrps/bhurwitz/mattmiller899/planet_microbe/planet-microbe-functional-annotation/singularity/bowtie.simg /bowtie2/bowtie2-2.4.2/bowtie2 -x $BOWTIE_IDX -U $INPUT --un-gz $OUTPUT -p $THREADS
