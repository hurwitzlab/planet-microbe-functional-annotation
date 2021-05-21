#!/bin/bash

set -u
export STDERR_DIR="./err"
export STDOUT_DIR="./out"
#init_dir "$STDERR_DIR" "$STDOUT_DIR"
ARGS="-q standard -W group_list=bhurwitz -M mattmiller899@email.arizona.edu -m a"
export GENOME_FILELIST="/xdisk/bhurwitz/mig2020/rsgrps/bhurwitz/mattmiller899/planet_microbe/planet-microbe-functional-annotation/bowtie_filelist.txt"
export INDEX_DIR="/xdisk/bhurwitz/mig2020/rsgrps/bhurwitz/mattmiller899/planet_microbe/planet-microbe-functional-annotation/data/bowtie_index"
init_dir "$INDEX_DIR"
export INDEX_NAME="human+phiX"
JOB_ID=`qsub $ARGS -v GENOME_FILELIST,INDEX_DIR,INDEX_NAME -N build_bowtie -e "$STDERR_DIR" -o "$STDOUT_DIR" ./run_generate_bowtie_index.sh`
if [ "${JOB_ID}x" != "x" ]; then
    echo "Job: \"$JOB_ID\""
else
    echo Problem submitting job. Job terminated.
    exit 1
fi
echo "job successfully submitted"
