#!/bin/bash
set -u

OUT_DIR="./out_bowtie"
ERR_DIR="./err_bowtie"
#init_dir "$OUT_DIR" "$ERR_DIR"

ARGS="-q standard -W group_list=bhurwitz -M mattmiller899@email.arizona.edu -m a"
#ARGS="--partition=standard --account=bhurwitz --mail-user=mattmiller899@email.arizona.edu --mail-type=FAIL"
JOB_ID=`qsub $ARGS -N annot_pipeline -e $ERR_DIR -o $OUT_DIR ./run_annot_pipeline.sh`
#JOB_ID=`sbatch $ARGS --job-name=annot_pipeline -e $ERR_DIR/bowtie.err -o $OUT_DIR/bowtie.out ./run_bowtie_mapping.sh`
if [ "${JOB_ID}x" != "x" ]; then
    echo Job: \"$JOB_ID\"
else
     echo Problem submitting job. Job terminated.
     exit 1
fi
echo "job successfully submitted"
