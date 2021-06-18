#!/bin/bash
set -u

OUT_DIR="./out"
ERR_DIR="./err"
cd /xdisk/bhurwitz/mig2020/rsgrps/bhurwitz/planet-microbe-functional-annotation

#ARGS="-q standard -W group_list=bhurwitz -M mattmiller899@email.arizona.edu -m a"
ARGS="--partition=standard --account=bhurwitz"
#JOB_ID=`qsub $ARGS -N annot_pipeline -e $ERR_DIR -o $OUT_DIR ./run_annot_pipeline.sh`
JOB_ID=`sbatch $ARGS --job-name=lookup_server -e $ERR_DIR/server.err -o $OUT_DIR/server.out ./bash/run_start_lookup_server.sh`
if [ "${JOB_ID}x" != "x" ]; then
    echo Job: \"$JOB_ID\"
else
     echo Problem submitting job. Job terminated.
     exit 1
fi
echo "job successfully submitted"
