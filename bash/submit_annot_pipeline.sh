#!/bin/bash
set -u

OUT_DIR="./out_test"
ERR_DIR="./err_test"
#init_dir "$OUT_DIR" "$ERR_DIR"
export TEST_TYPE="small"
#ARGS="-q standard -W group_list=bhurwitz -M mattmiller899@email.arizona.edu -m a"
ARGS="--partition=standard --account=bhurwitz --mail-user=mattmiller899@email.arizona.edu --mail-type=FAIL"
#JOB_ID=`qsub $ARGS -N annot_pipeline -e $ERR_DIR -o $OUT_DIR ./run_annot_pipeline.sh`
JOB_ID=`sbatch $ARGS --export=TEST_TYPE --job-name=annot_pipeline -e $ERR_DIR/${TEST_TYPE}.err -o $OUT_DIR/${TEST_TYPE}.out ./run_annot_pipeline.sh`
if [ "${JOB_ID}x" != "x" ]; then
    echo Job: \"$JOB_ID\"
else
     echo Problem submitting job. Job terminated.
     exit 1
fi
echo "job successfully submitted"
