#!/bin/bash

set -u
export STDERR_DIR="./err"
export STDOUT_DIR="./out"
JOBNAME="snakemake_name"
#init_dir "$STDERR_DIR" "$STDOUT_DIR"
ARGS="--partition=standard --account=bhurwitz"
JOB_ID=`sbatch $ARGS --job-name=${JOBNAME} -e $STDERR_DIR/${JOBNAME}.err -o $STDOUT_DIR/${JOBNAME}.out ./run_snakemake.sh`
if [ "${JOB_ID}x" != "x" ]; then
    echo Job: \"$JOB_ID\"
else
    echo Problem submitting job. Job terminated.
    exit 1
fi
echo "job successfully submitted"
