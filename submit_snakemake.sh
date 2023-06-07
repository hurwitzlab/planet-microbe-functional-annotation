#!/bin/bash

set -u

#SLURM jobname for the Snakemake job submitting all the other jobs
JOBNAME="snakemake_name"


# Parse config/cluster.yml for the user's info for SLURM
. ./bash/parse_yaml.sh
eval $(parse_yaml ./config/cluster.yml "config_")
tmp=$(sed -e "s/'//" -e "s/'//" <<<"$config___default___partition")
PARTITION="${tmp%%[[:space:]]*}"
tmp=$(sed -e "s/'//" -e "s/'//" <<<"$config___default___group")
ACCT="${tmp%%[[:space:]]*}"
ARGS="--partition=${PARTITION} --account=${ACCT}"


# Make error and output directories
export STDERR_DIR="./err"
export STDOUT_DIR="./out"

if [ ! -d "$STDERR_DIR" ]; then
    mkdir "$STDERR_DIR"
fi
if [ ! -d "$STDOUT_DIR" ]; then
    mkdir "$STDOUT_DIR"
fi


JOB_ID=`sbatch $ARGS --job-name=${JOBNAME} -e $STDERR_DIR/${JOBNAME}.err -o $STDOUT_DIR/${JOBNAME}.out ./run_snakemake.sh`
if [ "${JOB_ID}x" != "x" ]; then
    echo Job: \"$JOB_ID\"
else
    echo Problem submitting job. Job terminated.
    exit 1
fi
echo "job successfully submitted"
