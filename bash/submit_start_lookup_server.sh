#!/bin/bash
set -u

#cd /xdisk/bhurwitz/mig2020/rsgrps/bhurwitz/planet-microbe-functional-annotation
cd $SLURM_SUBMIT_DIR
OUT_DIR="./out"
ERR_DIR="./err"

# Parse config/cluster.yml for the user's info for SLURM
. ./bash/parse_yaml.sh
eval $(parse_yaml ./config/cluster.yml "config_")
tmp=$(sed -e "s/'//" -e "s/'//" <<<"$config___default___partition")
PARTITION="${tmp%%[[:space:]]*}"
tmp=$(sed -e "s/'//" -e "s/'//" <<<"$config___default___group")
ACCT="${tmp%%[[:space:]]*}"
ARGS="--partition=${PARTITION} --account=${ACCT}"

#JOB_ID=`qsub $ARGS -N annot_pipeline -e $ERR_DIR -o $OUT_DIR ./run_annot_pipeline.sh`
JOB_ID=`sbatch $ARGS --job-name=lookup_server -e $ERR_DIR/server.err -o $OUT_DIR/server.out ./bash/run_start_lookup_server.sh`
if [ "${JOB_ID}x" != "x" ]; then
    echo Job: \"$JOB_ID\"
else
     echo Problem submitting job. Job terminated.
     exit 1
fi
echo "job successfully submitted"
