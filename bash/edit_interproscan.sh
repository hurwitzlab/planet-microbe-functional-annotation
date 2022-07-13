#!/bin/bash

set -u
cd $SLURM_SUBMIT_DIR
#cd /xdisk/bhurwitz/mig2020/rsgrps/bhurwitz/planet-microbe-functional-annotation
#IPS="/groups/bhurwitz/tools/interproscan-*/interproscan.properties"
IPS="./tools/interproscan-*/interproscan.properties"
JOB=$(squeue -u $USER | grep "lookup" | head -1 | xargs)
if [[ $JOB == "" ]]; then
    echo "server not running, submitting a server"
    #sh submit_start_lookup_server.sh 
    sh ./bash/submit_start_lookup_server.sh
    sleep 5
    JOB=$(squeue -u $USER | grep "lookup" | head -1 | xargs)
fi

while [[ ${JOB: -1} == ")" ]]; do
    echo "waiting for server to start..."
    sleep 5
    JOB=$(squeue -u $USER | grep "lookup" | head -1 | xargs)
done

echo "job $JOB"
NODE=$(echo "$JOB" | tr -s " " | cut -d" " -f8)
#NODE=${JOB: -7}
echo "node $NODE"
sed -i "s/^precalculated\.match\.lookup\.service\.url=.*/precalculated\.match\.lookup\.service\.url=http:\/\/${NODE}:1234/" $IPS
touch $1
