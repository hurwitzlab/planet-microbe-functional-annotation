#!/bin/bash

IPS="/xdisk/bhurwitz/mig2020/rsgrps/bhurwitz/mattmiller899/planet_microbe/interproscan-5.46-81.0/interproscan.properties"

JOB=$(squeue -u $USER | grep "lookup" | head -1 | xargs)
if [[ $JOB == "" ]]; then
    echo "server not running, submitting a server"
    sh /xdisk/bhurwitz/mig2020/rsgrps/bhurwitz/mattmiller899/planet_microbe/planet-microbe-functional-annotation/submit_start_lookup_server.sh
    sleep 5
    JOB=$(squeue -u $USER | grep "lookup" | head -1 | xargs)
fi

while [[ ${JOB: -1} == ")" ]]; do
    echo "waiting for server to start..."
    sleep 5
    JOB=$(squeue -u $USER | grep "lookup" | head -1 | xargs)
done

NODE=${JOB: -7}
echo "node $NODE"
sed -i "s/^precalculated\.match\.lookup\.service\.url=.*/precalculated\.match\.lookup\.service\.url=http:\/\/${NODE}:1234/" $IPS
