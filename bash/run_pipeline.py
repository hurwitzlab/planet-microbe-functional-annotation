#!/bin/bash

set -u

python /xdisk/bhurwitz/mig2020/rsgrps/bhurwitz/mattmiller899/planet_microbe/planet-microbe-functional-annotation/pipeline/pipeline.py -c $1 -i $2 -o $3 
