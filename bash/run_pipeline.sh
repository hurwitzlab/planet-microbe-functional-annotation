#!/bin/bash

set -u

python ../pipeline/pipeline.py -c $1 -i $2 -o $3
#python /groups/bhurwitz/planet-microbe-functional-annotation/pipeline/pipeline.py -c $1 -i $2 -o $3 
