#!/bin/bash

#SBATCH --job-name=viral_hunt
#SBATCH --account=bhurwitz
#SBATCH --partition=standard
#SBATCH --ntasks=5
#SBATCH --mem=25gb
#SBATCH --time=24:00:00

source activate pyenv
source activate snakemake

cd /xdisk/bhurwitz/mig2020/rsgrps/bhurwitz/planet-microbe-functional-annotation

#dry run
#snakemake -n

echo "snakemake --cluster "sbatch -A {cluster.group} -p {cluster.partition} -n {cluster.n} -t {cluster.time} -mem={cluster.m}"  --cluster-config config/cluster.yml -j 10 --latency-wait 15"

#run in cluster
snakemake --cluster "sbatch -A {cluster.group} -p {cluster.partition} -n {cluster.n} -t {cluster.time} --mem={cluster.m}"  --cluster-config config/cluster.yml -j 30 --latency-wait 15

