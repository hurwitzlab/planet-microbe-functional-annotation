#!/bin/bash

#SBATCH --job-name=snakemake_king
#SBATCH --account=bhurwitz
#SBATCH --partition=standard
#SBATCH --ntasks=5
#SBATCH --mem=25gb
#SBATCH --time=24:00:00

source ~/.bashrc
source activate pm_env

cd /xdisk/bhurwitz/mig2020/rsgrps/bhurwitz/planet-microbe-functional-annotation

rm results/interproscan.txt results/killed_interproscan.txt
#dry run
#snakemake -n

snakemake --unlock

echo 'snakemake --cluster "sbatch -A {cluster.group} -p {cluster.partition} -n {cluster.n} -t {cluster.time} -N {cluster.N} -mem={cluster.m} -e {cluster.e} -o {cluster.o}"  --cluster-config config/cluster.yml -j 10 --latency-wait 15 --cores 5'

#run in cluster
snakemake --cluster "sbatch -A {cluster.group} -p {cluster.partition} -n {cluster.n} -t {cluster.time} -N {cluster.N} --mem={cluster.m} -e {cluster.e} -o {cluster.o}"  --cluster-config config/cluster.yml -j 30 --latency-wait 15 --cores 5
