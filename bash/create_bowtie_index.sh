#!/bin/bash


mkdir data/bowtie_index
cd data/bowtie
#wget https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_genomic.fna.gz
#gunzip GRCh38_latest_genomic.fna.gz

source ~/.bashrc
source activate pm_env
bowtie2-build -f GRCh38_latest_genomic.fna,phi-X174.fna human+phiX 
mv human+phiX* ../bowtie_index
