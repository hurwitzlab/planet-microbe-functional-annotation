#!/bin/bash

set -u

#OSD illumina example (lower bound)
#curl -O ftp.sra.ebi.ac.uk/vol1/fastq/ERR771/ERR771001/ERR771001_1.fastq.gz
# curl -O ftp.sra.ebi.ac.uk/vol1/fastq/ERR771/ERR771001/ERR771001_2.fastq.gz

# Tara 454 example (test for 454)
#curl -O ftp.sra.ebi.ac.uk/vol1/fastq/ERR321/ERR321017/ERR321017.fastq.gz


# Tara illumina example (upper bound)
#curl -O ftp.sra.ebi.ac.uk/vol1/fastq/ERR599/ERR599053/ERR599053_1.fastq.gz
#curl -O ftp.sra.ebi.ac.uk/vol1/fastq/ERR599/ERR599053/ERR599053_2.fastq.gz



### Round 2 testing Data ###
# Project #Method #Size (Gb) #Notes

## Datasets very likly for inclusion: ##

# Amazon Plume Metagenomes #Illumina #2.9 Gb #largest in dataset
#curl -O ftp.sra.ebi.ac.uk/vol1/fastq/SRR483/004/SRR4831664/SRR4831664.fastq.gz

# Amazon River Metagenomes #Illumina #3.6 Gb #largest in dataset
curl -O ftp.sra.ebi.ac.uk/vol1/fastq/SRR179/006/SRR1796116/SRR1796116_1.fastq.gz

# BATS Chisholm #Illumina #9.8 Gb #largest in dataset
curl -O ftp.sra.ebi.ac.uk/vol1/fastq/SRR650/000/SRR6507280/SRR6507280_1.fastq.gz

# HOT ALOHA time/depth series #Illumina #10.4 Gb #6th largest in dataset largest is 15.5 Gb
curl -O ftp.sra.ebi.ac.uk/vol1/fastq/SRR500/003/SRR5002313/SRR5002313.fastq.gz

# HOT Chisholm #Illumina #6.2 Gb #largest in dataset
curl -O ftp.sra.ebi.ac.uk/vol1/fastq/SRR650/007/SRR6507277/SRR6507277_1.fastq.gz

# Tara Oceans #Illumina #10 Gb # Around the median size, range is 0.0876-35.3 Gb
curl -O ftp.sra.ebi.ac.uk/vol1/fastq/ERR599/ERR599164/ERR599164_1.fastq.gz

## Datasets to maybe include: TBD ##

# GOS #454 #0.7302 Gb #largest in dataset
#curl -O ftp.sra.ebi.ac.uk/vol1/fastq/ERR986/ERR986609/ERR986609.fastq.gz

# HOT DeLong #454 #1.9 Gb #largest in dataset
curl -O ftp.sra.ebi.ac.uk/vol1/fastq/SRR130/001/SRR1303821/SRR1303821.fastq.gz

# HOT DeLong Metatranscriptomes #454 #0.2838 Gb #largest in dataset, our data only includes metagenomes not metatrancriptome, small and probably not worth including
#curl -O ftp.sra.ebi.ac.uk/vol1/fastq/SRR020/SRR020491/SRR020491.fastq.gz

# OSD #Illumina #1.1 Gb #largest in dataset
#curl -O ftp.sra.ebi.ac.uk/vol1/fastq/ERR771/ERR771040/ERR771040_1.fastq.gz
