#!/bin/bash

set -u

#OSD illumina example
curl -O ftp.sra.ebi.ac.uk/vol1/fastq/ERR771/ERR771001/ERR771001_1.fastq.gz

curl -O ftp.sra.ebi.ac.uk/vol1/fastq/ERR771/ERR771001/ERR771001_2.fastq.gz

# Tara 454 example
curl -O ftp.sra.ebi.ac.uk/vol1/fastq/ERR321/ERR321017/ERR321017.fastq.gz
