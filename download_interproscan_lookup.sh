#!/bin/bash
#PBS -N DownloadInterproscan
#PBS -m bea
#PBS -M mattmiller899@email.arizona.edu
#PBS -W group_list=bhurwitz
#PBS -q standard
#PBS -l select=1:ncpus=4:mem=24gb
#PBS -l walltime=72:00:00

set -u

cd $HOME/hurwitz/mattmiller899/my_databases/interproscan

#wget ftp://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/lookup_service/lookup_service_5.46-81.0.tar.gz
#wget ftp://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/lookup_service/lookup_service_5.46-81.0.tar.gz.md5

#md5sum -c lookup_service_5.46-81.0.tar.gz.md5

tar -pxvzf lookup_service_5.46-81.0.tar.gz
