# planet-microbe-functional-annotation
EBI Functional Annotation Pipeline For Planet Microbe
Kai Blumberg, Matthew Miller, Alise Ponsero, Bonnie Hurwitz

Based on the EBI functional pipeline version 4.1 (https://www.ebi.ac.uk/metagenomics/pipelines/4.1)

This pipeline analyzes metagenomic datasets on a SLURM HPC environment using Snakemake. It generates taxonomic classifications of reads using Kraken2 and Bracken and functional analyses of the reads using InterProScan.

Pipeline steps: 
1) Bowtie2 (http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)

2) Trimmomatic (http://www.usadellab.org/cms/?page=trimmomatic)

3) vsearch (https://github.com/torognes/vsearch)

(These next steps branch off from one another and are done separately, 4a -> 5a and 4b -> 5b)

4a) Kraken2 (https://ccb.jhu.edu/software/kraken2/)

5a) Bracken (https://ccb.jhu.edu/software/bracken/)

4b) FragGeneScan (https://sourceforge.net/projects/fraggenescan/files/)

5b) InterProScan (https://github.com/ebi-pf-team/interproscan)

## INSERT IMAGE HERE

## Required Software
Python3

Anaconda

SLURM 

Snakemake
## Tool and Database Installation

# WARNING: The InterProScan Lookup Service download will be very large (~2 TB). Ensure you have space for it
Installation directions:
```
sh install.sh
```

## Building the Bowtie Index
## TODO Fix Bowtie2 script, currently hardcoded
The pipeline uses Bowtie2 to clean the datasets for human and phi-X174 contamination. To build the bowtie index, run 
`sh bash/create_bowtie_index.sh`

## Running the pipeline
The pipeline runs SLURM jobs for each step. There are 3 files that users can/should make changes to.

1. The file `config/cluster.yml` is used by Snakemake to submit SLURM jobs, so you will need to add your credentials to the file. Change `partition`, `group`, and `M` to your HPC's user information. 

2. The file `config/config.yml` contains settings and parameters for the various tools used in the pipeline, and importantly, the list of samples the pipeline will be ran on. Under `samples:`, add your sample names following the format of
```
ID1: "/path/to/ID1"
ID2: "/path/to/ID2"
```

3. The file `submit_snakemake.sh` is used to submit the pipeline as a job to the HPC system. Line 6 contains a variable JOBNAME that can be changed to keep track of different jobs.


To submit your pipeline to the SLURM job scheduler, run `sh submit_snakemake.sh`.

The pipeline will output errors and outputs to files in the `err/` and `out/` directories. Some steps of the pipeline also output `log` files in their subdirectories that contain output/errors from the tool of that step.

The pipeline will output results into the `results` directory, with each sample having their own subdirectory (e.g. `results/ID1`). Each step of the pipeline will have separate directories within `results` for their output. The final functional analysis results will be found at `results/ID/step_07_combine_tsv` and the final taxonomic classifications will be found at `results/ID/kraken`.

## Contributors
Kai Blumberg developed the initial idea, fixed bugs, and performed all of the data analysis using the pipeline.

Matthew Miller developed the pipeline.

# TODO Fix hardcoded paths for Kraken/Bracken
