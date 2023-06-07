# Planet Microbe Functional Annotation
EBI Functional Annotation Pipeline For Planet Microbe
Kai Blumberg, Matthew Miller, Alise Ponsero, Bonnie Hurwitz

Based on the EBI functional pipeline version 4.1 (https://www.ebi.ac.uk/metagenomics/pipelines/4.1)

This pipeline analyzes metagenomic datasets in a SLURM HPC environment using Snakemake. It generates taxonomic classifications of reads using Kraken2 and Bracken and functional analyses of reads using InterProScan.

## Pipeline steps: 

1) [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) for removing human contamination

2) [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) for trimming

3) [vsearch](https://github.com/torognes/vsearch) for QC

(These next steps branch off from one another and are done separately, 4a -> 5a for taxonomic classification and 4b -> 5b for functional analysis)

4a) [Kraken2](https://ccb.jhu.edu/software/kraken2/)

5a) [Bracken](https://ccb.jhu.edu/software/bracken/)

4b) [FragGeneScan](https://sourceforge.net/projects/fraggenescan/files/)

5b) [InterProScan](https://github.com/ebi-pf-team/interproscan)


## Required Software
Python3.9

Anaconda

SLURM 

Snakemake
## Tool and Database Installation

## WARNING: The InterProScan Lookup Service download will be very large (~2 TB). Ensure you have space for it. It will take a while to download as well depending on your download speeds.
This script will create all of the necessary conda environments, install the InterProScan software and its lookup service, and download/build the indices for Kraken2/Bracken and Bowtie2.
Installation directions:
```
sh install.sh
```

## Configuring the pipeline
The pipeline runs SLURM jobs for each step. There are 3 files that users can/should make changes to.

1. The file `config/cluster.yml` is used by Snakemake to submit SLURM jobs, so you will need to add your credentials to the file. Change `partition`, `group`, and `M` to your HPC's user information. 

2. The file `config/config.yml` contains settings and parameters for the various tools used in the pipeline, and importantly, the list of samples the pipeline will be ran on. Under `samples:`, add your sample names following the format of
```
ID1: "data/ID1"
ID2: "data/ID2"
```
**NOTE**: Input files must be in the `data` directory and must be fastq format with one of these extensions:
* .fastq
* .fastq.gz
* .fq
* .fq.gz

Under `pipeline:` are named parameters for the tools used. The defaults are what we found worked best for our datasets, but they may need to be tweaked for other datasets depending on read length, which sequencer the reads came from, etc.
```
    debug: Debugging message, set to 0 to turn off
    adapter_fasta: Adapter file for Trimmomatic
    seed_mismatches: Max number of seed mismatches allowed for Trimmomatic
    palindrome_clip_thresh: palindromeClipThreshold parameter for Trimmomatic
    simple_clip_thresh: simpleClipThreshold parameter for Trimmomatic
    min_adapter_length: Minimum length of a sequence that can be considered an adapter for Trimmomatic
    min_quality: Minimum quality for leading and trailing bases, otherwise get trimmed, for Trimmomatic
    trim_min_length: Minimum length of reads to be kept after trimming, shorter reads are discarded, for Trimmomatic
    vsearch_filter_maxee: Discards sequences with more than specified number of errors, for vsearch
    vsearch_filter_minlen: Discards sequences with less than specified number of bases, for vsearch
    frag_train_file:  File name that contains model parameters for FragGeneScan
    delete_intermediates: Whether to keep intermediate results or delete them to save space (0 = keep, 1 = delete)
```
Under `snakemake:` are memory and thread parameters for the jobs. Bowtie2, Trimmomatic, vsearch, FragGeneScan, and InterProScan use `mem_small` and `threads_small` while Kraken2 and Bracken use `mem_big` and `threads_big`. Adjust these depending on the resources you need/have access to on your HPC.

3. The file `submit_snakemake.sh` is used to submit the pipeline as a job to the HPC system. Line 6 contains a variable JOBNAME that can be changed to keep track of different jobs.


## Running the pipeline

To submit your pipeline to the SLURM job scheduler, run `sh submit_snakemake.sh`.

The pipeline will output errors and outputs to files in the `err/` and `out/` directories. Some steps of the pipeline also output `log` files in their subdirectories that contain output/errors from the tool of that step.

The pipeline will output results into the `results` directory, with each sample having their own subdirectory (e.g. `results/ID1`). Each step of the pipeline will have separate directories within `results` for their output. The final functional analysis results will be found at `results/ID/step_07_combine_tsv` and the final taxonomic classifications will be found at `results/ID/kraken`.

## Example

Below are instructions for running the pipeline with a small example file.

#### Download the example data
```
cd data
sh get_test_data.sh
cd ..
```
#### Edit the config file with the input file
In `config/config.yml` under `samples:`, we will put `ERR771001_1: "data/ERR771001_1"`

#### Edit the cluster file with your HPC credentials
In `config/cluster.yml`, change the parameters to your credentials.

#### Edit the submit script with a job name
In `submit_snakemake.sh`, edit line 6 to change the job name, such as to `JOBNAME="snakemake_name"`.

#### Submitting and running the pipeline
Run the command `sh submit_snakemake.sh` to submit the job to your HPC.

#### Results


## Contributors
Matthew Miller developed the pipeline.

Kai Blumberg developed the initial idea, fixed bugs, and performed all of the data analysis using the pipeline.


