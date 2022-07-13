# planet-microbe-functional-annotation
EBI Functional Annotation Pipeline For Planet Microbe
Matt Miller, Kai Blumberg

Based on the EBI functional pipeline version 4.1 (https://www.ebi.ac.uk/metagenomics/pipelines/4.1)

Pipeline steps: 
1) Bowtie2 (http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)

2) Trimmomatic (http://www.usadellab.org/cms/?page=trimmomatic)

3) vsearch (https://github.com/torognes/vsearch)

(These next steps branch off from one another and are done separately, 4a -> 5a and 4b -> 5b)

4a) Kraken2 (https://ccb.jhu.edu/software/kraken2/)

5a) Bracken (https://ccb.jhu.edu/software/bracken/)

4b) FragGeneScan (https://sourceforge.net/projects/fraggenescan/files/)

5b) InterProScan (https://github.com/ebi-pf-team/interproscan)

## Required Software
Python3
Singularity **Check Singularity version**
Anaconda
SLURM **Prolly state this elsewhere**
Java (we used jdk-11.0.8)

## Tool and Database Installation

The pipeline uses Snakemake to submit SLURM jobs for each step of the pipeline. The instructions below install the tools to a default directory, `planet-microbe-functional-annotation/tools`. The tools `vsearch`, `FragGeneScan1.31`, and `Trimmomatic-0.39` can be installed elsewhere, but if they are not in `planet-microbe-functional-annotation/tools`, you will need to export environmental variables for them like so: 
```
export VSEARCH="/path/to/vsearch"
export FRAGGENESCAN="/path/to/run_FragGeneScan.pl"
export TRIMMOMATIC="/path/to/trimmomatic-0.39.jar"
```

**InterProScan and the lookup server MUST BE INSTALLED IN THE TOOLS DIRECTORY!!!**

Installation directions:
```
# Install the pipeline
git clone git@github.com:hurwitzlab/planet-microbe-functional-annotation.git
cd planet-microbe-functional-annotation

# Install the conda environments with bowtie2, kraken2, bracken, FragGeneScan, and snakemake
conda env create -f ./conda_envs/pm_env.yml
conda env create -f ./conda_envs/kraken2.yml
conda env create -f ./conda_envs/bracken.yml

# Install Trimmmomatic
cd tools
wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip
unzip Trimmomatic-0.39.zip
rm Trimmomatic-0.39.zip

# Install vsearch
wget https://github.com/torognes/vsearch/releases/download/v2.21.1/vsearch-2.21.1-linux-x86_64.tar.gz
tar -xvf vsearch-2.21.1-linux-x86_64.tar.gz
rm vsearch-2.21.1-linux-x86_64.tar.gz
cp vsearch-2.21.1-linux-x86_64/bin/vsearch .
rm -r vsearch-2.21.1-linux-x86_64

# Install InterProScan
wget http://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.56-89.0/interproscan-5.56-89.0-64-bit.tar.gz
wget https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.56-89.0/interproscan-5.56-89.0-64-bit.tar.gz.md5

# Recommended checksum to confirm the download was successful:
md5sum -c interproscan-5.56-89.0-64-bit.tar.gz.md5
# Must return *interproscan-5.56-89.0-64-bit.tar.gz: OK*
# If not - try downloading the file again as it may be a corrupted copy.

tar -pxvzf interproscan-5.56-89.0-64-bit.tar.gz
rm interproscan-5.56-89.0-64-bit.tar.gz
cd interproscan-5.56-89.0
python3 initial_setup.py
cd ..

# Install InterProScan's local match server
wget https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/lookup_service/lookup_service_5.56-89.0.tar.gz
wget https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/lookup_service/lookup_service_5.56-89.0.tar.gz.md5

# Recommended checksum to confirm the download was successful:
md5sum -c lookup_service_5.56-89.0.tar.gz.md5

tar -pxvzf lookup_service_5.56-89.0.tar.gz
rm lookup_service_5.56-89.0.tar.gz lookup_service_5.56-89.0.tar.gz.md5
```

**TODO: files that need paths fixed:
bash/edit_interproscan.sh
bash/run_start_lookup_server.sh**

**TODO: export the conda envs to ymls when finished**

## Running the pipeline
The pipeline runs SLURM jobs for each step. The file `config/cluster.yml` is used by Snakemake to submit SLURM jobs. Change `partition`, `group`, and `M` to your HPC's information. The config file `config/config.yml` contains settings and parameters for the various tools used in the pipeline, and the list of samples the pipeline will be ran on. Add your samples to the list.
