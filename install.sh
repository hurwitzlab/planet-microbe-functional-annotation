# Install the conda environments with bowtie2, kraken2, bracken, FragGeneScan, and snakemake
conda env create -f ./conda_envs/pm_env.yml
conda env create -f ./conda_envs/kraken2.yml
conda env create -f ./conda_envs/bracken.yml

mkdir tools
cd tools

#Get InterProScan Version
wget http://www.ebi.ac.uk/interpro/match-lookup/version
VERSION=$(head ./version | cut -d ':' -f2)
rm ./version

# Install InterProScan
wget http://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/${VERSION}/interproscan-${VERSION}-64-bit.tar.gz
wget https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/${VERSION}/interproscan-${VERSION}-64-bit.tar.gz.md5

# Recommended checksum to confirm the download was successful:
md5sum -c interproscan-${VERSION}-64-bit.tar.gz.md5
# Must return *interproscan-5.56-89.0-64-bit.tar.gz: OK*
# If not - try downloading the file again as it may be a corrupted copy.

tar -pxvzf interproscan-${VERSION}-64-bit.tar.gz
rm interproscan-${VERSION}-64-bit.tar.gz
cd interproscan-${VERSION}
python3 initial_setup.py
cd ..

# Install InterProScan's local match server
wget https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/${VERSION}/lookup_service_${VERSION}.tar.gz
wget https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/${VERSION}/lookup_service_${VERSION}.tar.gz.md5

# Recommended checksum to confirm the download was successful:
md5sum -c lookup_service_${VERSION}.tar.gz.md5

tar -pxvzf lookup_service_${VERSION}.tar.gz
rm lookup_service_${VERSION}.tar.gz lookup_service_${VERSION}.tar.gz.md5

# Add vsearch and FragGeneScan to PATH (add to $HOME/.bashrc, change to other location if they're installed elsewhere)
#export PATH="/path/to/planet-microbe-functional-annotation/tools":$PATH
#export PATH="/path/to/planet-microbe-functional-annotation/tools/FragGeneScan1.31":$PATH

# Download the Kraken2 database
wget https://genome-idx.s3.amazonaws.com/kraken/k2_pluspf_20230314.tar.gz
tar -xvf k2_pluspf_20230314.tar.gz

