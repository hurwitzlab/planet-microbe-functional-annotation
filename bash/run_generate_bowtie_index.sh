#PBS -l select=1:ncpus=28:mem=168gb
#PBS -l walltime=12:00:00

set -u
cd $INDEX_DIR
module load singularity

ALL_FILES=$(cat $GENOME_FILELIST | tr "\n" ",")
echo "singularity exec /xdisk/bhurwitz/mig2020/rsgrps/bhurwitz/mattmiller899/planet_microbe/planet-microbe-functional-annotation/singularity/bowtie.simg /bowtie2/bowtie2-2.4.2/bowtie2-build $ALL_FILES $INDEX_NAME"
singularity exec /xdisk/bhurwitz/mig2020/rsgrps/bhurwitz/mattmiller899/planet_microbe/planet-microbe-functional-annotation/singularity/bowtie.simg /bowtie2/bowtie2-2.4.2/bowtie2-build $ALL_FILES $INDEX_NAME
