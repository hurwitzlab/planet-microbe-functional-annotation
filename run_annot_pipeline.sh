#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mem=40gb
#SBATCH --time=24:00:00



## PBS -l select=1:ncpus=4:mem=24gb
## PBS -l walltime=12:00:00

set -u

#cd /xdisk/bhurwitz/mig2020/rsgrps/bhurwitz/mattmiller899/my_databases/interproscan/lookup_service_5.46-81.0
#java -Xmx3500m -jar server-5.46-81.0-jetty-console.war --headless --port 1234 &
#sleep 120
cd /xdisk/bhurwitz/mig2020/rsgrps/bhurwitz/mattmiller899/planet_microbe/planet-microbe-functional-annotation
source activate pyenv
start_time=`date +%s`
python pipeline/pipeline.py -c data/configs/config_${TEST_TYPE}.txt -i data/in_small/ERR771001_1.fastq.gz -o data/out_small/ 
#python pipeline/pipeline.py -c data/configs/config_big.txt
#python pipeline/pipeline.py -c data/configs/config_small.txt
#python pipeline/pipeline.py -c data/configs/config_454.txt
#python pipeline/pipeline.py -c data/configs/config_paired.txt
end_time=`date +%s`
run_time=$((end_time - start_time))
echo "run_time = ${run_time}s"
