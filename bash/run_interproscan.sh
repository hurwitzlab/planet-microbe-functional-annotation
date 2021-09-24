#!/bin/bash

set -u

cd $SLURM_SUBMIT_DIR
IN_FILE="$1"
OUT_FILE="$2"
export IPS="$3"
CHUNKS_DIR="$4"
IPS_DIR="$5"
mkdir $CHUNKS_DIR
mkdir $IPS_DIR
ARGS="--partition=standard --account=bhurwitz"
ERR_DIR="./err"
OUT_DIR="./out"


IN_BN=$(basename $IN_FILE)
FILE_BN=${IN_BN%.faa}
DATASET_NAME=$(echo ${FILE_BN} | cut -d"_" -f1)
echo "datasetname = ${DATASET_NAME}"
#echo "$FILE_BN"
NUM_LINES=$(wc -l $IN_FILE | cut -d" " -f1)
CHUNK_SIZE=200000
CURR_LC=0
CURR_CHUNK=0
echo "$NUM_LINES"
while [[ $CURR_LC -lt $NUM_LINES ]]; do
    export CURR_FILE="${CHUNKS_DIR}/${FILE_BN}_${CURR_CHUNK}.faa"
    #echo "$CURR_FILE"
    BATCH=$((CHUNK_SIZE * (1 + $CURR_CHUNK)))
    #echo "$BATCH"
    if [[ $BATCH -gt $NUM_LINES ]]; then
        PREV_BATCH=$((CHUNK_SIZE * CURR_CHUNK))
        TAIL_NUM=$((NUM_LINES - PREV_BATCH))
        tail -${TAIL_NUM} ${IN_FILE} > $CURR_FILE
    else
        head -${BATCH} ${IN_FILE} | tail -${CHUNK_SIZE} > $CURR_FILE
    fi
    sed -i 's/[*]/X/g' $CURR_FILE
    export IPS_OUT="${IPS_DIR}/${FILE_BN}_${CURR_CHUNK}_interpro"
    echo "$IPS_OUT"
    JOB_ID=`sbatch $ARGS --export=IPS,IPS_OUT,CURR_FILE --job-name=ips_${CURR_CHUNK} -e $ERR_DIR/ips_chunk${CURR_CHUNK}_${DATASET_NAME}.err -o $OUT_DIR/ips_chunk${CURR_CHUNK}_${DATASET_NAME}.out ./bash/run_interproscan_chunk.sh`
    if [ "${JOB_ID}x" != "x" ]; then
        echo Job: \"$JOB_ID\"
    else
        echo Problem submitting job. Job terminated.
        exit 1
    fi
    CURR_CHUNK=$((CURR_CHUNK + 1))
    #echo "old CURR_LC=${CURR_LC}"
    CURR_LC=$((CURR_LC + CHUNK_SIZE))
    #echo "new CURR_LC=${CURR_LC}"
done
COMPLETED_CHUNKS=0
TOTAL_CHUNKS=$CURR_CHUNK
ALL_CHUNKS=$( seq 0 $((TOTAL_CHUNKS - 1)) )
while [[ $COMPLETED_CHUNKS -lt $TOTAL_CHUNKS ]]; do
    sleep 60
    for i in ${ALL_CHUNKS[@]}; do
        STDOUT_FILE="$OUT_DIR/ips_chunk${i}_${DATASET_NAME}.out"
        if grep -Fxq "Done" $STDOUT_FILE; then
            COMPLETED_CHUNKS=$((COMPLETED_CHUNKS + 1))
            ALL_CHUNKS=("${ALL_CHUNKS[@]/$i}")
            echo "Chunk ${i} successfully completed"
        elif grep -Fq "Detailed performance metrics for this job" $STDOUT_FILE; then
            echo "ERROR: Chunk ${i} failed. Exiting..."
            exit 1
        fi
    done
done
ALL_CHUNKS=$( seq 0 $((TOTAL_CHUNKS - 1)) )
touch ${OUT_FILE}
for i in ${ALL_CHUNKS[@]}; do
    CURR_IPS_OUT="${IPS_DIR}/${FILE_BN}_${CURR_CHUNK}_interpro.tsv"
    cat $CURR_IPS_OUT >> $OUT_FILE
done
echo "Successfully finished interproscan for every chunk"

