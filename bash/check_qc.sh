#!/bin/bash
set -u
NUMS=$(cat $1 | grep "sequences kept" | cut -d" " -f1,6,8)
IFS=" " read -r -a NUMARR <<< "$NUMS"
kept=${NUMARR[0]}
trunc=${NUMARR[1]}
disc=${NUMARR[2]}
total=$((kept + disc))
mathstr="${kept}/${total}"
percent=$(bc -l <<< $mathstr)
thres=0.5
if (( $(echo "$thres > $percent" |bc -l) )); then
    echo "ERROR: percent kept by QC (${percent}) less than ${thres}, will not proceed
with rest of pipeline"
else
    echo "QC looks good, ${percent}% of reads kept"
    touch $2
fi
