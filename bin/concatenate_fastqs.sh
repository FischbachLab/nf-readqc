#!/usr/bin/bash -x

set -euo pipefail

PREFIX=$1
SINGLE_ENDED=$2

shift 2
LIST_OF_FILES=$@

echo $LIST_OF_FILES

if [ "$SINGLE_ENDED" = true ]; then
    
    ARRAY_OF_FILES=($LIST_OF_FILES)

    if [ ${#ARRAY_OF_FILES[@]} -eq 1 ]; then
        cp ${LIST_OF_FILES} ${PREFIX}.merged.R1.fastq.gz
    else
        zcat ${LIST_OF_FILES} | gzip > ${PREFIX}.merged.R1.fastq.gz
    fi
else
    # hack to separate forward and reverse files from the list
    # ToDo @sunitj: Do better!
    FWD_LIST=$(echo ${LIST_OF_FILES} |  tr ' ' '\n' | grep "_R1_" | tr '\n' ' ' | sort)
    REV_LIST=$(echo ${LIST_OF_FILES} |  tr ' ' '\n' | grep "_R2_" | tr '\n' ' ' | sort)

    FWD_LIST_ARR=($FWD_LIST)

    if [ ${#FWD_LIST_ARR[@]} -eq 1 ]; then
        cp ${FWD_LIST} ${PREFIX}.merged.R1.fastq.gz
        cp ${REV_LIST} ${PREFIX}.merged.R2.fastq.gz
    else
        zcat $FWD_LIST | gzip > ${PREFIX}.merged.R1.fastq.gz
        zcat $REV_LIST | gzip > ${PREFIX}.merged.R2.fastq.gz
    fi
fi