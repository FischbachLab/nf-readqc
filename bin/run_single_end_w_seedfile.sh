#!/bin/bash -x
## USAGE: bash -x run_single_end_w_seedfile.sh PROJECT_NAME PREFIX_NAME S3_PATH_2_SEEDFILE
## NOTE1: SEEDFILE must be an S3 Path.
## Note2: Convention is to store the seedfile in the expected output dir as:
##          s3://gwfcore-results/Results/ReadQC/PROJECT/00_seedfile/SEEDFILE.csv

set -euo pipefail

PROJECT="${1}"
PREFIX="${2}"
SEED="${3}"

# TODO: @sunitj
# if seedfile is local, upload to S3 and use that path

aws batch submit-job \
    --job-name nf-readqc-single-end-${PREFIX} \
    --job-queue priority-nextflow-omics \
    --job-definition nextflow-development \
    --container-overrides command=s3://nh-pipelines/nf-readqc,\
"--outdir","s3://gwfcore-results/Results/ReadQC",\
"--project","${PROJECT}",\
"--prefix","${PREFIX}",\
"--singleEnd","true",\
"--seedfile","${SEED}"