# nf-readqc

## Acknowledgement

This pipeline is based on the original [YAMP](https://github.com/alesssia/YAMP) repo. Modifications have been made to make use of our infrastrucutre more readily. If you're here for a more customizable and flexible pipeline, please consider taking a look at the original repo.

## Usage

```{bash}
aws batch submit-job \
    --profile maf \
    --job-name nf-readqc-0825-3 \
    --job-queue default-maf-pipelines \
    --job-definition nextflow-production \
    --container-overrides command=s3://nextflow-pipelines/nf-readqc,\
"--outdir","s3://genomics-workflow-core/Results/ReadQC/00_tests",\
"--prefix","paired_end_QC",\
"--singleEnd","false",\
"--dedupe","true",\
"--reads1","s3://nextflow-pipelines/nf-readqc/data/test_data/random_ncbi_reads_with_duplicated_and_contaminants_R1.fastq.gz",\
"--reads2","s3://nextflow-pipelines/nf-readqc/data/test_data/random_ncbi_reads_with_duplicated_and_contaminants_R2.fastq.gz"
```

### Local Test

```{bash}
nextflow run . \
--outdir s3://gwfcore-results/Results/ReadQC/00_tests \
--prefix paired_end_QC \
--singleEnd false \
--reads1 s3://nextflow-pipelines/nf-readqc/data/test_data/random_ncbi_reads_with_duplicated_and_contaminants_R1.fastq.gz \
--reads2 s3://nextflow-pipelines/nf-readqc/data/test_data/random_ncbi_reads_with_duplicated_and_contaminants_R2.fastq.gz
```
