# nf-readqc

## Acknowledgement

This pipeline is based on the original [YAMP](https://github.com/alesssia/YAMP) repo. Modifications have been made to make use of our infrastrucutre more readily. If you're here for a more customizable and flexible pipeline, please consider taking a look at the original repo.

## Databases

### Contaminants

You'll need an instance with at least 24GB of RAM.

```{bash}
cd /mnt/efs/databases/contaminants
wget https://zenodo.org/record/4629921/files/hg19_main_mask_ribo_animal_allplant_allfungus.fa.gz # or latest url
docker container run \
    --volume $PWD:$PWD \
    --workdir $PWD \
    -it \
    --rm \
    quay.io/biocontainers/bbmap:38.87--h1296035_0 \
    bbmap.sh -Xmx24G ref=hg19_main_mask_ribo_animal_allplant_allfungus.fa.gz
```

Update the nextflow config to indicate that a pre-indexed contaminant genome is available, by making the following update

```{bash}
foreign_genome = ""
foreign_genome_ref = "/mnt/efs/databases/contaminants"
```

## Usage

```{bash}
aws batch submit-job \
    --profile maf \
    --job-name nf-readqc-1210-1 \
    --job-queue priority-maf-pipelines \
    --job-definition nextflow-production \
    --container-overrides command=s3://nextflow-pipelines/nf-readqc,\
"--outdir","s3://genomics-workflow-core/Results/ReadQC/00_tests",\
"--prefix","paired_end_QC",\
"--singleEnd","false",\
"--dedupe","true",\
"--reads1","s3://czb-seqbot/fastqs/200817_NB501938_0185_AH23FNBGXG/MITI_Purification_Healthy/E8_SH0000236_0619-Cult-2-481_S22_R1_001.fastq.gz",\
"--reads2","s3://czb-seqbot/fastqs/200817_NB501938_0185_AH23FNBGXG/MITI_Purification_Healthy/E8_SH0000236_0619-Cult-2-481_S22_R2_001.fastq.gz"
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
