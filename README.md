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
    --rm \
    quay.io/biocontainers/bbmap:38.87--h1296035_0 \
    bbmap.sh -Xmx24G ref=hg19_main_mask_ribo_animal_allplant_allfungus.fa.gz
```

Update the nextflow config to indicate that a pre-indexed contaminant genome is available, by making the following updates

```{bash}
foreign_genome = ""
foreign_genome_ref = "/mnt/efs/databases/contaminants"
```

### Seedfiles

>The `seedfile` should be a __TWO__ column csv file with the following headers.

```bash
sampleName,Reads
P3_H07_C-14-W3_S224,s3://czb-seqbot/fastqs/20241204_LH00572_0143_A22HMG3LT4/AshleyC/MI_MITI001-Preclinical-MouseStudy/P3_H07_C-14-W3_S224_R1_001.fastq.gz
P3_H07_C-14-W3_S224,s3://czb-seqbot/fastqs/20241204_LH00572_0143_A22HMG3LT4/AshleyC/MI_MITI001-Preclinical-MouseStudy/P3_H07_C-14-W3_S224_R2_001.fastq.gz

```

It must also be stored on S3. Convention is to store the `seedfile` at the following location.

```bash
s3://genomics-workflow-core/Results/ReadQC/<PROJECT_ID>/seedfile/seedfile.csv
```

Using a `seedfile` is the preferred mode of using the pipeline as it allows us to save the input for the pipeline alongside the outputs.

## Usage

- Preferred usage (`--seedfile`) with the ustomized quality threshold and minimum length after trimming

```{bash}
aws batch submit-job \
    --job-name nf-readqc_mlpe_test \
    --job-queue priority-maf-pipelines \
    --job-definition nextflow-production \
    --container-overrides command="FischbachLab/nf-readqc,\
"--project","00_TEST",\
"--prefix","20220215_multilane_PE",\
"--singleEnd","false",\
"--phred", "30",\
"--minlength", "50",\
"--seedfile","s3://genomics-workflow-core/Results/ReadQC/seedfiles/20241216_MI_MITI001-Preclinical-MouseStudy.seedfile.csv" "
```

__NOTE: the seedfile MUST be present on S3 before executing the above command.__

- Using `--reads` flag for a multilane single ended sample

```{bash}
aws batch submit-job \
    --job-name nf-readqc_mlse_test \
    --job-queue priority-maf-pipelines \
    --job-definition nextflow-production \
    --container-overrides command="FischbachLab/nf-readqc,\
"--project","00_TEST",\
"--prefix","20220215_multilane_SE",\
"--singleEnd","true",\
"--reads","s3://nextflow-pipelines/nf-readqc/data/test_data/random_ncbi_reads_with_duplicated_and_contaminants*_R1_*.fastq.gz" "
```

- Using `--reads` flag for single lane paired end sample

```{bash}
aws batch submit-job \
    --job-name nf-readqc_slpe_test \
    --job-queue priority-maf-pipelines \
    --job-definition nextflow-production \
    --container-overrides command="FischbachLab/nf-readqc,\
"--project","00_TEST",\
"--prefix","20220215_singlelane_PE",\
"--singleEnd","false",\
"--reads","'s3://czb-seqbot/fastqs/200817_NB501938_0185_AH23FNBGXG/MITI_Purification_Healthy/E8_SH0000236_0619-Cult-2-481*_R{1,2}_*.fastq.gz'" "
```

## Updating the pipeline

Simply updating the github repo will NOT update the pipeline. The repo also needs to be synced with it's S3 location:

```{bash}
cd nf-readqc
aws s3 sync . s3://nextflow-pipelines/nf-readqc --exclude ".git/*" --delete
```
