/**
    Gets software version.

    This process ensures that software version are included in the logs.
*/
process get_software_versions {

    container params.docker_container_multiqc

    input:
    val (some_value)

    output:
    path "software_versions_mqc.yaml"

    script:
    """
    echo $workflow.manifest.version > v_pipeline.txt
    echo $workflow.nextflow.version > v_nextflow.txt

    echo $params.docker_container_fastqc | cut -d: -f 2 > v_fastqc.txt
    echo $params.docker_container_bbmap | cut -d: -f 2 > v_bbmap.txt

    echo $params.docker_container_multiqc | cut -d: -f 2 > v_multiqc.txt

    scrape_software_versions.py > software_versions_mqc.yaml
    """
}


process concatenate {

    tag "${id}"

    container params.docker_container_multiqc

    // publishDir "${workingpath}", mode: 'copy', pattern: "*.merged.R{1,2}.fastq.gz"

    input:
    tuple val(id), path(read_files)

    output:
    tuple val(id), path("${id}.merged.R{1,2}.fastq.gz")

    script:
    """
    concatenate_fastqs.py ${id} ${params.singleEnd} $read_files
    """
}

// ------------------------------------------------------------------------------
//    MULTIQC LOGGING
// ------------------------------------------------------------------------------


/**
    Generate Logs.

    Logs generate at each analysis step are collected and processed with MultiQC
*/

process log {

    publishDir "${params.outdir}/${params.project}", mode: 'copy'

    container params.docker_container_multiqc

    input:
    path multiqc_config
    path workflow_summary
    path "software_versions_mqc.yaml"
    path "fastqc/*"
    path "dedup_mqc.yaml"
    path "synthetic_contaminants_mqc.yaml"
    path "trimming_mqc.yaml"
    path "foreign_genome_indexing_mqc.yaml"
    path "decontamination_mqc.yaml"

    output:
    path "*multiqc_report*.html", emit: multiqc_report
    path "*multiqc_data*"

    script:
    """
    multiqc --config $multiqc_config . -f
    mv multiqc_report.html ${params.prefix}_multiqc_report.html
    mv multiqc_data ${params.prefix}_multiqc_data
    """
}
