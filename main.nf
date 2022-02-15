#!/usr/bin/env nextflow

def versionMessage() 
{
    log.info"""
     
    nf-readqc - Version: ${workflow.manifest.version} 
    """.stripIndent()
}

def helpMessage() 
{
    log.info"""

nf-readqc - Version: ${workflow.manifest.version} 
  
  Mandatory arguments:
    --reads   path     A glob pattern like '/some/path/SampleName_*R{1,2}.fastq.gz' or '/some/path/SampleName_*_R1_*.fastq.gz'
    --seedfile path     csv file with headers in the format "sampleName,Reads", 
                        fwd and rev is interpreted based on the expression "_R{1,2}_" in the file name.
                        lines with the same 'sampleName' are grouped together.
    --prefix   prefix  Prefix used to name the result files (only used to name outputs when --reads flag is used)
    --outdir   path    Output directory (will be outdir/project/)
    --project   path    Project directory (will be used to group outputs of the same project)
  
  Main options:
    --singleEnd  <true|false>   whether the layout is single-end
    --dedup      <true|false>   whether to perform de-duplication
  
  Other options:
  BBduk parameters for removing synthetic contaminants and trimming:
    --qin                 <33|64> Input quality offset 
    --kcontaminants       value   kmer length used for identifying contaminants
    --phred               value   regions with average quality BELOW this will be trimmed 
    --minlength           value   reads shorter than this after trimming will be discarded
    --mink                value   shorter kmer at read tips to look for 
    --hdist               value   maximum Hamming distance for ref kmer
    --artefacts           path    FASTA file with artefacts
    --phix174ill          path    FASTA file with phix174_ill
    --adapters            path    FASTA file with adapters         
  
  BBwrap parameters for decontamination:
    --foreign_genome      path    FASTA file for contaminant (pan)genome
    --foreign_genome_ref  path    folder for for contaminant (pan)genome (pre indexed)
    --mind                value   approximate minimum alignment identity to look for
    --maxindel            value   longest indel to look for
    --bwr                 value   restrict alignment band to this

nf-readQC supports FASTQ and compressed FASTQ files.
"""
}

/**
Prints version when asked for
*/
if (params.version) {
    versionMessage()
    exit 0
}

/**
Prints help when asked for
*/

if (params.help) {
    helpMessage()
    exit 0
}

/**
STEP 0. 
    
Checks input parameters and (if it does not exists) creates the directory 
where the results will be stored (aka working directory). 
Initialises the log file.
    
The working directory is named after the prefix and located in the outdir 
folder. The log file, that will save summary statistics, execution time,
and warnings generated during the pipeline execution, will be saved in the 
working directory as "prefix.log".
*/


if (params.qin != 33 && params.qin != 64) {  
    exit 1, "Input quality offset (qin) not available. Choose either 33 (ASCII+33) or 64 (ASCII+64)" 
}   

// //--reads2 can be omitted when the library layout is "single"
// if (!params.singleEnd && (params.reads2 == "null") ) {
//     exit 1, "If dealing with paired-end reads, please set the reads2 parameters\nif dealing with single-end reads, please set the '--singleEnd' option to 'true'"
// }

//Creates working dir
// fixes #1
workingpath = params.outdir + "/" + params.project
workingdir = file(workingpath)

if( !workingdir.exists() ) {
    if( !workingdir.mkdirs() )     {
        exit 1, "Cannot create working directory: $workingpath"
    } 
}    

if(params.prefix){
    workingpath = workingpath + "/" + params.prefix
}


//
if (params.reads && params.seedfile){
   exit 1, "Input reads must be defined using either '--reads' or '--seedfile' parameter. Please choose one way"
}

if(params.seedfile){
    Channel
        .fromPath(params.seedfile)
        .ifEmpty { exit 1, "Cannot find any seed file matching: ${params.seedfile}." }
        .splitCsv(header:true)
        .map{ row -> tuple(row.sampleName, file(row.Reads))}
        .groupTuple(sort:true)
        .set { reads_concat }
} else {
    Channel
        .fromPath(params.reads)
        .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}." }
        .map {it -> tuple(params.prefix, it)}
        .groupTuple(sort:true)
        .set { reads_concat }
}

// reads_concat.view()

process concatenate {
    tag "$id"

    container params.docker_container_multiqc

    // publishDir "${workingpath}", mode: 'copy', pattern: "*.merged.R{1,2}.fastq.gz"

    input:
    set val(id), file(read_files) from reads_concat

    output:
    tuple(id, file("${id}.merged.R{1,2}.fastq.gz")) into concat_fq 
    
    script:
    """
    concatenate_fastqs.py ${id} ${params.singleEnd} $read_files
    """
}

concat_fq.into {read_files_fastqc; read_files_dedup; read_files_synthetic_contaminants}

// Header log info
log.info """---------------------------------------------
nf-readQC
---------------------------------------------

Analysis introspection:

"""

def summary = [:]

summary['Starting time'] = new java.util.Date() 
//Environment
summary['Environment'] = ""
summary['Pipeline Name'] = 'nf-readQC'
summary['Pipeline Version'] = workflow.manifest.version

summary['Config Profile'] = workflow.profile
summary['Resumed'] = workflow.resume
        
summary['Nextflow version'] = nextflow.version.toString() + " build " + nextflow.build.toString() + " (" + nextflow.timestamp + ")"

summary['Java version'] = System.getProperty("java.version")
summary['Java Virtual Machine'] = System.getProperty("java.vm.name") + "(" + System.getProperty("java.vm.version") + ")"

summary['Operating system'] = System.getProperty("os.name") + " " + System.getProperty("os.arch") + " v" +  System.getProperty("os.version")
summary['User name'] = System.getProperty("user.name") //User's account name

summary['BBmap'] = "quay.io/biocontainers/bbmap:38.87--h1296035_0"
summary['FastQC'] = "quay.io/biocontainers/fastqc:0.11.9--0"

summary['MultiQC'] = "quay.io/biocontainers/multiqc:1.9--py_1"


//General
summary['Running parameters'] = ""
summary['Reads'] = params.seedfile ? params.seedfile : params.reads
summary['Prefix'] = params.prefix
summary['Layout'] = params.singleEnd ? 'Single-End' : 'Paired-End'
summary['Performing de-duplication'] = params.dedup

//remove_synthetic_contaminants 
summary['Synthetic contaminants'] = ""
summary['Artefacts'] = params.artefacts
summary['Phix174ill'] = params.phix174ill

//Trimming
summary['Adapters'] = params.adapters
summary['Trimming parameters'] = ""
summary['Input quality offset'] = params.qin == 33 ? 'ASCII+33' : 'ASCII+64'
summary['Min phred score'] = params.phred
summary['Min length'] = params.minlength
summary['kmer lenght'] = params.kcontaminants
summary['Shorter kmer'] = params.mink 
summary['Max Hamming distance'] = params.hdist 

//Decontamination
summary['Decontamination parameters'] = ""
if (params.foreign_genome_ref != "") {
    summary['Contaminant (pan)genome'] = params.foreign_genome_ref + " (indexed)"
} else if (    params.foreign_genome_ref == "") {
    summary['Contaminant (pan)genome'] = params.foreign_genome
}    
summary['Min alignment identity'] = params.mind
summary['Max indel length'] = params.maxindel
summary['Max alignment band'] = params.bwr


//Folders
summary['Folders'] = ""
summary['Output dir'] = workingpath
summary['Working dir'] = workflow.workDir
summary['Output dir'] = params.outdir
summary['Script dir'] = workflow.projectDir
summary['Lunching dir'] = workflow.launchDir

log.info summary.collect { k,v -> "${k.padRight(27)}: $v" }.join("\n")
log.info ""


/**
    Prepare workflow introspection

    This process adds the workflow introspection (also printed at runtime) in the logs
    This is NF-CORE code.
*/

def create_workflow_summary(summary) {
    def yaml_file = workDir.resolve('workflow_summary_mqc.yaml')
    yaml_file.text  = """
    id: 'workflow-summary'
    description: "This information is collected when the pipeline is started."
    section_name: 'nf-readQC Workflow Summary'
    section_href: 'https://github.com/nuancedhealth/nf-readqc'
    plot_type: 'html'
    data: |
        <dl class=\"dl-horizontal\">
${summary.collect { k,v -> "            <dt>$k</dt><dd>$v</dd>" }.join("\n")}
        </dl>
    """.stripIndent()

   return yaml_file
}

/**
    Gets software version. 

    This process ensures that software version are included in the logs.
*/
process get_software_versions {

    container params.docker_container_multiqc

    output:
    file "software_versions_mqc.yaml" into software_versions_yaml

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

/**
    Creates a set of channels for input read files.
    - read_files_fastqc is used for the first QC assessment (on the raw reads)
    - read_files_dedup  is used for the deduplication step (which is optional and may skip to trimming)
    - read_files_trim   is used for the decontamination from synthetic contaminants (used only if
      deduplication is not run)
*/

// ------------------------------------------------------------------------------   
//    QUALITY CONTROL 
// ------------------------------------------------------------------------------   

/**
    Quality Control - STEP 1. De-duplication. Only exact duplicates are removed.

    This step is OPTIONAL. De-duplication should be carried on iff you are
    using PCR amplification (in this case identical reads are technical artefacts)
    but not otherwise (identical reads will identify natural duplicates).
*/

process dedup {
    
    tag "$name"
    container params.docker_container_bbmap
        
    input:
    tuple val(name), file(reads) from read_files_dedup

    output:
    tuple val(name), path("${name}_dedup*.fq.gz") into to_synthetic_contaminants
    file "dedup_mqc.yaml" into dedup_log
    
    when:
    params.dedup

    script:
    // This is to deal with single and paired end reads
    def input = params.singleEnd ? "in=\"${reads[0]}\"" :  "in1=\"${reads[0]}\" in2=\"${reads[1]}\""
    def output = params.singleEnd ? "out=\"${name}_dedup.fq.gz\"" :  "out1=\"${name}_dedup_R1.fq.gz\" out2=\"${name}_dedup_R2.fq.gz\""
    
    """
    #Sets the maximum memory to the value requested in the config file
    maxmem=\$(echo \"$task.memory\" | sed 's/ //g' | sed 's/B//g')
    echo \"$reads\"
    clumpify.sh -Xmx\"\$maxmem\" $input $output qin=$params.qin dedupe subs=0 threads=${task.cpus} &> dedup_mqc.txt
    
    # MultiQC doesn't have a module for clumpify yet. As a consequence, I
    # had to create a YAML file with all the info I need via a bash script
    bash scrape_dedup_log.sh > dedup_mqc.yaml
    """
}

/**
    Quality control - STEP 2. A decontamination of synthetic sequences and artefacts 
    is performed.
*/

//When the de-duplication is not done, the raw file should be pushed in the correct channel
//FIXME: make this also optional?
if (!params.dedup) {
    to_synthetic_contaminants = read_files_synthetic_contaminants
    dedup_log = Channel.from(file("$baseDir/assets/no_dedup.yaml"))
}

// Defines channels for resources file 
Channel.fromPath( "${params.artefacts}", checkIfExists: true ).set { artefacts }
Channel.fromPath( "${params.phix174ill}", checkIfExists: true ).set { phix174ill }

process remove_synthetic_contaminants {
    tag "$name"
    
    container params.docker_container_bbmap

    input:
    tuple file(artefacts), file(phix174ill), val(name), file(reads) from artefacts.combine(phix174ill).combine(to_synthetic_contaminants)
   
    output:
    tuple val(name), path("${name}_no_synthetic_contaminants*.fq.gz") into to_trim
    file "synthetic_contaminants_mqc.yaml" into synthetic_contaminants_log

       script:
    def input = params.singleEnd ? "in=\"${reads[0]}\"" :  "in1=\"${reads[0]}\" in2=\"${reads[1]}\""
    def output = params.singleEnd ? "out=\"${name}_no_synthetic_contaminants.fq.gz\"" :  "out=\"${name}_no_synthetic_contaminants_R1.fq.gz\" out2=\"${name}_no_synthetic_contaminants_R2.fq.gz\""
    """
    #Sets the maximum memory to the value requested in the config file
    maxmem=\$(echo ${task.memory} | sed 's/ //g' | sed 's/B//g')
    bbduk.sh -Xmx\"\$maxmem\" $input $output k=31 ref=$phix174ill,$artefacts qin=$params.qin threads=${task.cpus} ow &> synthetic_contaminants_mqc.txt
    
    # MultiQC doesn't have a module for bbduk yet. As a consequence, I
    # had to create a YAML file with all the info I need via a bash script
    bash scrape_remove_synthetic_contaminants_log.sh > synthetic_contaminants_mqc.yaml
    """
}


/**
    Quality control - STEP 3. Trimming of low quality bases and of adapter sequences. 
    Short reads are discarded. 
    
    If dealing with paired-end reads, when either forward or reverse of a paired-read
    are discarded, the surviving read is saved on a file of singleton reads.
*/

// Defines channels for resources file 
Channel.fromPath( "${params.adapters}", checkIfExists: true ).set { adapters }

process trim {
    tag "$name"
    
    container params.docker_container_bbmap
    
    input:
    tuple file(adapters), val(name), file(reads) from adapters.combine(to_trim) 
    
    output:
    tuple val(name), path("${name}_trimmed*.fq.gz") into to_decontaminate
    file "trimming_mqc.yaml" into trimming_log

       script:
    def input = params.singleEnd ? "in=\"${reads[0]}\"" :  "in1=\"${reads[0]}\" in2=\"${reads[1]}\""
    def output = params.singleEnd ? "out=\"${name}_trimmed.fq.gz\"" :  "out=\"${name}_trimmed_R1.fq.gz\" out2=\"${name}_trimmed_R2.fq.gz\" outs=\"${name}_trimmed_singletons.fq.gz\""
    """
    #Sets the maximum memory to the value requested in the config file
    maxmem=\$(echo ${task.memory} | sed 's/ //g' | sed 's/B//g')

    bbduk.sh -Xmx\"\$maxmem\" $input $output ktrim=r k=$params.kcontaminants mink=$params.mink hdist=$params.hdist qtrim=rl trimq=$params.phred  minlength=$params.minlength ref=$adapters qin=$params.qin threads=${task.cpus} tbo tpe ow &> trimming_mqc.txt

    # MultiQC doesn't have a module for bbduk yet. As a consequence, I
    # had to create a YAML file with all the info I need via a bash script
    bash scrape_trimming_log.sh > trimming_mqc.yaml
    """
}


/**
    Quality control - STEP 4. Decontamination. Removes external organisms' contamination, 
    using given genomes. 

    When an indexed contaminant (pan)genome is not provided, the index_foreign_genome process is run 
    before the decontamination process. This process require the FASTA file of the contaminant (pan)genome.
*/

// Defines channels for foreign_genome file 
/* sunitj: pass a bogus channel to initiate the process
*       because the database is on the EFS and not on S3
*       we can pass the path to the database directly.
*       Not doing so and passing it via a channel, Nextflow tries to 
*       "stage" the file by creating a symbolic link and is unable to find
*       the file on EFS. Not sure why this happens though.
*/
foreign_genome_ch = Channel.value(1)

//Stage boilerplate log when the contaminant (pan)genome is indexed
if (params.foreign_genome_ref == "") {
    index_foreign_genome_log = Channel.from(file("$baseDir/assets/foreign_genome_indexing_mqc.yaml"))
} else {
    index_foreign_genome_log = Channel.empty()
}

process index_foreign_genome {

    container params.docker_container_bbmap

    // publishDir "${params.outdir}/${params.prefix}/foreign_genome_index", mode="copy"

    input:
    val(some_value) from foreign_genome_ch 

    output:
    path("ref/", type: 'dir') into ref_foreign_genome
    
    when:
    params.foreign_genome_ref == ""

    script:
    """
    #Sets the maximum memory to the value requested in the config file
    maxmem=\$(echo ${task.memory} | sed 's/ //g' | sed 's/B//g')

    # This step will have a boilerplate log because the information saved by bbmap are not relevant
    bbmap.sh -Xmx\"\$maxmem\" ref=${params.foreign_genome} &> foreign_genome_index_mqc.txt    
    """
}

//Channel.fromPath( "${params.foreign_genome_ref}", checkIfExists: true ).set { ref_foreign_genome }

//When the indexed contaminant (pan)genome is already available, its path should be pushed in the correct channel
if (params.foreign_genome_ref != "") {
    ref_foreign_genome = Channel.from(file(params.foreign_genome_ref))
}

process decontaminate {
    tag "$name"

    container params.docker_container_bbmap

    publishDir "${workingpath}/02_ReadQC_Output", mode: 'copy', pattern: "*.qcd.fq.gz"

    input:
    tuple path(ref_foreign_genome), val(name), file(reads) from ref_foreign_genome.combine(to_decontaminate)

    output:
    tuple val(name), path("*.qcd.fq.gz") into qcd_reads
    file "decontamination_mqc.yaml" into decontaminate_log

    script:
    // When paired-end are used, decontamination is carried on independently on paired reads
    // and on singleton reads thanks to BBwrap, that calls BBmap once on the paired reads
    // and once on the singleton ones, merging the results on a single output file
    def input = params.singleEnd ? "in=\"${reads[0]}\"" :  "in1=\"${reads[0]}\",\"${reads[2]}\" in2=\"${reads[1]}\",null"
    def output = "outu=\"${name}_QCd.fq.gz\" outm=\"${name}_contamination.fq.gz\""
    """
    #Sets the maximum memory to the value requested in the config file
    maxmem=\$(echo ${task.memory} | sed 's/ //g' | sed 's/B//g')

    bbwrap.sh \\
        -Xmx\"\$maxmem\"  \\
        mapper=bbmap \\
        append=t \\
        $input \\
        $output \\
        minid=$params.mind \\
        maxindel=$params.maxindel \\
        bwr=$params.bwr \\
        bw=12 \\
        minhits=2 \\
        qtrim=rl \\
        trimq=$params.phred \\
        path="./" \\
        qin=$params.qin \\
        threads=${task.cpus} \\
        untrim \\
        quickmatch \\
        fast \\
        ow &> decontamination_mqc.txt

    if [ $params.singleEnd = true ]; then
        # Rename
        mv ${name}_QCd.fq.gz ${name}_R1.qcd.fq.gz
    else
        # Deinterleave
        reformat.sh in=${name}_QCd.fq.gz out1=${name}_R1.qcd.fq.gz out2=${name}_R2.qcd.fq.gz
    fi

    # MultiQC doesn't have a module for bbwrap yet. As a consequence, I
    # had to create a YAML file with all the info I need via a bash script
    bash scrape_decontamination_log.sh > decontamination_mqc.yaml
    """
}

// ------------------------------------------------------------------------------   
//    QUALITY ASSESSMENT 
// ------------------------------------------------------------------------------   


process quality_assessment {
    tag "$name"
    
    container params.docker_container_fastqc
    
    publishDir "${workingpath}/01_fastqc"

    input:
    set val(name), file(reads) from read_files_fastqc.mix(qcd_reads)

    output:
    path "*_fastqc.{zip,html}" into fastqc_log

    script:
    """
    fastqc -q $reads
    """
}

// ------------------------------------------------------------------------------   
//    MULTIQC LOGGING
// ------------------------------------------------------------------------------   


/**
    Generate Logs. 

    Logs generate at each analysis step are collected and processed with MultiQC 
*/

// Stage config files
multiqc_config = file(params.multiqc_config)

process log {
    
    publishDir "${workingpath}", mode: 'copy'

    container params.docker_container_multiqc

    input:
    file multiqc_config
    file workflow_summary from create_workflow_summary(summary)
    file "software_versions_mqc.yaml" from software_versions_yaml
    path "fastqc/*" from fastqc_log.collect().ifEmpty([])
    file "dedup_mqc.yaml" from dedup_log.ifEmpty([])
    file "synthetic_contaminants_mqc.yaml" from synthetic_contaminants_log.ifEmpty([])
    file "trimming_mqc.yaml" from trimming_log.ifEmpty([])
    file "foreign_genome_indexing_mqc.yaml" from index_foreign_genome_log.ifEmpty([])
    file "decontamination_mqc.yaml" from decontaminate_log.ifEmpty([])
    
    output:
    path "*multiqc_report*.html" into multiqc_report
    path "*multiqc_data*"

    script:
    """
    multiqc --config $multiqc_config . -f
    mv multiqc_report.html ${params.prefix}_multiqc_report.html
    mv multiqc_data ${params.prefix}_multiqc_data
    """
}

