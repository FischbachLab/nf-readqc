#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { dedup; remove_synthetic_contaminants; trim; index_foreign_genome; decontaminate; quality_assessment } from './modules/quality_control'
include { concatenate; get_software_versions; log } from './modules/house_keeping'


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
if (params.reads && params.seedfile){
   exit 1, "Input reads must be defined using either '--reads' or '--seedfile' parameter. Please choose one way"
}

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



workflow {

  if(params.seedfile){
      Channel
          .fromPath(params.seedfile)
          .ifEmpty { exit 1, "Cannot find any seed file matching: ${params.seedfile}." }
          .splitCsv(header:true)
          .map{ row -> tuple(row.sampleName, row.Reads)}
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


  Channel.of(1) | get_software_versions
  workflow_summary = create_workflow_summary(summary)

  reads_concat | concatenate

  concatenate.out | dedup


  if(!params.dedup) {
      to_synthetic_contaminants = concatenate.out
      dedup_log = Channel.fromPath( "$baseDir/assets/no_dedup.yaml" )
  }else{
      to_synthetic_contaminants = dedup.out.to_synthetic_contaminants_out
      dedup_log = dedup.out.dedup_log_out
  }

  // Defines channels for resources file
  artefacts  = Channel.fromPath( params.artefacts) //"${params.artefacts}" , checkIfExists: true )
  phix174ill = Channel.fromPath( params.phix174ill)  //"${params.phix174ill}", checkIfExists: true )
  adapters = Channel.fromPath( params.adapters) //"${params.adapters}", checkIfExists: true )

  artefacts.combine(phix174ill).combine(to_synthetic_contaminants) | remove_synthetic_contaminants

  adapters.combine(remove_synthetic_contaminants.out.to_trim) | trim


  //When the indexed contaminant (pan)genome is already available, its path should be pushed in the correct channel
  foreign_genome_ch = Channel.value(1)
  //Stage boilerplate log when the contaminant (pan)genome is indexed
  if ("${params.foreign_genome_ref}" == "") {
      index_foreign_genome_log = Channel.fromPath("$projectDir/assets/foreign_genome_indexing_mqc.yaml")
  } else {
      index_foreign_genome_log = Channel.empty()
  }
  foreign_genome_ch | index_foreign_genome

  //When the indexed contaminant (pan)genome is already available, its path should be pushed in the correct channel
  if ("${params.foreign_genome_ref}" != "") {
      ref_foreign_genome = Channel.fromPath("${params.foreign_genome_ref}")
      ref_foreign_genome.combine(trim.out.to_decontaminate) | decontaminate
  }else {
      index_foreign_genome.out.combine(trim.out.to_decontaminate) | decontaminate
  }

  concatenate.out.mix(decontaminate.out.qcd_reads) | quality_assessment

}

/*

// Stage config files
multiqc_config = Channel.fromPath("${params.multiqc_config}")

log ( multiqc_config,
      workflow_summary,
      get_software_versions.out,
      quality_assessment.out.collect().ifEmpty([]),
      dedup_log.ifEmpty([]),
      remove_synthetic_contaminants.out.synthetic_contaminants_log.ifEmpty([]),
      trim.out.trimming_log.ifEmpty([]),
      index_foreign_genome_log.ifEmpty([]),
      decontaminate.out.decontaminate_log.ifEmpty([])
  )

*/
