////////////////////////////////////////////////////
/* --         LOCAL PARAMETER VALUES           -- */
////////////////////////////////////////////////////

params.summary_params = [:]

////////////////////////////////////////////////////
/* --          VALIDATE INPUTS                 -- */
////////////////////////////////////////////////////

// Check input path parameters to see if they exist
checkPathParamList = [ params.input,  params.multiqc_config, params.fasta, params.annotationfile ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }
if (params.fasta) { ch_fasta = file(params.fasta) } else { exit 1, 'Adapter file not specified!' }
if (params.annotationfile) { ch_annotationfile = file(params.annotationfile) } else { exit 1, 'Annotation file not specified!' }


////////////////////////////////////////////////////
/* --          CONFIG FILES                    -- */
////////////////////////////////////////////////////

ch_multiqc_config        = file("$projectDir/assets/multiqc_config.yaml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()

////////////////////////////////////////////////////
/* --       IMPORT MODULES / SUBWORKFLOWS      -- */
////////////////////////////////////////////////////

// Don't overwrite global params.modules, create a copy instead and use that within the main script.
def modules = params.modules.clone()

def multiqc_options   = modules['multiqc']
//multiqc_options.args += params.multiqc_title ? Utils.joinModuleArgs(["--title \"$params.multiqc_title\""]) : ''

// Modules: local
include { GET_SOFTWARE_VERSIONS     } from '../modules/local/get_software_versions'       addParams( options: [publish_files : ['csv':'']] )

// Modules: nf-core/modules
include { FASTQC as FASTQC_RAW      } from '../modules/nf-core/software/fastqc/main'             addParams( options: modules['fastqc_raw']        )
include { CUTADAPT                  } from '../modules/nf-core/software/cutadapt/main'           addParams( options: modules['cutadapt']          )
include { FASTQC as FASTQC_TRIMMED  } from '../modules/nf-core/software/fastqc/main'             addParams( options: modules['fastqc_trimmed']    )
include { BWAMEM2_INDEX             } from '../modules/nf-core/software/bwamem2/index/main'      addParams( options: modules['bwamem2_index']     )
include { BWAMEM2_MEM               } from '../modules/nf-core/software/bwamem2/mem/main'        addParams( options: modules['bwamem2_mem']       )
include { MULTIQC                   } from '../modules/nf-core/software/multiqc/main'            addParams( options: multiqc_options              )
include { CNVKIT                    } from '../modules/nf-core/software/cnvkit/main'             addParams( options: modules['cnvkit' ]           )

// Subworkflows: local
include { INPUT_CHECK               } from '../subworkflows/local/input_check'                   addParams( options: [:] )
//include { CHECK_INPUT             } from '../subworkflows/check_input'                         addParams( options: [:] )

// Subworkflows: nf-core/subworkflows
def picard_markduplicates_samtools   = modules['picard_markduplicates_samtools']
//picard_markduplicates_samtools.args += params.bam_csi_index ? Utils.joinModuleArgs(['-c']) : ''

include { BAM_SORT_SAMTOOLS      } from '../subworkflows/nf-core/bam_sort_samtools'      addParams( sort_options: modules['samtools_sort'], index_options: modules['samtools_index'], stats_options: modules['samtools_stats']                                )
include { MARK_DUPLICATES_PICARD } from '../subworkflows/nf-core/mark_duplicates_picard' addParams( markduplicates_options: modules['picard_markduplicates'], samtools_index_options: picard_markduplicates_samtools, samtools_stats_options:  picard_markduplicates_samtools )

////////////////////////////////////////////////////
/* --           RUN MAIN WORKFLOW              -- */
////////////////////////////////////////////////////

// Info required for completion email and summary
def multiqc_report      = []


workflow CONVA {

    ch_software_versions = Channel.empty()

    /*
     * SUBWORKFLOW: Read in samplesheet, validate and stage input files
     */
    INPUT_CHECK ( 
        ch_input
    )
    .map {
        meta, fastq ->
            meta.id = meta.id.split('_')[0..-2].join('_')
            [ meta, fastq ] }
    .set { ch_fastq }

    /*
     * MODULE: Read QC
     */
    FASTQC_RAW (
        ch_fastq
    )
    ch_software_versions = ch_software_versions.mix(FASTQC_RAW.out.version.first().ifEmpty(null))

    
    /*
     * MODULE: Quality trimming
     */
    CUTADAPT (
        ch_fastq
    )
    ch_trimmed_reads = CUTADAPT.out.reads
    ch_trimmed_multiqc = CUTADAPT.out.log
    ch_software_versions = ch_software_versions.mix(CUTADAPT.out.version.first().ifEmpty(null))


    /*
     * MODULE: Trimmed-read QC
     */
    FASTQC_TRIMMED (
        ch_trimmed_reads
    )
    ch_software_versions = ch_software_versions.mix(FASTQC_TRIMMED.out.version.first().ifEmpty(null))


    /*
     * MODULE: Index FASTA file with BWAMEM2
     */
    BWAMEM2_INDEX (
        ch_fasta
    )
    ch_index = BWAMEM2_INDEX.out.index
	ch_software_versions = ch_software_versions.mix(BWAMEM2_INDEX.out.version.first().ifEmpty(null))


    /*
     * MODULE: Alignment with BWAMEM2
     */
    BWAMEM2_MEM (
        ch_trimmed_reads, ch_index
    )
	ch_bam = BWAMEM2_MEM.out.bam
    ch_software_versions = ch_software_versions.mix(BWAMEM2_MEM.out.version.first().ifEmpty(null))


    /*
     * SUBWORKFLOW: Sort, index and run stats on bam
     */
    BAM_SORT_SAMTOOLS (
        ch_bam
    )
	bam = BAM_SORT_SAMTOOLS.out.bam
    //ch_software_versions = ch_software_versions.mix(BAM_SORT_SAMTOOLS.out.samtools_version.first().ifEmpty(null))


    /*
     * SUBWORKFLOW: Mark duplicate reads and run stats on it
     */
    MARK_DUPLICATES_PICARD (
        bam
    )
	ch_markduplicates_multiqc = MARK_DUPLICATES_PICARD.out.metrics
    ch_software_versions = ch_software_versions.mix(MARK_DUPLICATES_PICARD.out.picard_version.first().ifEmpty(null))


    /*
     * MODULE: Copy number variation identification and quantification
     */
	lst = MARK_DUPLICATES_PICARD.out.bam.flatten().toList()
    rmlst = lst.remove(2)
    bamlst = lst.minus(rmlst) 
	
    CNVKIT (
        bamlst, ch_fasta, ch_annotationfile
	)
    ch_software_versions = ch_software_versions.mix(CNVKIT.out.version.first().ifEmpty(null))


    /*
     * MODULE: Pipeline reporting
     */
    // Get unique list of files containing version information
    ch_software_versions
        .map { it -> if (it) [ it.baseName, it ] }
        .groupTuple()
        .map { it[1][0] }
        .flatten()
        .collect()
        .set { ch_software_versions }
		
    GET_SOFTWARE_VERSIONS ( 
        ch_software_versions
    )

    /*
     * MODULE: MultiQC
     */
    workflow_summary    = Workflow.paramsSummaryMultiqc(workflow, params.summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(Channel.from(ch_multiqc_config))
    ch_multiqc_files = ch_multiqc_files.mix(ch_multiqc_custom_config.collect().ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(GET_SOFTWARE_VERSIONS.out.yaml.collect())
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC_RAW.out.zip.collect{it[1]}.ifEmpty([]))
    
    MULTIQC (
        ch_multiqc_files.collect()
    )
    multiqc_report       = MULTIQC.out.report.toList()
    ch_software_versions = ch_software_versions.mix(MULTIQC.out.version.ifEmpty(null))
}

////////////////////////////////////////////////////
/* --              COMPLETION EMAIL            -- */
////////////////////////////////////////////////////

workflow.onComplete {
    Completion.email(workflow, params, params.summary_params, projectDir, log, multiqc_report)
    Completion.summary(workflow, params, log)
}

////////////////////////////////////////////////////
/* --                  THE END                 -- */
////////////////////////////////////////////////////
