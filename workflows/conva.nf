////////////////////////////////////////////////////
/* --         LOCAL PARAMETER VALUES           -- */
////////////////////////////////////////////////////

params.summary_params = [:]

////////////////////////////////////////////////////
/* --          VALIDATE INPUTS                 -- */
////////////////////////////////////////////////////

// Check input path parameters to see if they exist
checkPathParamList = [ params.input,  params.fasta, params.annotationfile ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }
if (params.fasta) { ch_fasta = file(params.fasta) } else { exit 1, 'Adapter file not specified!' }
if (params.annotationfile) { ch_annotationfile = file(params.annotationfile) } else { exit 1, 'Annotation file not specified!' }
//if (params.vcf) { ch_vcf = file(params.vcf) } else { exit 1, 'VCF file not specified!' }
//if (params.tbi) { ch_tbi = file(params.tbi) } else { exit 1, 'TBI (VCF index) file not specified!' }

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
multiqc_options.args += params.multiqc_title ? Utils.joinModuleArgs(["--title \"$params.multiqc_title\""]) : ''

// Modules: local
include { GET_SOFTWARE_VERSIONS     } from '../modules/local/get_software_versions'         addParams( options: [publish_files : ['csv':'']] )
include { SEQUENZAUTILS_SEQZBINNING } from '../modules/local/sequenzautils_seqzbinning' addParams( options: modules['sequenzautils_seqzbinning' ]  )
include { CNVFACETS                 } from '../modules/local/cnvfacets'                 addParams( options: modules['cnvfacets'] )

// Modules: nf-core/modules
include { FASTQC as FASTQC_RAW      } from '../modules/nf-core/software/fastqc/main'        addParams( options: modules['fastqc_raw']        )
include { CUTADAPT                  } from '../modules/nf-core/software/cutadapt/main'      addParams( options: modules['cutadapt']          )
include { FASTQC as FASTQC_TRIMMED  } from '../modules/nf-core/software/fastqc/main'        addParams( options: modules['fastqc_trimmed']    )
include { BWAMEM2_INDEX             } from '../modules/nf-core/software/bwamem2/index/main' addParams( options: modules['bwamem2_index']     )
include { BWAMEM2_MEM               } from '../modules/nf-core/software/bwamem2/mem/main'   addParams( options: modules['bwamem2_mem']       )
include { MULTIQC                   } from '../modules/nf-core/software/multiqc/main'       addParams( options: multiqc_options              )
include { CNVKIT                    } from '../modules/nf-core/software/cnvkit/main'        addParams( options: modules['cnvkit' ]           )
include { SEQUENZAUTILS_GCWIGGLE    } from '../modules/nf-core/software/sequenzautils/gcwiggle/main' addParams( options: modules['sequenzautils_gcwiggle' ]  )
include { SEQUENZAUTILS_BAM2SEQZ    } from '../modules/nf-core/software/sequenzautils/bam2seqz/main' addParams( options: modules['sequenzautils_bam2seqz' ]  )

// Subworkflows: local
include { INPUT_CHECK               } from '../subworkflows/local/input_check'              addParams( options: [:]                          )

// Subworkflows: nf-core/subworkflows
def picard_markduplicates_samtools   = modules['picard_markduplicates_samtools']

include { BAM_SORT_SAMTOOLS      } from '../subworkflows/nf-core/bam_sort_samtools'      addParams( sort_options: modules['samtools_sort'], index_options: modules['samtools_index'], stats_options: modules['samtools_stats']                                )
include { MARK_DUPLICATES_PICARD } from '../subworkflows/nf-core/mark_duplicates_picard' addParams( markduplicates_options: modules['picard_markduplicates'], samtools_index_options: picard_markduplicates_samtools, samtools_stats_options:  picard_markduplicates_samtools )

////////////////////////////////////////////////////
/* --           RUN MAIN WORKFLOW              -- */
////////////////////////////////////////////////////

// Info required for completion email and summary
def multiqc_report = []

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
    ch_software_versions  = ch_software_versions.mix(FASTQC_RAW.out.version.first().ifEmpty(null))

    
    /*
     * MODULE: Quality trimming
     */
    CUTADAPT (
        ch_fastq
    )
    ch_trimmed_reads = CUTADAPT.out.reads
    //ch_trimmed_reads.view()
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
    ch_software_versions = ch_software_versions.mix(BWAMEM2_INDEX.out.version.last().ifEmpty(null))


    /*
     * MODULE: Alignment with BWAMEM2
     */
    BWAMEM2_MEM (
        ch_trimmed_reads, ch_index
    )
	ch_bam = BWAMEM2_MEM.out.bam


    /*
     * SUBWORKFLOW: Sort, index and run stats on bam
     */
    
    BAM_SORT_SAMTOOLS (
        ch_bam
    )
    bam = BAM_SORT_SAMTOOLS.out.bam
    ch_software_versions = ch_software_versions.mix(BAM_SORT_SAMTOOLS.out.version.first().ifEmpty(null))


    /*
     * SUBWORKFLOW: Mark duplicate reads and run stats on it
     */
    MARK_DUPLICATES_PICARD (
        bam
    )
    ch_samtools_stats         = MARK_DUPLICATES_PICARD.out.stats
    ch_samtools_flagstat      = MARK_DUPLICATES_PICARD.out.flagstat
    ch_samtools_idxstats      = MARK_DUPLICATES_PICARD.out.idxstats
    ch_software_versions      = ch_software_versions.mix(MARK_DUPLICATES_PICARD.out.picard_version.first().ifEmpty(null))
    //MARK_DUPLICATES_PICARD.out.bam.view()
   
    /*
     * MODULE: Copy number variation identification and quantification
     */
    // Remove  meta.ids for cnvkit and sequenza modules
    //lst1 = MARK_DUPLICATES_PICARD.out.bam.flatten().filter( ~/^.*normal.*bam/ ).toList()
    //lst1.view()
    //lst2 = MARK_DUPLICATES_PICARD.out.bam.flatten().filter( ~/^.*tumour.*bam/ ).toList()
    //lst2.view()
    //lst_bams = lst1.add(lst2)
    //println(lst_bams)
    //lst_bams.view()
    //lst.sort { it.size() }.view()
    //rmlst = lst.remove(0)
    //bamlst = lst.minus(rmlst)
    //normal_index = lst.indexOf('normal.markdup.sorted.bam')
    
    //bamlst.view()
    
    //lst2 = bamlst.remove(1)
    //ch_bams = bamlst.minus(lst2).sort()
    //ch_bams.view()

    ch_bams = MARK_DUPLICATES_PICARD.out.bam.flatten().filter( ~/^.*bam/ ).toList().sort { it.size() }.groupTuple()
    //lsts.view()
    //rmlsts = lsts.remove(0)
    //bamlsts = lsts.minus(rmlsts)
    //bamlsts.view()
    //lsts = bamlsts.remove(1)
    //ch_bamss = bamlsts.minus(lst2).toSortedList()
    //ch_bams.view()

    // Run CNVkit tool	
    CNVKIT (
        ch_bams, ch_fasta, ch_annotationfile
    )
    ch_software_versions = ch_software_versions.mix(CNVKIT.out.version.first().ifEmpty(null))

    // Run Sequenza tool
    //SEQUENZAUTILS_GCWIGGLE (
    //    ch_fasta
    //)
    //ch_wig = SEQUENZAUTILS_GCWIGGLE.out.wig

    // Remove meta.id
    //lst2 = bamlst.remove(0)
    //ch_bam4seqz = bamlst.minus(bamlst.remove(0))

    //SEQUENZAUTILS_BAM2SEQZ (
    //    ch_bams, ch_fasta, ch_wig
    //)
    //ch_seqz = SEQUENZAUTILS_BAM2SEQZ.out.seqz    

    
    //SEQUENZAUTILS_SEQZBINNING (
    //    ch_seqz
    //)


    /*
     * MODULE: CNV_FACETS
     */
    ch_merged_bam_bai = MARK_DUPLICATES_PICARD.out.bam.join(MARK_DUPLICATES_PICARD.out.bai)
    //ch_merged_bam_bai.view() 
    ch_normalbam_bai = ch_merged_bam_bai.filter( ~/^.*normal.markdup.sorted.bam.bai]/ )
    //ch_normalbam_bai.view()
    ch_tumourbam_bai = ch_merged_bam_bai.filter( ~/^.*tumour.markdup.sorted.bam.bai]/ )      
    ch_tumourbam_bai.view()


    //CNVFACETS (
    //    ch_normalbam_bai, ch_tumourbam_bai, ch_vcf, ch_tbi    
    //)

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
    if (!params.skip_multiqc) {
        workflow_summary    = Workflow.paramsSummaryMultiqc(workflow, params.summary_params)
        ch_workflow_summary = Channel.value(workflow_summary)

        ch_multiqc_files = Channel.empty()
        ch_multiqc_files = ch_multiqc_files.mix(Channel.from(ch_multiqc_config))
        ch_multiqc_files = ch_multiqc_files.mix(ch_multiqc_custom_config.collect().ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
        ch_multiqc_files = ch_multiqc_files.mix(GET_SOFTWARE_VERSIONS.out.yaml.collect())
        ch_multiqc_files = ch_multiqc_files.mix(FASTQC_RAW.out.zip.collect{it[1]}.ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(FASTQC_TRIMMED.out.zip.collect{it[1]}.ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(MARK_DUPLICATES_PICARD.out.metrics.collect{it[1]}.ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(ch_samtools_stats.collect{it[1]}.ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(ch_samtools_flagstat.collect{it[1]}.ifEmpty([]))

        MULTIQC (
                ch_multiqc_files.collect()
        )
            
        multiqc_report       = MULTIQC.out.report.toList()
        ch_software_versions = ch_software_versions.mix(MULTIQC.out.version.first().ifEmpty(null))

        }
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
