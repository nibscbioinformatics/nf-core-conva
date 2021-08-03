// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process CNVFACETS {
    tag "$meta.id"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::cnv_facets=0.16.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/cnv_facets:0.16.0--py38r36h4b26f60_1"
    } else {
        container "quay.io/biocontainers/cnv_facets:0.16.0--py38r36h4b26f60_1"
    }

    input:
    tuple val(meta), path(normalbam), path(normalbai)
    tuple val(meta), path(tumourbam), path(tumourbai)
    path vcf 
    path tbi

    output:
    path "*", emit: outdir
    path "*.version.txt"          , emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    cnv_facets.R \\
        -n $normalbam \\
        -t $tumourbam \\
        -vcf $vcf \\
        $options.args \\
        -o ${prefix}

    echo \$(cnv_facets --version 2>&1) | sed 's/^.*cnv_facets //; s/Using.*\$//' > ${software}.version.txt
    """
}
