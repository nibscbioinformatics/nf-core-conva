// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process SEQUENZAUTILS_BAM2SEQZ {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::sequenza-utils=3.0.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/sequenza-utils:3.0.0--py38h6ed170a_2"
    } else {
        container "quay.io/biocontainers/sequenza-utils:3.0.0--py38h6ed170a_2"
    }

    input:
    tuple val(meta), path(normalbam), path(tumourbam)
    path fasta
    path wigfile

    output:
    tuple val(meta), path("*.seqz.gz")  , emit: seqz
    path "*.version.txt", emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    sequenza-utils \\
        bam2seqz \\
        $options.args \\
        -t $tumourbam \\
        -n $normalbam \\
        --fasta $fasta \\
        -gc $wigfile \\
        -o ${prefix}.seqz.gz

    echo \$(sequenzautils --version 2>&1) | sed 's/^.*sequenzautils //; s/Using.*\$//' > ${software}.version.txt
    """
}
