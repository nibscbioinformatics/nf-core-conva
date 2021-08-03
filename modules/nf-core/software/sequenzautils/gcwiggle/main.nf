// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process SEQUENZAUTILS_GCWIGGLE {
    tag "$fasta"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:[:], publish_by_meta:[]) }

    conda (params.enable_conda ? "bioconda::sequenza-utils=3.0.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/sequenza-utils:3.0.0--py38h6ed170a_2"
    } else {
        container "quay.io/biocontainers/sequenza-utils:3.0.0--py38h6ed170a_2"
    }

    input:
    path fasta

    output:
    path "*.wig.gz", emit: wig
    path "*.version.txt", emit: version

    script:
    def software = getSoftwareName(task.process)
    """
    sequenza-utils \\
        gc_wiggle \\
        $options.args \\
        --fasta $fasta \\
        -o ${fasta}.wig.gz

    echo \$(sequenzautils --version 2>&1) | sed 's/^.*sequenzautils //; s/Using.*\$//' > ${software}.version.txt
    """
}
