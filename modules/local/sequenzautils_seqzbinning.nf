// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process SEQUENZAUTILS_SEQZBINNING {
    tag "out"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:[:], publish_by_meta:[]) }

    conda (params.enable_conda ? "bioconda::sequenza-utils=3.0.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/https://depot.galaxyproject.org/singularity/sequenza-utils:3.0.0--py39h8f06ca0_4"
    } else {
        container "quay.io/biocontainers/quay.io/biocontainers/sequenza-utils:3.0.0--py27h304d29a_4"
    }

    input:
    path seqz

    output:
    path "out.small.seqz.gz", emit: seqzbinned
    path "*.version.txt"    , emit: version

    script:
    def software = getSoftwareName(task.process)
    """
    sequenza-utils \\
        seqz_binning \\
        $options.args \\
        --seqz $seqz \\
        -o out.small.seqz.gz

    echo \$(sequenza-utils --version 2>&1) | sed 's/^.*seqenza-utils //; s/Using.*\$//' > ${software}.version.txt
    """
}
