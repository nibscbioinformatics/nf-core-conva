// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process CNVKIT {
    tag "$tumourbam"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:[:], publish_by_meta:[]) }

    conda (params.enable_conda ? "bioconda::cnvkit=0.9.8" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/cnvkit:0.9.8--py_0"
    } else {
        container "quay.io/biocontainers/cnvkit:0.9.8--py_0"
    }

    input:
    tuple path(normalbam), path(tumourbam)
    path fasta
    path annotationfile

    output:
    path "output", emit: calls
    path "*.version.txt", emit: version

    script:
    def software = getSoftwareName(task.process)
    """
    cnvkit.py batch \\
        $tumourbam \\
        --normal $normalbam\\
        --fasta $fasta \\
        --annotate $annotationfile \\
        $options.args \\
        --output-dir output/

    cnvkit.py version | sed -e "s/cnvkit v//g" > ${software}.version.txt
    """
}
