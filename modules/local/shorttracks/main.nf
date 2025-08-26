process SHORTTRACKS {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'quay.io/biocontainers/shorttracks:1.2--hdfd78af_0' :
        'quay.io/biocontainers/shorttracks:1.2--hdfd78af_0' }"

    input:
    tuple val(meta), path(bam), path(csi)

    output:
    tuple val(meta), path("*.bw")                 , emit: bigwig
    path  "versions.yml"                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args          = task.ext.args ?: ''
    def prefix        = task.ext.prefix ?: "${meta.id}"

    """
    ShortTracks \\
        ${args} \\
        --bamfile ${bam}\\

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        shortstack: \$( ShortTracks --version | sed 's/^ShortTracks//g')
    END_VERSIONS
    """

}
