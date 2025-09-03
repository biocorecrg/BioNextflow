process GTF_TO_TX2GENE {
    tag "${meta.id}"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'biocorecrg/almalinux-perlbrew-pyenv3':
        'biocorecrg/almalinux-perlbrew-pyenv3' }"

    input:
    tuple val(meta), path(gtffile)

    output:
    tuple val(meta), path("tx2gene.csv")         , emit: tx2gene

    when:
    task.ext.when == null || task.ext.when

    script:
    def args          = task.ext.args ?: ''
    def prefix        = task.ext.prefix ?: "${meta.id}"

    """
        conv_ens_gtf_to_tx2gene.sh ${gtffile}
    """

}
