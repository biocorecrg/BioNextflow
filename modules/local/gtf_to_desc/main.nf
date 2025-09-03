process GTF_TO_DESC {
    tag "${meta.id}"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'biocorecrg/almalinux-perlbrew-pyenv3':
        'biocorecrg/almalinux-perlbrew-pyenv3' }"

    input:
    tuple val(meta), path(gtffile)

    output:
    tuple val(meta), path("gene_desc.txt")         , emit: gene_desc

    when:
    task.ext.when == null || task.ext.when

    script:
    def args          = task.ext.args ?: ''
    def prefix        = task.ext.prefix ?: "${meta.id}"

    """
        conv_ens_gtf_to_desc.sh ${gtffile}
    """

}
