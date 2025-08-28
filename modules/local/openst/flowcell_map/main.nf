process FLOWCELL_MAP {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'biocorecrg/openst:0.2.3'  :
        'biocorecrg/openst:0.2.3' }"

    input:
    tuple val(meta), path(bcl)

    output:
    tuple val(meta), path("*_fc_tiles"), emit: tiles
    path  "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args          = task.ext.args ?: ''
    def prefix        = task.ext.prefix ?: "${meta.id}"


    """
    openst flowcell_map \
       --bcl-in ${bcl} \
       --tiles-out ./${prefix}_fc_tiles \
       --parallel-processes ${task.cpus} ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        openst: \$(  openst --version  )
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_fc_tiles

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        openst: \$(  openst --version  )
    END_VERSIONS
    """
}
