process SPACEMAKE_ADD_SPECIES {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'biocorecrg/spacemake:0.9.1' :
        'biocorecrg/spacemake:0.9.1' }"

    input:
    tuple val(meta), path(genome)
    tuple val(meta2), path(annotation)

    output:
    tuple val(meta), path("puck_data")   , emit: puck
    tuple val(meta), path("config.yaml") , emit: config
    path  "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args          = task.ext.args ?: ''
    def prefix        = task.ext.prefix ?: "${meta.id}"

    """
    spacemake init --dropseq-tools /opt/conda/bin/ --temp-dir ./
    spacemake config add_species  --name ${meta.id}  --sequence ${genome}  --annotation ${annotation} 

	
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        spacemake: \$(  spacemake --version --version  )
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_fc_tiles

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        spacemake: \$(  spacemake --version --version  )
    END_VERSIONS
    """
}
