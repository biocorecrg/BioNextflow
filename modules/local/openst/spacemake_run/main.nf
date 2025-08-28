process SPACEMAKE_RUN {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'biocorecrg/spacemake:0.9.1' :
        'biocorecrg/spacemake:0.9.1' }"

    input:
    tuple val(meta),  path(reads)
    tuple val(meta2), path(comb_csv)
    tuple val(meta3), path(puck)
    tuple val(meta4), path(config)
    tuple val(meta5), path(barcode_folder)
    tuple val(meta6), path(genome)
    tuple val(meta7), path(annotation)

    output:
    path  "versions.yml"                    , emit: versions
    tuple val(meta), path("project_df.csv") , emit: csv

    when:
    task.ext.when == null || task.ext.when

    script:

    """
	spacemake run --cores ${task.cpus}
	
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
