process SPACEMAKE_ADD_SAMPLE {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'biocorecrg/spacemake:0.9.1' :
        'biocorecrg/spacemake:0.9.1' }"

    input:
    tuple val(meta),  path(reads)
    tuple val(meta2), path(puck)
    tuple val(meta3), path(config)
    tuple val(meta4), path(barcode_folder)

    output:
    path  "versions.yml"                    , emit: versions
    tuple val(meta), path("project_df.csv") , emit: csv

    when:
    task.ext.when == null || task.ext.when

    script:
    def args          = task.ext.args ?: ''
    def prefix        = task.ext.prefix ?: "${meta.id}"
    def reads1 = []
    def reads2 = []
    meta.single_end ? [reads].flatten().each { reads1 << it } : reads.eachWithIndex { v, ix -> (ix & 1 ? reads2 : reads1) << v }
    def input_reads = meta.single_end ? "--R1 ${reads1.join(" ")}" : "--R1 ${reads1.join(" ")} --R2 ${reads2.join(" ")}"

    """
	spacemake projects add_sample \
	   --project-id project \
	   --sample-id ${meta.id} ${input_reads} \
	   --species ${meta2.id}  \
	   --puck openst \
	   --puck-barcode-file ${barcode_folder}/*.txt.gz \
	   --run_mode openst \
	   --barcode-flavor openst

	
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
