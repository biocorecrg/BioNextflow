process STAR_GENOMEGENERATE {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'biocorecrg/spacemake:0.9.1' :
        'biocorecrg/spacemake:0.9.1' }"

    input:
    tuple val(meta), path(fasta)
    tuple val(meta2), path(gtf)

    output:
    tuple val(meta), path("star_index")  , emit: species_folder
    path  "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args          = task.ext.args ?: ''
    def prefix        = task.ext.prefix ?: "${meta.id}"


    """
 	/opt/conda/bin/STAR-avx2 --runMode genomeGenerate --runThreadN ${task.cpus} --genomeDir star_index --genomeFastaFiles ${fasta} --sjdbGTFfile ${gtf}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        STAR: \$(  /opt/conda/bin/STAR-avx2 --version  )
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch species_data

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        STAR: \$(  /opt/conda/bin/STAR-avx2 --version  )
    END_VERSIONS
    """
}
