process SHORTSTACK_FROM_BAM {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'community.wave.seqera.io/library/shortstack:4.1.1--e0fb243de40539fc' :
        'community.wave.seqera.io/library/shortstack:4.1.1--e0fb243de40539fc' }"

    input:
    tuple val(meta), path(bams)
    path(fasta)
    path(index_fold)
    path(faidx)
    path(micros)

    output:
    tuple val(meta), path("shortstack_output")                      , emit: output_folder
    tuple val(meta), path("shortstack_output/merged.bam")           , emit: bams
    tuple val(meta), path("shortstack_output/merged.bam.csi")       , emit: indexes
    tuple val(meta), path("shortstack_output/Counts.txt")           , emit: counts
    tuple val(meta), path("shortstack_output/Results.txt")          , emit: results_txt
    path  "versions.yml"                                            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args          = task.ext.args ?: ''
    def unzip         = ""

    if ("${micros}".endsWith(".gz")) {
        unzip = "zcat ${micros} > mymicrornas.fa"
	    args = args + " --known_miRNAs mymicrornas.fa"
    } else if (micros) {
        args = args + " --known_miRNAs ${micros}"
    }

    def prefix        = task.ext.prefix ?: "${meta.id}"

    """
    ln -s ${index_fold}/* .
    ${unzip}
    ShortStack \\
        ${args} \\
        --genomefile ${fasta}\\
        --bamfile ${bams}\\
        --outdir shortstack_output \\
        --threads ${task.cpus} \\

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        shortstack: \$( ShortStack --version | sed 's/^ShortStack//g')
    END_VERSIONS
    """

}
