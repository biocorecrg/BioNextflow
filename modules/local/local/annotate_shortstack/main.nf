process ANNOTATE_SHORTSTACK {
    tag "${meta.id}"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'community.wave.seqera.io/library/bedtools_pip_grep:fe694158e08503cb':
        'community.wave.seqera.io/library/bedtools_pip_grep:fe694158e08503cb' }"

    input:
    tuple val(meta), path(short_res)
    tuple val(meta2), path(short_counts)
    path(annotation)
    val(strand)

    output:
    tuple val(meta), path("annotated_counts.txt.gz")            , emit: annotated_counts
    path  "versions.yml"                                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args          = task.ext.args ?: ''
    def bedopt        = ""
    if (strand == "stranded") {
    	bedopt        = "-s"
    }

    if ("${annotation}".endsWith(".gz")) {
        exec = "zcat ${annotation}"
    } else {
        exec = "cat ${annotation}"
    }

    def prefix        = task.ext.prefix ?: "${meta.id}"

    """
       ${exec} > mygenes.gtf
       awk '{print \$3"\t"\$4"\t"\$5"\t"\$1"\t"\$19"-"\$20"\t"\$10}' ${short_res} | grep -v "Start" > Results.bed
       bedtools intersect -a Results.bed -b mygenes.gtf -loj ${bedopt} | gzip -c > intersect.txt.gz
       parse_bed_anno_to_tsv.py -i intersect.txt.gz -o half_table.txt
	   cut -f4-  Counts.txt| sed s/_trimmed_condensed//g  > raw_counts.txt
       paste half_table.txt raw_counts.txt | gzip -c > annotated_counts.txt.gz
       rm mygenes.gtf Results.bed raw_counts.txt half_table.txt intersect.txt.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$( bedtools --version | sed 's/^bedtools//g')
    END_VERSIONS
    """

}
