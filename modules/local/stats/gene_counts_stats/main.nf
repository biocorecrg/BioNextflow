process GENE_COUNTS_STATS {

    container 'biocorecrg/multiqc:1.28'
    label     'low'

    
    input:
    path(desc_file)         // desc.txt file 
    path(annotated_file)    // norm_counts.genes file output by multiq_pca as emit: norm_counts
    val(experiment)         // "rnaseq" or "smallrnaseq"
    val(rna_type)           // If rna_type = "" is set a default list for smallrnaseq and rnaseq

    output:
    path("*.csv")               , emit: csv          // file for multiqc to create RNA stats table
    path("versions.yml")        , emit: versions

    script:
    def rna_type_flag = rna_type != "" ? "-r ${rna_type}" : ""
    
    """
    RNA_summary.py --desc ${desc_file} --annotation ${annotated_file} --experiment ${experiment} ${rna_type_flag}
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        container: "${task.container}"
    END_VERSIONS
    """
}
