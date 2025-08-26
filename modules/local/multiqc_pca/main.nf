process MULTIQC_PCA {
    
    container 'biocorecrg/deseq2:2022_2'

    input:
    tuple val(meta), path(input_files)
    path(desc)
    path(dgenes)
    val(wtype)
    val(extrapars)
	
    output:
    path("*_data.tsv"), emit: data
    path("*_variance.tsv"), emit: variance
	path("norm_counts.genes"), emit: norm_counts
	path("raw_counts.genes"), emit: raw_counts
	
	
    script:
    """
    multiqc_pca.R -type ${wtype} -dgenes ${dgenes} -desc ${desc}  ${extrapars} 
    """
}

