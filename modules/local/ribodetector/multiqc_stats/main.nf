process STATS {

    container 'quay.io/biocontainers/ribodetector:0.3.1--pyhdfd78af_0'
    label     'low'

    
    input:
    tuple val(meta), path(reports)         
 
    output:
    tuple val(meta), path("ribo_stats_mqc.txt")     , emit: log         


    script:
    """
    echo '# id: ribodetector
    # plot_type: bargraph
    # section_name: Ribosome contamination
    # description: Percent of ribosomal reads 
    Filename	reads' > ribo_stats_mqc.txt;
	    awk '{print \$1"\t"\$2}' ${reports} >> ribo_stats_mqc.txt
    """
}
