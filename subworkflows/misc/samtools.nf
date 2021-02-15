/*
*  samtools 
*/

params.LABEL = ""
params.EXTRAPARS = ""

params.OUTPUT = "samtools_out"
params.CONTAINER = "quay.io/biocontainers/mulled-v2-8a9a988fff4785176b70ce7d14ff00adccf8a5b8:aeac8200e5c50c5acf4dd14792fd8453255af835-0"

process getVersion {
    container params.CONTAINER

    output:
	stdout emit: out    
    
    shell:
    """
    samtools --version | grep samtools
    """
}




process sortAln {
    label (params.LABEL)
    tag { pair_id }
    container params.CONTAINER
    publishDir(params.OUTPUT, mode:'copy')

    input:
    tuple val(pair_id), path(reads)

    output:
    tuple val(pair_id), path("${pair_id}_s.bam") 
    
	script:
    """    
    samtools sort -t ${task.cpus} ${params.EXTRAPARS} -o ${pair_id}_s.bam  ${reads}
    """
}

workflow SAMTOOLS_SORT {
    take: 
    reads
    
    main:
		out = sortAln(reads)
    emit:
    	out
}




workflow GET_VERSION {
    main:
		getVersion()
    emit:
    	getVersion.out
}    
