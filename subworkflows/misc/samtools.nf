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
    if (params.OUTPUT != "") { publishDir(params.OUTPUT, mode:'copy') }

    input:
    tuple val(pair_id), path(reads)

    output:
    tuple val(pair_id), path("${pair_id}_s.bam") 
    
	script:
    """    
    samtools sort -t ${task.cpus} ${params.EXTRAPARS} -o ${pair_id}_s.bam  ${reads}
    """
}

process indexBam {
    label (params.LABEL)
    tag { pair_id }
    container params.CONTAINER
    if (params.OUTPUT != "") { publishDir(params.OUTPUT, mode:'copy') }

    input:
    tuple val(pair_id), path(reads)

    output:
    tuple val(pair_id), path("*.bai") 
    
	script:
    """    
    samtools index ${params.EXTRAPARS} ${reads}
    """
}

process viewBam {
    label (params.LABEL)
    tag { pair_id }
    container params.CONTAINER
    if (params.OUTPUT != "") { publishDir(params.OUTPUT, mode:'copy') }

    input:
    tuple val(pair_id), path(reads)

    output:
    tuple val(pair_id), path("${pair_id}_f.bam") 
    
	script:
    """    
	samtools view ${params.EXTRAPARS} ${reads} > ${pair_id}_f.bam
    """
}

workflow INDEX {
    take: 
    reads
    
    main:
		out = indexBam(reads)
    emit:
    	out
}

workflow SORT {
    take: 
    reads
    
    main:
		out = sortAln(reads)
    emit:
    	out
}

workflow BVIEW {
    take: 
    reads
    
    main:
		out = viewBam(reads)
    emit:
    	out
}




workflow GET_VERSION {
    main:
		getVersion()
    emit:
    	getVersion.out
}    
