/*
*  picard
*/

params.LABEL = ""
params.EXTRAPARS = ""

params.OUTPUT = "picard_out"
params.CONTAINER = "quay.io/biocontainers/picard:2.23.2--0"

process getVersion {
    container params.CONTAINER

    output:
	stdout emit: out    
    
    shell:
    """
    picard BamIndexStats 2>&1| grep Version | awk '{print "Picard\t"\$2}'
    """
}

process sortSamCoord {
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
	picard SortSam I=${reads} O=${pair_id}_s.bam SORT_ORDER=coordinate
    """
}

process removeDuplicates {
    label (params.LABEL)
    tag { pair_id }
    container params.CONTAINER
    publishDir(params.OUTPUT, mode:'copy')

    input:
    tuple val(pair_id), path(reads)

    output:
    tuple val(pair_id), path("${pair_id}_dedup.bam") 
    
	script:

    """    
	picard MarkDuplicates ${params.EXTRAPARS} REMOVE_SEQUENCING_DUPLICATES=TRUE I=${reads} O=${pair_id}_dedup.bam M=${pair_id}.dupmet.txt 
    """
}

workflow PICARD_SORT_COORD {
    take: 
    input
    
    main:
		out = sortSamCoord(input)
    emit:
    	out
}


workflow PICARD_REM_DUP {
    take: 
    input
    
    main:
		out = removeDuplicates(input)
    emit:
    	out
}



workflow GET_VERSION {
    main:
		getVersion()
    emit:
    	getVersion.out
}    
