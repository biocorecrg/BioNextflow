/*
*  picard
*/

params.LABEL = ""
params.EXTRAPARS = ""

params.OUTPUT = ""
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
    if (params.OUTPUT != "") { publishDir(params.OUTPUT, mode:'copy') }

    input:
    tuple val(pair_id), path(reads)

    output:
    tuple val(pair_id), path("${pair_id}_s.bam") 
    
	script:

    """    
	picard -Xmx${task.memory.giga}g SortSam I=${reads} TMP_DIR=`pwd`/tmp O=${pair_id}_s.bam SORT_ORDER=coordinate
	rm -fr tmp
    """
}

process sortSamName {
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
	picard -Xmx${task.memory.giga}g SortSam I=${reads} TMP_DIR=`pwd`/tmp O=${pair_id}_s.bam SORT_ORDER=queryname
	rm -fr tmp
    """
}

process markDuplicates {
    label (params.LABEL)
    tag { pair_id }
    container params.CONTAINER
    if (params.OUTPUT != "") { publishDir(params.OUTPUT, mode:'copy') }

    input:
    tuple val(pair_id), path(reads)
	val(remove)

    output:
    tuple val(pair_id), path("${pair_id}_dedup.bam") 
    

	script:
    def remove_cmd = ""
    if (remove == "remove") {
    	remove_cmd = "REMOVE_SEQUENCING_DUPLICATES=TRUE"
    }     
    """    
	picard -Xmx${task.memory.giga}g MarkDuplicates ${params.EXTRAPARS} TMP_DIR=`pwd`/tmp ${remove_cmd} I=${reads} O=${pair_id}_dedup.bam M=${pair_id}.dupmet.txt 
    """
}

workflow SORT_COORD {
    take: 
    input
    
    main:
		out = sortSamCoord(input)
    emit:
    	out
}

workflow SORT_NAME {
    take: 
    input
    
    main:
		out = sortSamName(input)
    emit:
    	out
}

workflow REM_DUP {
    take: 
    input
    
    main:
		out = markDuplicates(input, "remove")
    emit:
    	out
}

workflow MARK_DUP {
    take: 
    input
    
    main:
		out =  markDuplicates(input, "no")
    emit:
    	out
}



workflow GET_VERSION {
    main:
		getVersion()
    emit:
    	getVersion.out
}    
