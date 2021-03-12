/*
*  Bedtools 
*/

params.LABEL = ""
params.EXTRAPARS = ""

params.OUTPUT = ""
params.CONTAINER = "quay.io/biocontainers/bedtools:2.30.0--hc088bd4_0"

include { unzipCmd } from '../global_functions.nf'

process getVersion {
    container params.CONTAINER

    output:
	stdout emit: out    
    
    shell:
    """
	bedtools --version
    """
}


process sortBed {
    label (params.LABEL)
    tag { id }
    container params.CONTAINER
    if (params.OUTPUT != "") { publishDir(params.OUTPUT, mode:'copy') }

    input:
    tuple val(id), path(input)

    output:
	tuple val(id), path("${id}_sorted.bed")
    
	script:
    """
	bedtools sort -i ${input} ${params.EXTRAPARS} > ${id}_sorted.bed
    """
}

process mergeBed {
    label (params.LABEL)
    tag { id }
    container params.CONTAINER
    if (params.OUTPUT != "") { publishDir(params.OUTPUT, mode:'copy') }

    input:
    tuple val(id), path(input)

    output:
	tuple val(id), path("${id}_merged.bed")
    
	script:
    """
	bedtools sort -i ${input} | bedtools merge ${params.EXTRAPARS} -i - > ${id}_merged.bed
    """
}

process concatBed {
    label (params.LABEL)
    tag { id }
    container params.CONTAINER

    input:
    tuple val(id), path(input)

    output:
	tuple val(id), path("${id}_concat.bed")
    
	script:
    """
	cat ${input} | awk '{OFS="\t"; print \$1,\$2,\$3,\$4,\$5,\$6}' > ${id}_concat.bed
    """
}

process multiInter {
    label (params.LABEL)
    tag { id }
    container params.CONTAINER
    if (params.OUTPUT != "") { publishDir(params.OUTPUT, mode:'copy') }

    input:
    tuple val(id), path(input)

    output:
	tuple val(id), path("${id}_multiinter.bed")
    
	script:
    """
	bedtools multiinter ${params.EXTRAPARS} -header -i ${input} > ${id}_multiinter.bed
    """

}

workflow BEDTOOLS_MULTIINTER {
    take: 
    input
    
    main:
		out = multiInter(input)
    emit:
    	out
}

workflow BEDTOOLS_SORT {
    take: 
    input
    
    main:
		out = sortBed(input)
    emit:
    	out
}

workflow BEDTOOLS_MERGE {
    take: 
    input
    
    main:
		out = mergeBed(input)
    emit:
    	out
}

workflow BEDTOOLS_MERGE_MULTI {
    take: 
    input
    
    main:
    	concat = concatBed(input)
		out = mergeBed(concat)
    emit:
    	out
}



workflow GET_VERSION {
    main:
		getVersion()
    emit:
    	getVersion.out
}    
