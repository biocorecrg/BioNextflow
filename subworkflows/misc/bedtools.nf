/*
*  Bedtools 
*/

params.LABEL = ""
params.EXTRAPARS = ""

params.OUTPUT = ""
params.CONTAINER = "quay.io/biocontainers/bedtools:2.30.0--hc088bd4_0"

include { unzipCmd } from '../global_functions.nf'
include { PossiblyUnzipGenome } from '../misc/misc.nf'

process getVersion {
    container params.CONTAINER

    output:
	stdout emit: out    
    
    shell:
    """
	bedtools --version
    """
}

process shuffleBed {
    label (params.LABEL)
    tag { id }
    container params.CONTAINER
    if (params.OUTPUT != "") { publishDir(params.OUTPUT, mode:'copy') }

    input:
    tuple val(id), path(bed)
	path(fasta) 
	
    output:
	tuple val(id), path("${id}_shuffle.bed")
    
	script:
    """
    bedtools shuffle -i ${bed} ${params.EXTRAPARS} -g ${fasta} > ${id}_shuffle.bed
    """
}

process bedToFasta {
    label (params.LABEL)
    tag { id }
    container params.CONTAINER
    if (params.OUTPUT != "") { publishDir(params.OUTPUT, mode:'copy') }

    input:
    tuple val(id), path(bed)
	path(fasta) 
	
    output:
	tuple val(id), path("${id}.fasta")
    
	script:
    """
	bedtools getfasta -fi ${fasta} -bed ${bed} ${params.EXTRAPARS} > ${id}.fasta
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

process multiCov {
    label (params.LABEL)
    tag { id }
    container params.CONTAINER
    if (params.OUTPUT != "") { publishDir(params.OUTPUT, mode:'copy') }

    input:
    tuple val(id), path(bed), path(bams), path(bais)

    output:
	tuple val(id), path("${id}.count")
    
	script:
	def bam_list = bams.join(" ")
	def bam_header = bams.join("\t")

    """
    echo "#chr	start	end	${bam_header}" > ${id}.count 
	bedtools multicov ${params.EXTRAPARS} -bams ${bam_list} -bed ${bed} >> ${id}.count
    """

}

process unionBedG {
    label (params.LABEL)
    tag { id }

    container params.CONTAINER
    if (params.OUTPUT != "") { publishDir(params.OUTPUT, mode:'copy') }

    input:
    tuple val(id), path(bed_files), path(genome_index)

    output:
	path("${id}.union.bed")
    
	script:
	def bed_list = bed_files.join(" ")

    """
	bedtools unionbedg ${params.EXTRAPARS} -empty -g ${genome_index}  -names ${bed_list} -i ${bed_list} > ${id}.union.bed
    """

}

workflow BEDTOFASTA {
    take: 
    bed
    fasta
    
    main:
        fasta_ch = PossiblyUnzipGenome(fasta)
		out = bedToFasta(bed, fasta_ch)
    emit:
    	out
}

workflow SHUFFLEBED {
    take: 
    bed
    fasta
    
    main:
        out = shuffleBed(bed, fasta)
    emit:
    	out
}



workflow UNIONBEDG {
    take: 
    input
    genomeindex
    
    main:
		out = unionBedG(input.combine(genomeindex))
    emit:
    	out
}


workflow MULTIINTER {
    take: 
    input
    
    main:
		out = multiInter(input)
    emit:
    	out
}

workflow SORT {
    take: 
    input
    
    main:
		out = sortBed(input)
    emit:
    	out
}

workflow MERGE {
    take: 
    input
    
    main:
		out = mergeBed(input)
    emit:
    	out
}

workflow MERGE_MULTI {
    take: 
    input
    
    main:
    	concat = concatBed(input)
		out = mergeBed(concat)
    emit:
    	out
}

workflow MULTICOV {
    take: 
    intervals
    bams
    bais
    
    main:
    	bams_bais_coll = bams.toList().map{[it]}.combine(bais.toList().map{[it]})
        comb = intervals.combine(bams_bais_coll)
    	out = multiCov(comb)

    emit:
    	out
}


workflow GET_VERSION {
    main:
		getVersion()
    emit:
    	getVersion.out
}    
