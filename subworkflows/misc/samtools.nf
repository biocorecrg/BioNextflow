/*
*  samtools 
*/

params.LABEL = ""
params.EXTRAPARS = ""

params.OUTPUT = ""
params.CONTAINER = "quay.io/biocontainers/samtools:1.16.1--h6899075_1"
params.OUTPUTMODE = "copy"

include { unzipCmd } from '../global_functions.nf'


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
    if (params.OUTPUT != "") { publishDir(params.OUTPUT, mode:params.OUTPUTMODE) }

    input:
    tuple val(pair_id), path(reads)

    output:
    tuple val(pair_id), path("${pair_id}_s.bam") 
    
	script:
    """    
    samtools sort -@ ${task.cpus} ${params.EXTRAPARS} -o ${pair_id}_s.bam  ${reads}
    """
}

process faidx {
    label (params.LABEL)
    tag { genome }
    container params.CONTAINER
    if (params.OUTPUT != "") { publishDir(params.OUTPUT, mode:params.OUTPUTMODE) }

    input:
    path(genome)

    output:
    path("*.fai") 
    
	script:
    def unzip_ref = unzipCmd(genome)
    def cmd_ref = unzip_ref[1]
    def ref_name = unzip_ref[0]
    def clean_ref = unzip_ref[2]
    """    
 		${cmd_ref} 
	    samtools faidx ${params.EXTRAPARS} ${ref_name}
	    ${clean_ref}
    """
}


process indexBam {
    label (params.LABEL)
    tag { pair_id }
    container params.CONTAINER
    if (params.OUTPUT != "") { publishDir(params.OUTPUT, mode:params.OUTPUTMODE) }

    input:
    tuple val(pair_id), path(reads)

    output:
    tuple val(pair_id), path("*.bai") 
    
	script:
    """    
    samtools index ${params.EXTRAPARS} ${reads}
    """
}

process catAln {
    label (params.LABEL)
    tag { pair_id }
    container params.CONTAINER
    if (params.OUTPUT != "") { publishDir(params.OUTPUT, mode:params.OUTPUTMODE) }

    input:
    tuple val(pair_id), path(reads)

    output:
    tuple val(pair_id), path("*_cat.bam") 
    
	script:
    """    
    samtools cat ${params.EXTRAPARS} -o ${pair_id}_cat.bam ${reads}
    """
}

process viewBam {
    label (params.LABEL)
    tag { pair_id }
    container params.CONTAINER
    if (params.OUTPUT != "") { publishDir(params.OUTPUT, mode:params.OUTPUTMODE) }

    input:
    tuple val(pair_id), path(reads)

    output:
    tuple val(pair_id), path("${pair_id}_f.bam") 
    
	script:
    """    
	samtools view -@ ${task.cpus} ${params.EXTRAPARS} ${reads} > ${pair_id}_f.bam
    """
}

process viewBam_two {
    label (params.LABEL)
    tag { pair_id }
    
    container params.CONTAINER
    if (params.OUTPUT != "") { publishDir(params.OUTPUT, mode:params.OUTPUTMODE) }

    input:
    tuple val(pair_id), path(reads), path(extrafile)

    output:
    tuple val(pair_id), path("${pair_id}_f.bam") 
    
	script:
    """    
	samtools view -@ ${task.cpus} ${params.EXTRAPARS} ${extrafile} ${reads} > ${pair_id}_f.bam
    """
}

process viewBam_exclude {
    label (params.LABEL)
    tag { pair_id }
    
    container params.CONTAINER
    if (params.OUTPUT != "") { publishDir(params.OUTPUT, mode:params.OUTPUTMODE) }
    
    input:
    tuple val(pair_id), path(reads), path(extrafile)
    
    output:
    tuple val(pair_id), path("${pair_id}_f.bam"), emit: filtered_bam
    tuple val(pair_id), path("${pair_id}_ex.bam"), emit: excluded_bam
    
        script:
    """    
        samtools view -@ ${task.cpus} -U ${pair_id}_ex.bam ${params.EXTRAPARS} ${extrafile} ${reads} > ${pair_id}_f.bam
    """ 
}  



process statBam {
    label (params.LABEL)
    tag { pair_id }
    container params.CONTAINER
    if (params.OUTPUT != "") { publishDir(params.OUTPUT, mode:params.OUTPUTMODE) }

    input:
    tuple val(pair_id), path(reads)

    output:
    tuple val(pair_id), path("${pair_id}.stat") 
    
	script:
    """    
	samtools flagstat -@ ${task.cpus} ${params.EXTRAPARS} ${reads} > ${pair_id}.stat
    """
}

process primaryCountFromStat {
    label (params.LABEL)
    tag { pair_id }
    container params.CONTAINER
    if (params.OUTPUT != "") { publishDir(params.OUTPUT, mode:params.OUTPUTMODE) }

    input:
    tuple val(pair_id), path(reads)

    output:
    tuple val(pair_id), stdout
    
	script:
    """    
	samtools flagstat -@ ${task.cpus} ${params.EXTRAPARS} ${reads} |  grep 'primary mapped' | cut -d " " -f 1 
    """
}

process primaryCountPercStat {
    label (params.LABEL)
    tag { pair_id }
    container params.CONTAINER
    if (params.OUTPUT != "") { publishDir(params.OUTPUT, mode:params.OUTPUTMODE) }

    input:
    tuple val(pair_id), path(reads)

    output:
    tuple val(pair_id), stdout
    
	script:
    """    
	samtools flagstat -@ ${task.cpus} ${params.EXTRAPARS} ${reads} | grep "primary mapped" | cut -d " " -f 6 | sed s/\\(//g |  sed s/\\%//g
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

workflow FAIDX {
    take: 
    genome
    
    main:
		out = faidx(genome)
    emit:
    	out
}

workflow STAT {
    take: 
    reads
    
    main:
		out = statBam(reads)
    emit:
    	out
}

workflow PRIM_COUNT {
    take: 
    reads
    
    main:
		out = primaryCountFromStat(reads)
    emit:
    	out
}

workflow PRIM_PERC {
    take: 
    reads
    
    main:
		out = primaryCountPercStat(reads)
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

workflow CAT {
    take: 
    reads
    
    main:
		out = catAln(reads)
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

workflow BVIEW2 {
    take: 
    reads
    
    main:
		out = viewBam_two(reads)
    emit:
    	out
}

workflow BVIEW_EXCLUDE {
   take:
   reads

   main:
	viewBam_exclude(reads)

   emit:

   filtered_bam = viewBam_exclude.out.filtered_bam
   excluded_bam = viewBam_exclude.out.excluded_bam


}



workflow GET_VERSION {
    main:
		getVersion()
    emit:
    	getVersion.out
}    
