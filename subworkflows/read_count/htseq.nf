/*
*  nanoCount module
*  This workflow allows to make count on input data
*  It needs input fastq
*/

params.LABEL = ""
params.CONTAINER = "quay.io/biocontainers/htseq:0.13.5--py39h70b41aa_1"
params.OUTPUT = ""
params.EXTRAPARS = ""

include { unzipNamedPipe } from '../global_functions.nf'



process HtseqCount {
    tag "${id}"
    label (params.LABEL)
    container params.CONTAINER
    if (params.OUTPUT != "") { publishDir(params.OUTPUT, pattern:'*.counts', mode:'copy') }
   
   
    input:
    path(annotation_file)
    tuple val(id), path(bamfile)
    val(doanno)

    output:
    tuple val(id), path("${id}.counts"), emit: counts
    tuple val(id), path("${id}_anno.bam"), emit: bam optional true
  
	script:    
	def anno = unzipNamedPipe(annotation_file)
	def annopar = ""
	if (doanno=="yes") {
		annopar = "-p bam -o ${id}_anno.bam"
	} 
	"""
	    htseq-count ${annopar} ${params.EXTRAPARS} -n ${task.cpus} ${bamfile} ${anno} > ${id}.counts
	"""

}

workflow COUNT {
    take: 
    annotation
    fastq
    
    main:
	anno_file = file(annotation)
	if( !anno_file.exists() ) exit 1, "Missing ${annotation} file!"
    out = HtseqCount(annotation, fastq, "no")
    
    emit:
	counts = out.counts

}

workflow COUNT_AND_ANNO {
    take: 
    annotation
    fastq
    
    main:
	anno_file = file(annotation)
	if( !anno_file.exists() ) exit 1, "Missing ${annotation} file!"
    out = HtseqCount(annotation, fastq, "yes")
    
    emit:
	counts = out.counts
	bam = out.bam

}

process getVersion {
    container params.CONTAINER

    output:
	stdout emit: out    
    
    shell:
    """
    htseq-count --version
    """
}


workflow GET_VERSION {
    main:
		getVersion()
    emit:
    	getVersion.out
}   
 
