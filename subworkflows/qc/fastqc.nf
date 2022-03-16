/*
*  QC module
*  This workflow allows to make QC on input data
*  It needs input fastq
*/

params.LABEL = ""
params.CONTAINER = "quay.io/biocontainers/fastqc:0.11.9--0"
params.OUTPUT = ""

process fastQC {
    tag "${fastq}"
    label (params.LABEL)
    container params.CONTAINER
    if (params.OUTPUT != "") { publishDir(params.OUTPUT, mode:'copy') }

    input:
    path(fastq)

    output:
    path("*_fastqc.*") 

    script:
	"""
	fastqc -t ${task.cpus} ${fastq} 
	"""
}

process fastQC2 {
    tag "${fastq}"
    label (params.LABEL)
    container params.CONTAINER
    if (params.OUTPUT != "") { publishDir(params.OUTPUT, mode:'copy') }

    input:
    tuple val(id), path(fastq)

    output:
    tuple val(id), path("*_fastqc.*") 

    script:
	"""
	fastqc -t ${task.cpus} ${fastq} 
	"""
}

workflow FASTQC {
    take: 
    fastq
    
    main:
    out = fastQC(fastq)
    
    emit:
    out

}

process getVersion {
    container params.CONTAINER

    output:
	stdout emit: out    
    
    shell:
    """
    fastqc --version
    """
}


workflow FASTQCP {
    take: 
    fastqp
    
    main:
    fastqp.map{
		[it[1]]
	}.flatten().set{fastq}	
    out = fastQC(fastq)
    
    emit:
    out

}

workflow FASTQC_ID {
    take: 
    fastqp
    
    main:
    fastqp.map{
		[it[0], it[1]]
	}.set{fastq}	
    out = fastQC2(fastq)
    
    emit:
    out

}


workflow GET_VERSION {
    main:
		getVersion()
    emit:
    	getVersion.out
}   
 
