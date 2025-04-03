/*
*  QC module
*  This workflow allows to make QC on input data
*  It needs input fastq
*/

params.LABEL = ""
params.CONTAINER = "quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0"
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
    def memory = task.memory.toMega() > 10000 ? 10000 : task.memory.toMega()

    """
	fastqc --dir ./ --memory ${memory} -t ${task.cpus} ${fastq} 
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
    fastqp.transpose().set{fastq}	
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
 
