/*
*  QC module
*  This workflow allows to make QC on input data
*  It needs input fastq
*/

params.LABEL = ""
params.CONTAINER = "quay.io/biocontainers/fastqc:0.11.9--0"

process fastQC {
    tag { fastq }
    label (params.LABEL)
    container params.CONTAINER

    input:
    path(fastq)

    output:
    path("*_fastqc.*") 

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
