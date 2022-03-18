/*
*  QC module
*  This workflow allows to make QC on input data
*  It needs input fastq
*/

params.LABEL = ""
params.EXTRAPARS = ""
params.CONTAINER = "quay.io/biocontainers/fastq-screen:0.14.0--pl5321hdfd78af_2"
params.OUTPUT = ""

process fastqScreen {
    tag "${fastq}"
    label (params.LABEL)
    container params.CONTAINER
    if (params.OUTPUT != "") { publishDir(params.OUTPUT, mode:'copy') }

    input:
    path(conf)
    tuple val(id), path(fastq)

    output:
    tuple val(id), path("*_screen.*") 

    script:
	"""
	fastq_screen ${params.EXTRAPARS} --conf ${conf} --threads ${task.cpus} ${fastq} 
	"""
}

workflow FASTQ_SCREEN {
    take: 
    conf
    fastq
    
    main:
    out = fastqScreen(conf, fastq)
    
    emit:
    out

}

process getVersion {
    container params.CONTAINER

    output:
	stdout emit: out    
    
    shell:
    """
    fastq_screen --version
    """
}



workflow GET_VERSION {
    main:
		getVersion()
    emit:
    	getVersion.out
}   
 
