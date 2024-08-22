/*
*/

params.LABEL = ""
params.EXTRAPARS = ""
params.OUTPUT = ""
params.CONTAINER = "quay.io/biocontainers/nanoq:0.10.0--h031d066_2"


process getVersion {
    container params.CONTAINER

    output:
	stdout emit: out    
    
    shell:
    """
		nanoq -V    
    """
}


process filter {
    tag { idfile }
    label (params.LABEL)
    if (params.OUTPUT != "") {publishDir(params.OUTPUT, mode: 'copy') }

    container params.CONTAINER
             
    input:
    tuple val(idfile), path(fastq)
 
    
    output:
    tuple val(idfile), path("*-filt.fastq*")

    script:
    def output = "${idfile}-filt.fastq.gz"
	"""
		nanoq -i ${fastq} ${params.EXTRAPARS} -O g -o ${output}
	"""

}

process report {
    tag { idfile }
    label (params.LABEL)
    if (params.OUTPUT != "") {publishDir(params.OUTPUT, mode: 'copy') }

    container params.CONTAINER
             
    input:
    tuple val(idfile), path(fastq)
 
    
    output:
    tuple val(idfile), path("*-report.txt")

    script:

	"""
		nanoq -i ${fastq} -s ${params.EXTRAPARS} -r ${idfile}-report.txt
	"""

}


 workflow REPORT {

    take: 
    input_fastq
    
    main:
     out = report(input_fastq)

	emit:
		out
		
}

 workflow FILTER {
 
    take: 
    input_fastq
    
    main:
     out = filter(input_fastq)

	emit:
		out
		
}

workflow GET_VERSION {
    main:
		getVersion()
    emit:
    	getVersion.out
} 


