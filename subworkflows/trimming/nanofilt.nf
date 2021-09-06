/*
*/

params.LABEL = ""
params.EXTRAPARS = ""
params.OUTPUT = ""
params.CONTAINER = "quay.io/biocontainers/nanofilt:2.8.0--py_0"

include { zcatOrCat } from '../global_functions.nf'

process getVersion {
    container params.CONTAINER

    output:
	stdout emit: out    
    
    shell:
    """
		NanoFilt --version       
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
    def cmd = zcatOrCat(fastq)
    def output = "${idfile}-filt.fastq.gz"
    
	"""
		${cmd} | NanoFilt ${params.EXTRAPARS} | gzip > ${output}
	"""

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


