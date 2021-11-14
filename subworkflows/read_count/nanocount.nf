/*
*  nanoCount module
*  This workflow allows to make count on input data
*  It needs input fastq
*/

params.LABEL = ""
params.CONTAINER = "biocorecrg/nanocount:0.1"
params.OUTPUT = ""
params.EXTRAPARS = ""

process nanoCount {
    tag { id }
    label (params.LABEL)
    container params.CONTAINER
    if (params.OUTPUT != "") { publishDir(params.OUTPUT, mode:'copy') }
   
    input:
    tuple val(id), path(bamfile)

    output:
    tuple val(id), path("${id}.count")
    
	script:    
	"""
		NanoCount -i ${bamfile} ${params.EXTRAPARS} -o ${id}.count
	"""

}



workflow COUNT {
    take: 
    fastq
    
    main:
    out = nanoCount(fastq)
    
    emit:
	out

}

process getVersion {
    container params.CONTAINER

    output:
	stdout emit: out    
    
    shell:
    """
    NanoCount --version
    """
}


workflow GET_VERSION {
    main:
		getVersion()
    emit:
    	getVersion.out
}   
 
