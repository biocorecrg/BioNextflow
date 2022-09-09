/*
*  QC module
*  This workflow allows to make QC on input data
*  It needs input fastq
*/

params.LABEL = ""
params.CONTAINER = "biocorecrg/falco:1.1.0"
params.OUTPUT = ""
params.EXTRAPARS = ""

process falco {

    tag "${fastq}"
    label (params.LABEL)
    container params.CONTAINER
    if (params.OUTPUT != "") { publishDir(params.OUTPUT, mode:'copy') }

    input:
    tuple val(id), path(fastq)

    output:
    tuple val(id), path("*_fastqc") 

    script:
    """
	falco -o ${id}_fastqc -t ${task.cpus} ${params.EXTRAPARS} ${fastq} 
    """
}


process getVersion {
    container params.CONTAINER

    output:
	stdout emit: out    
    
    shell:
    """
    falco --version
    """
}



workflow FALCOQC {
    take: 
    fastq
    
    main:
    out = falco(fastq)
    
    emit:
    out

}


workflow GET_VERSION {
    main:
		getVersion()
    emit:
    	getVersion.out
}   
 
