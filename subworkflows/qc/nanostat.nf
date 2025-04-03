/*
*  nanoplot module
*  This workflow allows to make QC on input data
*  It needs input fastq
*/

params.LABEL = ""
params.CONTAINER = "quay.io/biocontainers/nanostat:1.6.0--pyhdfd78af_0"
params.OUTPUT = ""

process nanoStat {
    tag { id }
    
    label (params.LABEL)
    container params.CONTAINER
    if (params.OUTPUT != "") { publishDir(params.OUTPUT, mode:'copy') }
   
    input:
    tuple val(id), path(bamfile)

    output:
    tuple val(id), path("${id}*")
    
    script:
    """
    NanoStat --bam ${bamfile} --tsv -o ./ -t ${task.cpus} -n ${id}_nanostat.txt
    """
}



workflow QC {
    take: 
    bam
    
    main:
    out = nanoStat(bam)
    
    emit:
   	out

}

process getVersion {
    container params.CONTAINER

    output:
	stdout emit: out    
    
    shell:
    """
    NanoPlot --version
    """
}


workflow GET_VERSION {
    main:
		getVersion()
    emit:
    	getVersion.out
}   
 
