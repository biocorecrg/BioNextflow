/*
*  nanoplot module
*  This workflow allows to make QC on input data
*  It needs input fastq
*/

params.LABEL = ""
params.CONTAINER = "quay.io/biocontainers/qualimap:2.2.2d--hdfd78af_2"
params.OUTPUT = ""

process qualimap_bamqc {
    tag { id }
    label (params.LABEL)
    container params.CONTAINER
    if (params.OUTPUT != "") { publishDir(params.OUTPUT, mode:'copy') }
   
    input:
    tuple val(id), path(bamfile)

    output:
    tuple val(id), path("*_stats"), emit: out optional true
    
    script:
    memory = task.memory.toGiga()-1
    """
    	qualimap bamqc --java-mem-size=${memory}G -nt ${task.cpus} -bam ${bamfile} 
    """
}



workflow BAM_QC {
    take: 
    bam
    
    main:
    out = qualimap_bamqc(bam)
    
    emit:
   	out = out

}

process getVersion {
    container params.CONTAINER

    output:
	stdout emit: out    
    
    shell:
    """
    qualimap -version | grep QualiMap
    """
}


workflow GET_VERSION {
    main:
		getVersion()
    emit:
    	getVersion.out
}   
 
