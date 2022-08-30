/*
* Graphmap subworkflows 
*
* The parameters are: 
*	LABEL that allows connecting labels specified in nextflow.config with the subworkflows
*	EXTRAPARS only for mapping step for adding custom command line parameters for star
*	OUTPUT for storing the final sub-workflow output 
*	CONTAINER that can be eventually overridden for feeding a custom container from the main.nf file
*/

params.LABEL = ""
params.LABELINDEX = ""
params.EXTRAPARS = ""
params.EXTRAPARSINDEX = ""
params.OUTPUT = ""
params.CONTAINER = "biocorecrg/winnowmap:2.03"

process getVersion {
    container params.CONTAINER

    output:
	stdout emit: out    
    
    shell:
    """
    echo "winnowmap "`winnowmap --version`
    """
}

process index {
    label (params.LABEL)
    tag "${reference}"
    container params.CONTAINER

    input:
	path(reference)  

    output:
	path("repetitive_k.txt")
    
    script:
    """
		meryl count ${params.EXTRAPARSINDEX} output merylDB ${reference}
  		meryl print greater-than distinct=0.9998 merylDB > repetitive_k.txt	
    """
}

process map {
    label (params.LABEL)
    tag "${idfile}"
    container params.CONTAINER

    input:
    tuple val(idfile), path(fastq_file)
	path(reference)  
	path(repetitive)  

    output:
    tuple val(idfile), path("${idfile}.bam")
    
    script:
    """
    winnowmap -W ${repetitive} -t ${task.cpus} -a ${params.EXTRAPARS} ${reference} ${fastq_file} | samtools view -@ ${task.cpus} -F4 -hSb - > ${idfile}.bam
    """
}


workflow MAP {
    take: 
    input
    reference
    repetitive
    
    main:
		map(input, reference, repetitive)
	emit:
    	bams = map.out
	
}

workflow INDEX {
    take: 
    reference
    
    main:
		index(reference)
	emit:
    	repetitive = index.out
	
}

workflow ALL {
    take: 
    input
    reference
    
    main:
		repetitive = INDEX(reference)
		bams = MAP(input, reference, repetitive)
		
	emit:
    	bams = bams
	
}


workflow GET_VERSION {
    main:
		getVersion()
    emit:
    	getVersion.out
}    
