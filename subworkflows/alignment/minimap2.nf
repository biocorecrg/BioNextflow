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
params.EXTRAPARS = ""
params.OUTPUT = ""
params.CONTAINER = "quay.io/biocontainers/mulled-v2-66534bcbb7031a148b13e2ad42583020b9cd25c4:7e6194c85b2f194e301c71cdda1c002a754e8cc1-0"
params.INCLUDEUNMAP = ""

process getVersion {
    container params.CONTAINER

    output:
	stdout emit: out    
    
    shell:
    """
    echo "minimap2 "`minimap2 --version`
    """
}

process map {
    label (params.LABEL)
    tag "${idfile}"
    container params.CONTAINER
    if (params.OUTPUT != "") { publishDir(params.OUTPUT, mode:'copy') }

    input:
    tuple val(idfile), path(fastq_file)
	path(reference)  

    output:
    tuple val(idfile), path("${idfile}.bam")
    
    script:
    if (params.INCLUDEUNMAP=="")
    """
        minimap2 -t ${task.cpus} -a ${params.EXTRAPARS} ${reference} ${fastq_file} | samtools view -@ ${task.cpus} -F4 -hSb - > ${idfile}.bam
    """
    else
    """
        minimap2 -t ${task.cpus} -a ${params.EXTRAPARS} ${reference} ${fastq_file} | samtools view -@ ${task.cpus} -hSb - > ${idfile}.bam
    """
}


workflow MAP {
    take: 
    input
    reference
    
    main:
		map(input, reference)
	emit:
    	bams = map.out
	
}


workflow GET_VERSION {
    main:
		getVersion()
    emit:
    	getVersion.out
}    
