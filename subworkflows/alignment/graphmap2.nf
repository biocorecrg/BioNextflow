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
params.CONTAINER = "biocorecrg/mopprepr:0.7"

process getVersion {
    container params.CONTAINER

    output:
	stdout emit: out    
    
    shell:
    """
    graphmap2 2>&1| grep Versio| awk '{print "graphmap2 "\$2}'
    """
}

process map {
    label (params.LABEL)
    tag "${idfile}"
    container params.CONTAINER

    input:
    tuple val(idfile), path(fastq_file)
	path(reference)  

    output:
    tuple val(idfile), path("${idfile}.bam")
    
    script:
    """
        graphmap2 align -t ${task.cpus} ${params.EXTRAPARS} -r ${reference} -d ${fastq_file}  | samtools view -@ ${task.cpus} -F4 -hSb - > ${idfile}.bam
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
