/*
* bwa subworkflows (samtools embedded for converting the output) 
* The accessible subworkflows are:
* GET_VERSION that emits the version of bwa and samtools as stdout
* INDEX that takes:
*	a channel with an optionally gzipped fasta file
*   it emits a list of files as index
* ALIGN that takes:
*	a channel list with index files as produced by BWA_INDEX
*	a channel containing a tuple with id, and one or two (gzipped) fastq files 
*	it emits a channel containing a tuple of id, bam file
* ALL (INDEX + MAP) that takes:
*	a channel with an optionally gzipped fasta file
*	a channel containing a tuple with id, and one or two (gzipped) fastq files 
*   it emits a channel containing a tuple of id, bam file
*
* The parameters are: 
*	LABEL that allows connecting labels specified in nextflow.config with the subworkflows
*	EXTRAPARS only for mapping step for adding custom command line parameters for bwa
*	OUTPUT for storing the final sub-workflow output 
*	CONTAINER that can be eventually overridden for feeding a custom container from the main.nf file
*/

params.LABEL = ""
params.EXTRAPARS = ""
params.OUTPUT = ""
params.CONTAINER = "quay.io/biocontainers/mulled-v2-8a9a988fff4785176b70ce7d14ff00adccf8a5b8:aeac8200e5c50c5acf4dd14792fd8453255af835-0"
params.STOREINDEX = ""

process getVersion {
    container params.CONTAINER

    output:
	stdout emit: out    
    
    shell:
    """
    bwa 2>&1| grep Version | awk '{print "bwa "\$0}';
    samtools --version | grep samtools
    """
}


process index {
    label (params.LABEL)
    tag { "${reference}" }
    container params.CONTAINER
    if (params.STOREINDEX != "") { storeDir(params.STOREINDEX) }
 

    input:
    path(reference)

    output:
    path("${reference}.*")
    
    """
    bwa index ${reference}
    """
}

process map {
    label (params.LABEL)
    tag { pair_id }
    container params.CONTAINER
    if (params.OUTPUT != "") { publishDir(params.OUTPUT, mode:'copy') }

    input:
    tuple val(pair_id), path(reads)
    path(indexes)

    output:
    tuple val(pair_id), path("${pair_id}.bam") 
    
	script:
    def indexname = indexes[0].baseName

    """    
    bwa mem -t ${task.cpus} ${params.EXTRAPARS} ${indexname} ${reads} | samtools view -@ ${task.cpus} -Sb > ${pair_id}.bam
    """
}

process map_with_params {

    label (params.LABEL)
    tag { pair_id }
    container params.CONTAINER
    if (params.OUTPUT != "") { publishDir(params.OUTPUT, mode:'copy') }

    input:
    tuple val(pair_id), path(reads), val(extrapars)
    path(indexes)

    output:
    tuple val(pair_id), path("${pair_id}.bam") 
    
	script:
    def indexname = indexes[0].baseName

    """    
    bwa mem -t ${task.cpus} ${params.EXTRAPARS} ${extrapars} ${indexname} ${reads} | samtools view -@ ${task.cpus} -Sb > ${pair_id}.bam
    """
}

workflow MAP {
    take: 
    input
    indexes
    
    main:
		out = map(input, indexes)
    emit:
    	out
}

workflow INDEX {
    take: 
    reference
    
    main:
		ref_file = file(reference)
		if( !ref_file.exists() ) exit 1, "Missing ${reference} file!"
		out = index(reference)

    emit:
    	out
}

workflow ALL {
    take: 
    reference
    input
    
    main:
		index = INDEX(reference)
		out = MAP(input, index)
    emit:
    	out
}

workflow GET_VERSION {
    main:
		getVersion()
    emit:
    	getVersion.out
}    
