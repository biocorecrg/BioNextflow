/*
* bwa subworkflows (samtools embedded for converting the output) 
* The accessible subworkflows are:
* GET_VERSION that emits the version of bwa and samtools as stdout
* INDEX that takes:
*	a channel with an optionally gzipped fasta file
*   it emits a list of files as index
* ALIGN that takes:
*	a channel list with index files as produced by INDEX
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
params.CONTAINER = "quay.io/biocontainers/biscuit:1.4.0.20240108--h0be9327_0"

process getVersion {
    container params.CONTAINER

    output:
	stdout emit: out    
    
    shell:
    """
		biscuit --version
    """
}


process index {
    label (params.LABEL)
    tag { "${reference}" }
    container params.CONTAINER

    input:
    path(reference)

    output:
    path("${reference}.*")
    
    """
    biscuit index ${reference}
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
    def indexname = indexes[0].baseName - ~/\.\w+$/

    """
    biscuit align -@ ${task.cpus} ${params.EXTRAPARS} ${indexname} ${reads} > ${pair_id}.bam
    """
}

process pileup {

    label (params.LABEL)
    tag { id }
    container params.CONTAINER
    if (params.OUTPUT != "") { publishDir(params.OUTPUT, mode:'copy') }

    input:
    tuple val(id), path(bam), path(bai), path(reference), path(reffai)

    output:
    tuple val(id), path("${id}.vcf") 
    
	script:

    """
    biscuit pileup -@ ${task.cpus} ${params.EXTRAPARS} ${reference} ${bam} -o ${id}.vcf
    """

}

workflow PILEUP {

    take: 
    input
    reference
    
    main:
		out = pileup(input.combine(reference))
		
    emit:
    	out
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
