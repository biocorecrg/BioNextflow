/*
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
params.CONTAINER = "quay.io/biocontainers/bismark_0.24.2hdfd78af_0"

process getVersion {
    container params.CONTAINER

    output:
	stdout emit: out    
    
    shell:
    """
		bismark --version|grep Version
    """
}


process index {
    label (params.LABEL)
    tag { "${reference}" }
    container params.CONTAINER

    input:
    path(reference)

    output:
    path("Bisulfite_Genome")
    
    """
    bismark_genome_preparation --single_fasta ./
    """
}

process map {
    label (params.LABEL)
    tag { pair_id }
    container params.CONTAINER
    if (params.OUTPUT != "") { publishDir(params.OUTPUT, mode:'copy') }

    input:
    tuple val(pair_id), path(reads), path(index), path(genome)

    output:
    tuple val(pair_id), path("*.bam") 
    
	script:
    def cmd = "${reads[0]}"
	if (reads[1]) {
		cmd = "-1 ${reads[0]} -2 ${reads[1]}"
	}

    """
    bismark ${params.EXTRAPARS} --genome ./ ${cmd} --bam --parallel ${task.cpus}
    """
}


workflow MAP {
    take: 
    input
    
    main:
		out = map(input)
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
		out = MAP(input.combine(index).combine(reference))
    emit:
    	out
}

workflow GET_VERSION {
    main:
		getVersion()
    emit:
    	getVersion.out
}    
