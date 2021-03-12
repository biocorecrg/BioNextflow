/*
* samblaster subworkflows (samtools and bwa embedded)
* It requires the subworkflow bwa index 
* The accessible subworkflows are:
* GET_VERSION that emits the version of bwa, samtools and samblaster as stdout
* SAMBLASTER_MAP that takes:
*	a channel list with index files as produced by BWA_INDEX
*	a channel containing a tuple with id and one or two (gzipped) fastq files
*	it emits a channel containing a tuple of id, bam file
* SAMBLASTER_ALL (BWA_INDEX + SAMBLASTER_MAP) that takes:
*	a channel with an optionally gzipped fasta file
*   a channel containing one or two (gzipped) fastq files
*   it emits a channel containing a tuple of id, bam file
* The parameters are: 
*	LABEL that allows connecting labels specified in nextflow.config with the subworkflows
*	EXTRAPARS only for mapping step for adding custom command line parameters for bwa
*	OUTPUT for storing the final sub-workflow output 
*	CONTAINER that can be eventually overridden for feeding a custom container from the main.nf file
*/


include { BWA_INDEX as BWA_INDEX } from "./bwa"

params.LABEL = ""
params.EXTRAPARS = ""

params.OUTPUT = "bwa_sv_out"
params.CONTAINER = "quay.io/biocontainers/mulled-v2-8a9a988fff4785176b70ce7d14ff00adccf8a5b8:aeac8200e5c50c5acf4dd14792fd8453255af835-0"

process getVersion {
    container params.CONTAINER

    output:
	stdout emit: out    
    
    shell:
    """
    bwa 2>&1| grep Version | awk '{print "bwa "\$0}';
    samblaster -h 2>&1 | grep Version;
    samtools --version | grep samtools
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
	if (reads[1]) {
	    """
	    bwa mem -R "@RG\\tID:id\\tSM:sample\\tLB:lib" -t ${task.cpus} ${indexname} ${reads} \
	    | samblaster --excludeDups --addMateTags --maxSplitCount 2 --minNonOverlap 20 \
	    | samtools view  -@ ${task.cpus} -S -b - > ${pair_id}.bam;
	    """
	} else {
    	""" 
   		bwa mem -R "@RG\\tID:id\\tSM:sample\\tLB:lib" -t ${task.cpus} ${indexname} ${reads} \
    	| samblaster --excludeDups --ignoreUnmated --maxSplitCount 2 --minNonOverlap 20 \
    	| samtools view  -@ ${task.cpus} -S -b - > ${pair_id}.bam;
    	"""		
	}
}

workflow SAMBLASTER_MAP {
    take: 
    input
    indexes
    
    main:
		out = map(input, indexes)
    emit:
    	out
}


workflow SAMBLASTER_ALL {
    take: 
    reference
    input
    
    main:
		index = BWA_INDEX(reference)
		out = SAMBLASTER_MAP(input, index)
    emit:
    	out
}


workflow GET_VERSION {
    main:
		getVersion()
    emit:
    	getVersion.out
}    
