/*
* isoquant
*/

params.LABEL = ""
params.EXTRAPARS = ""
params.OUTPUT = ""
params.OUTPUTMODE = "copy"
params.CONTAINER = "quay.io/biocontainers/isoquant:3.6.3--hdfd78af_0"

include { unzipCmd } from '../global_functions.nf'

/*
*/
process getVersion {
    container params.CONTAINER

    output:
	stdout emit: out    
    
    shell:
    """
    	
    """
}


/*
*/

process assembleTranscripts {
    if (params.OUTPUT != "") { publishDir(params.OUTPUT, mode:'copy') }

    container params.CONTAINER
    label (params.LABEL)
    tag "${sampleID}" 
 	
    input:
    path(genome)
    path(annotation)
    tuple val(sampleID), path(bamfiles), path(indexes)
    
    output:
    tuple val(sampleID), path("${sampleID}_isoquant")
    
    script:

	"""
	isoquant.py -r ${genome} --threads ${task.cpus} ${params.EXTRAPARS} --genedb ${annotation} --bam ${bamfiles}  -o ${sampleID}_isoquant
	"""
}

/*
*/

workflow ASSEMBLE {

    take: 
	genome
    annotation
    bamfiles
    
    main:
    	out = assembleTranscripts(genome, annotation, bamfiles)

	emit:
    	out
}

/*
*/

workflow GET_VERSION {
    main:
		getVersion()
    emit:
    	getVersion.out
}    
