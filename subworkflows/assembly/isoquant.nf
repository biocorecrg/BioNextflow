/*
* isoquant
*/

params.LABEL = ""
params.EXTRAPARS = ""
params.OUTPUT = ""
params.OUTPUTMODE = "copy"
params.CONTAINER = "quay.io/biocontainers/isoquant:2.0.0--hdfd78af_0"

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
	def unzipGen  	 = unzipCmd(genome)
	def genome_name  = unzipGen[0]
	def cmd_g_unzip  = unzipGen[1]
	def cmd_g_clean  = unzipGen[2]

	def unzipAnno  		 = unzipCmd(annotation)
	def annotation_name  = unzipAnno[0]
	def cmd_a_unzip  	 = unzipAnno[1]
	def cmd_a_clean  	 = unzipAnno[2]

	"""
	${cmd_g_unzip}
	${cmd_a_unzip}
	isoquant.py -r ${genome_name} --threads ${task.cpus} ${params.EXTRAPARS} --genedb ${annotation_name} --bam ${bamfiles}  -o ${sampleID}_isoquant
	${cmd_a_clean}
	${cmd_g_clean}
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
