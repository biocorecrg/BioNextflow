/*
* Epinano
*/

params.LABEL = ""
params.EXTRAPARS = ""
params.OUTPUT = ""
params.OUTPUTMODE = "copy"
params.CONTAINER = "biocorecrg/mopexpress:0.1"

include { unzipCmd } from '../global_functions.nf'

/*
*/
process getVersion {
    container params.CONTAINER

    output:
	stdout emit: out    
    
    shell:
    """
    	echo bambu' '`Rscript -e "library('bambu'); packageVersion('bambu')"` 2>/dev/null
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
    tuple val(sampleID), path(bamfiles)
    
    output:
    tuple val(sampleID), path("${sampleID}_bambu")
    
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
	R --vanilla --slave -e "library(bambu)
bamfiles <- list.files(path='./', full.names = TRUE, pattern = '*.bam')
bamFiles <- Rsamtools::BamFileList(bamfiles)
bambuAnnotations <- prepareAnnotations(\'./${annotation_name}\')
Rsamtools::indexFa(\'./${genome_name}\')
se <- bambu(reads = bamFiles, annotations = bambuAnnotations, genome = \'./${genome_name}\', \
ncore = ${task.cpus} ${params.EXTRAPARS})
writeBambuOutput(se, \'${sampleID}_bambu\')"
	${cmd_g_clean}
	${cmd_a_clean}
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
