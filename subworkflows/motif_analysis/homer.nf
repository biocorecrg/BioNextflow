/*
*  Homer 
*/

params.LABEL = ""
params.EXTRAPARS = ""

params.OUTPUT = ""
params.CONTAINER = "biocorecrg/chipanno:0.3"

include { unzipCmd } from '../global_functions.nf'

process getVersion {
    container params.CONTAINER

    output:
	stdout emit: out    
    
    shell:
    """
    echo 'homer version 4.9.1'
    """
}


process annotatePeaks {
    label (params.LABEL)
    tag { pair_id }
    container params.CONTAINER
    
    if (params.OUTPUT != "") { publishDir(params.OUTPUT, mode:'copy', pattern: '*.anno') }

    input:
    tuple val(pair_id), path(bedfile)
    path(annofile)

    output:
    tuple val(pair_id), path("*.anno") 
    
	script:
	def unzip      = unzipCmd(annofile)
	def file_name  = unzip[0]
	def cmd_unzip  = unzip[1]
	def cmd_clean  = unzip[2]
	"""
	${cmd_unzip}
	annotatePeaks.pl ${bedfile} none -gtf ${file_name} > ${pair_id}.anno
	${cmd_clean}
    """    
}


workflow ANNOTATE_PEAKS {
    take: 
    bedfile
    annotation
    
    main:
    	
		out = annotatePeaks(bedfile, annotation)
    emit:
    	out
}


workflow GET_VERSION {
    main:
		getVersion()
    emit:
    	getVersion.out
}    
