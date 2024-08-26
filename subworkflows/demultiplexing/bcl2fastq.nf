/*
*  This workflow allows to wrap the GUPPY program for basecalling on ONT fast5 data
*  when included you can specify the GPU param ON or OFF for using the GPU
*/

params.LABEL = ""
params.EXTRAPARS = ""
params.OUTPUTMODE = "copy"
params.OUTPUT = ""
params.CONTAINER = "centos:centos6.10"
params.PARSING = "false"


process getVersion {
    container params.CONTAINER
    label (params.LABEL)

    output:
	stdout emit: out    
    
    shell:
    """
		bcl2fastq --version 2>&1 | grep bcl2fastq  
    """
}

process demultiplex {

    tag {"${idfile} with ${samplesheet.getName()}" }
    label (params.LABEL)

    if (params.OUTPUT != "") { publishDir(params.OUTPUT,  mode: params.OUTPUTMODE) }
    container params.CONTAINER
             
    input:
    tuple val(idfile), path(samplesheet), path(infolder)
    
    output:
    tuple val(idfile), path("${idfile}_output"), emit: outfolder
    tuple val(idfile), path("${idfile}_output/*.fastq.gz"), emit: undet_fastqs
    tuple val(idfile), path("${idfile}_output/*/*/*.fastq.gz"), emit: demux_fastqs
    tuple val(idfile), path("${idfile}_output/Stats"), emit: statfolder
    tuple val(idfile), path("IndexMetricsOut.bin", optional: true), emit: indexmetrics
    tuple val(idfile), path("samplesheet*.csv", optional: true, includeInputs: true), emit: samplesheet

    script:

	def parsing_cmd = "par_cmd=''"
	
	if (params.PARSING == "true") {
		parsing_cmd = "par_cmd=`grep OverrideCycles ${samplesheet} | awk -F ',' '{print \"--use-bases-mask \" \$2 }' | sed s@';'@,@g`"
	} 

    """
    	${parsing_cmd}
        bcl2fastq ${params.EXTRAPARS} \$par_cmd --sample-sheet ${samplesheet} \
        --output-dir ${idfile}_output \
        --runfolder-dir ${infolder} \
        -p ${task.cpus} -r ${task.cpus} -w ${task.cpus} 
    """
}


 workflow CONV_DEMUX {
    take: 
    input_data
    
    
    main:
    	demultiplex(input_data)
    	
	emit:
    	outfolder = demultiplex.out.outfolder
    	undet_fastqs = demultiplex.out.undet_fastqs
    	demux_fastqs = demultiplex.out.demux_fastqs
    	statfolder = demultiplex.out.statfolder
    	indexmetrics = demultiplex.out.indexmetrics
 
}


workflow GET_VERSION {
    main:
		getVersion()
    emit:
    	getVersion.out
} 

