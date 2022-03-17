/*
*  This workflow allows to wrap the GUPPY program for basecalling on ONT fast5 data
*  when included you can specify the GPU param ON or OFF for using the GPU
*/

params.LABEL = ""
params.EXTRAPARS = ""
params.OUTPUTMODE = "copy"
params.OUTPUT = ""
params.CONTAINER = "centos:centos6.10"


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

process baseCall {

    tag { idfile }
    label (params.LABEL)

	if (params.OUTPUT != "") { publishDir(params.OUTPUT,  mode: params.OUTPUTMODE) }
    container params.CONTAINER
             
    input:
    tuple val(idfile), path(samplesheet), path(infolder)
    
    output:
    tuple val(idfile), path("${idfile}_ouput"), emit: outfolder
    tuple val(idfile), path("${idfile}_ouput/*.fastq.gz"), emit: undet_fastqs
    tuple val(idfile), path("${idfile}_ouput/*/*/*.fastq.gz"), emit: demux_fastqs
    tuple val(idfile), path("${idfile}_ouput/Stats"), emit: statfolder

    script:
    """
        bcl2fastq ${params.EXTRAPARS} --sample-sheet ${samplesheet} \
        --output-dir ${idfile}_ouput \
        --runfolder-dir ${infolder} \
        -p ${task.cpus} -r ${task.cpus} -w ${task.cpus} 
    """
}

 workflow BASECALL {
    take: 
    input_data
    
    main:
    	baseCall(input_data)
    	
	emit:
    	outfolder = baseCall.out.outfolder
    	undet_fastqs = baseCall.out.undet_fastqs
    	demux_fastqs = baseCall.out.demux_fastqs
    	statfolder = baseCall.out.statfolder
 
}


workflow GET_VERSION {
    main:
		getVersion()
    emit:
    	getVersion.out
} 

