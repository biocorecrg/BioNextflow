/*
*  This workflow allows to wrap the deeplexicon program for demuliplexing ONT fast5 data
*  when included you can specify the GPU param ON or OFF for using the GPU
*/

params.LABEL = ""
params.EXTRAPARS = ""
params.OUTPUT = ""
params.GPU = ""
params.CONTAINER = (params.GPU == "ON" ? 'lpryszcz/deeplexicon:1.2.0-gpu': 'lpryszcz/deeplexicon:1.2.0')

process getVersion {
    container params.CONTAINER
    label (params.LABEL)

    output:
	stdout emit: out    
    
    shell:
    """
    deeplexicon_sub.py --version | grep Deep
    """
}

process demultiplex {
    tag { idfile }
    label (params.LABEL)

    container params.CONTAINER
             
    input:
    tuple val(idfile), path(fast5), path(models)

    output:
	tuple val(idfile), path("${idfile}_demux.tsv"), emit: demux_files
 
    script:  
    """	
    	deeplexicon_sub.py dmux ${params.EXTRAPARS} -p ./ > ${idfile}_demux.tsv
    """
}

 workflow DEMULTIPLEX {
    take: 
    models
    input_fast5
    
    main:
    	demultiplex(input_fast5.combine(models))

	emit:
    	demultiplex.out.demux_files
  
}

workflow GET_VERSION {
    main:
		getVersion()
    emit:
    	getVersion.out
}   
