/*
*  This workflow allows to wrap the deeplexicon program for demuliplexing ONT fast5 data
*  when included you can specify the GPU param ON or OFF for using the GPU
*/

params.LABEL = ""
params.EXTRAPARS = ""
params.OUTPUT = ""
params.GPU = ""
//params.CONTAINER = 'lpryszcz/deeplexicon:1.2.0'
params.CONTAINER = (params.GPU == "ON" ? 'lpryszcz/deeplexicon:1.2.0-gpu': 'lpryszcz/deeplexicon:1.2.0')

process demultiplex {
    tag { idfile }
    label (params.LABEL)

    container params.CONTAINER
             
    input:
    tuple val(idfile), path(fast5)

    output:
	tuple val(idfile), path("${idfile}_demux.tsv"), emit: demux_files
 
    script:  
    """	
    	deeplexicon_sub.py dmux ${params.EXTRAPARS} -p ./ > ${idfile}_demux.tsv
    """
}

 workflow DEMULTIPLEX {
    take: 
    input_fast5
    
    main:
    	demultiplex(input_fast5)

	emit:
    	demultiplex.out.demux_files
  
}

 
