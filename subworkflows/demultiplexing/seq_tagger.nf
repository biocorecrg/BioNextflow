/*
*  This workflow allows to wrap the seq_tagger program for demuliplexing ONT fast5 data
*  when included you can specify the GPU param ON or OFF for using the GPU
*/

params.LABEL = ""
params.EXTRAPARS = ""
params.OUTPUT = ""
params.CONTAINER = 'lpryszcz/seqtagger:1.0a'

process getVersion {
    container params.CONTAINER
    label (params.LABEL)

    output:
	stdout emit: out    
    
    shell:
    """
	mRNA --version
    """
}

process demultiplex {
    tag { idfile }
    label (params.LABEL)

    container params.CONTAINER
             
    input:
    tuple val(idfile), path(fast5), path("models")

    output:
	tuple val(idfile), path("${idfile}_demux.tsv.gz"), emit: demux_files
	tuple val(idfile), path("${idfile}.boxplot.pdf"), emit: demux_boxplot
 
    script:  
    """	
        mkdir tmp
        export MPLCONFIGDIR=$PWD/tmp
    	mRNA ${params.EXTRAPARS} -r -i ./ -o temp_output -t ${task.cpus}
    	mv temp_output/..demux.tsv.gz ${idfile}_demux.tsv.gz
    	mv temp_output/..demux.tsv.gz.boxplot.pdf ${idfile}.boxplot.pdf
    """
}



 workflow DEMULTIPLEX {
 
    take: 
    	input_fast5
    	model_folder
    
    main:
		demultiplex(input_fast5.combine(model_folder))
        

	emit:
    	demultiplex.out.demux_files
  
}

workflow GET_VERSION {
    main:
		getVersion()
    emit:
    	getVersion.out
}   
