/*
*  This workflow allows to wrap the deeplexicon program for demuliplexing ONT fast5 data
*  when included you can specify the GPU param ON or OFF for using the GPU
*/

params.LABEL = ""
params.EXTRAPARS = ""
params.OUTPUT = ""
params.CONTAINER = 'quay.io/biocontainers/readucks:0.0.3--py_0'

process getVersion {
    container params.CONTAINER
    label (params.LABEL)

    output:
	stdout emit: out    
    
    shell:
    """
    readucks --version
    """
}

process demultiplex {
    tag { idfile }
    label (params.LABEL)

    container params.CONTAINER
             
    input:
    tuple val(idfile), path(fastq)

    output:
	tuple val(idfile), path("${idfile}.barcode*_rd.fastq.gz"), emit: demux_fastq
 
    script:  
    """	
    readucks ${params.EXTRAPARS} -p ${idfile}. -b -t ${task.cpus} -o ./ -i *.unclassified.fastq.gz

	if ls ${idfile}.barcode*.fastq.gz  1> /dev/null 2>&1; then \
		for i in ${idfile}.barcode*.fastq.gz; do cat \$i > `basename \$i .fastq.gz`_rd.fastq.gz; done \
	fi
	
	if ls ${idfile}.NB*.fastq  1> /dev/null 2>&1; then \
		for i in ${idfile}.NB*.fastq; do gzip -c \$i >> `echo \$i | sed s@NB@barcode@g | sed s@.fastq@_rd.fastq.gz@g`; done \
	fi

	
	if ls ${idfile}.unassigned.fastq 1> /dev/null 2>&1; then \
		gzip -c ${idfile}.unassigned.fastq > ${idfile}.unassigned_rd.fastq.gz; \
	fi

	rm *.fastq
    	


    """
}

 workflow DEMULTIPLEX {
    take: 
    input_fastq
    
    main:        
	demultiplex(input_fastq)
        

	emit:
    	demultiplex.out.demux_fastq
  
}

workflow GET_VERSION {
    main:
		getVersion()
    emit:
    	getVersion.out
}   
