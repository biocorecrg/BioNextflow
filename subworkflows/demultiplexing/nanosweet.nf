/*
*  This workflow allows to wrap the seq_tagger program for demuliplexing ONT fast5 data
*  when included you can specify the GPU param ON or OFF for using the GPU
*/

params.LABEL = ""
params.EXTRAPARS = ""
params.OUTPUT = ""
params.CONTAINER = 'biocorecrg/nanosweet:0.1'

process getVersion {
    container params.CONTAINER
    label (params.LABEL)

    output:
	stdout emit: out    
    
    shell:
    """
	echo nothing
    """
}

process demultiplex {
    tag { idfile }
    label (params.LABEL)

    container params.CONTAINER
             
    input:
    tuple val(idfile), path(fastqs), path(barcodes)

    output:
	tuple val(idfile), path("demux/*.fq.gz"), emit: demux_files
 
    script:  
    
    """	
		nanomux -b ${barcodes} ${params.EXTRAPARS} -f ${fastqs} -o ./demux  -j ${task.cpus}
    """
}




 workflow DEMULTIPLEX {
 
    take: 
    	fastq
    	barcode
    
    main:
		demultiplex(fastq.combine(barcode))
		
		out = demultiplex.out.demux_files.transpose().map {
			def new_id = "${it[0]}---${it[1].getSimpleName()}"
			[ new_id, it[1] ]
		}        

	emit:
    	out
  
}

workflow GET_VERSION {
    main:
		getVersion()
    emit:
    	getVersion.out
}   
