params.OUTPUT = ""
params.OUTPUTMODE = "copy"
params.TYPE = "seqtagger"
params.CONTAINER = "biocorecrg/demxufastq:0.2"

process extract_deeplexicon_fastq {

    tag "${ idfile }"

    container params.CONTAINER
			
	input:
	tuple val(idfile), path(demux), path(fastq) 
	
	
	output:
	tuple val(idfile), path ("*.fastq.gz")

	script:
	"""
		extract_sequence_from_fastq.py ${demux} ${fastq}
		for i in *.fastq; do gzip \$i; done
	"""
}

process extract_seqtagger_fastq {

    tag { idfile }

    container params.CONTAINER
             
 	input:
	tuple val(idfile), path(demux), path(fastq) 
	
	
	output:
	tuple val(idfile), path ("*.fastq.gz")

	script:
	"""
		fastq_split_by_barcode.py -b 50 -i ${demux} -f ${fastq} -o ${idfile}
	"""
	
}

workflow DEMULTI_FASTQ {

    take: 
    basecalled_tsv
    basecalled_fastq
    
    main:
       out = channel.empty()
       switch(params.TYPE) {                      
           case "seqtagger":
                out = extract_seqtagger_fastq(basecalled_tsv.join(basecalled_fastq))
           break;
           case "deeplexicon":
				out = extract_deeplexicon_fastq(basecalled_tsv.join(basecalled_fastq))
    	}
    	
    emit:
    	out
 
 }
