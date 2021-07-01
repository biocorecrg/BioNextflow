/*
*  lumpy module + samtools 
*/

params.LABEL = ""
params.EXTRAPARS_IND = ""
params.EXTRAPARS = ""

params.OUTPUT = "cnvkit_out"
params.CONTAINER = "quay.io/biocontainers/cnvkit:0.9.7--pyh9f0ad1d_0"
include { unzipCmd } from '../global_functions.nf'

process getVersion {
    container params.CONTAINER

    output:
	stdout emit: out    
    
    shell:
    """
    cnvkit.py version| awk '{print "cnvkit\t"\$0}'
    """
}

process estimateBinSize {
    tag { "${genomefile}" }
    
    label (params.LABEL)
    container params.CONTAINER

    input:
    path(genomefile)
    path(bamfiles)
    path(indexes)

    output:
    stdout
 
 	script:
 	def unzip_data = unzipCmd(genomefile)
 	def cmd = unzip_data[1]
 	def fname = unzip_data[0]
    """
    ${cmd}
	cnvkit.py access ${fname} -s 10000 -o access-10kb.bed
	cnvkit.py autobin *.bam -m wgs -b 50000 -g access-10kb.bed > stats.txt
	cat stats.txt | tail -n 1|cut -f 3 
    """    
}

//

process doIndexWGS_UCSC_NONORM {
    tag { genomefile }
    label (params.LABEL)
    container params.CONTAINER

    input:
    path(genomefile)
    path(annotation)
    val(breaksize)

    output:
    path("reference") 
 
 	script:
 	def unzip_data = unzipCmd(genomefile)
 	def cmd = unzip_data[1]
 	def fname = unzip_data[0]
    """	
    ${cmd}
    cnvkit.py batch ${params.EXTRAPARS_IND} --target-avg-size ${breaksize} --annotate ${annotation} -m wgs -n -f ${fname} -d reference
    grep "chr" reference/reference.cnn > reference/reference2.cnn
    mv reference/reference.cnn reference/reference.cnn.ori
    mv reference/reference2.cnn reference/reference.cnn
    """    
}

process doIndexWGS_UCSC {
    tag { genomefile }
    label (params.LABEL)
    container params.CONTAINER

    input:
    path(genomefile)
    path(annotation)
    path(norm_samples)
    val(breaksize)
  

    output:
    path("reference") 
 
 	script:
 	def unzip_data = unzipCmd(genomefile)
 	def cmd = unzip_data[1]
 	def fname = unzip_data[0]
 	def normals = norm_samples.join(' ')
    """	
    ${cmd}
    cnvkit.py batch ${params.EXTRAPARS_IND} --target-avg-size ${breaksize} --annotate ${annotation} --segment-method hmm -m wgs -n ${normals} -f ${fname} -d reference
    grep "chr" reference/reference.cnn > reference/reference2.cnn
    mv reference/reference.cnn reference/reference.cnn.ori
    mv reference/reference2.cnn reference/reference.cnn
    """    
}


process doCNV {
    tag { pair_id }
    label (params.LABEL)
    container params.CONTAINER
    if (params.OUTPUT != "") { publishDir(params.OUTPUT, mode:'copy') }

    input:
    path(reference)
    tuple val(pair_id), path(reads)

    output:
    tuple val(pair_id), path("${pair_id}_out/${pair_id}.cns"), path("${pair_id}_out/${pair_id}.cnr") 
 
 	script:
    """	
    cnvkit.py batch ${params.EXTRAPARS} --segment-method hmm -m wgs -d ${pair_id}_out -r ${reference}/reference.cnn ${reads}
    mv ${pair_id}_out/`basename ${reads} .bam`.cns ${pair_id}_out/${pair_id}.cns 
    mv ${pair_id}_out/`basename ${reads} .bam`.cnr ${pair_id}_out/${pair_id}.cnr 
    """    
}

process plotDiagram {
    tag { pair_id }
    label (params.LABEL)
    container params.CONTAINER
    publishDir(params.OUTPUT, mode:'copy')

    input:
    tuple val(pair_id), path(cnsfile), path(cnrfile)


    output:
    tuple val(pair_id), path("${pair_id}.diag.pdf") 
 
 	script:
    """	
	cnvkit.py diagram  -s ${cnsfile} ${cnrfile} -o ${pair_id}.diag.pdf
    """    
}

process plotScatter {
    tag { pair_id }
    label (params.LABEL)
    container params.CONTAINER
    publishDir(params.OUTPUT, mode:'copy')

    input:
    tuple val(pair_id), path(cnsfile), path(cnrfile)


    output:
    tuple val(pair_id), path("${pair_id}.scatter.pdf") 
 
 	script:
    """	
	cnvkit.py scatter  -s ${cnsfile} ${cnrfile} -o ${pair_id}.scatter.pdf
    """    
}

process calcBreaks {
    tag { pair_id }
    label (params.LABEL)
    container params.CONTAINER
    publishDir(params.OUTPUT, mode:'copy')

    input:
    tuple val(pair_id), path(cnsfile), path(cnrfile)


    output:
    tuple val(pair_id), path("${pair_id}.breaks.txt") 
 
 	script:
    """	
    cnvkit.py breaks ${cnrfile} ${cnsfile} > ${pair_id}.breaks.txt
    """    
}



workflow WGS_NORM_ALL {
    take: 
    reference
    annotation
    normals_files
    sorted_aln
    bamfiles
    indexes_files
    
    
    main:
    	binsize_ch = estimateBinSize(reference, bamfiles, indexes_files)
    	binsize_ch.map{it.trim().toInteger()}.set{binsize}
		ref = doIndexWGS_UCSC(reference, annotation, normals_files, binsize)
		out = doCNV(ref, sorted_aln)
		plotDiagram(out)
		plotScatter(out)
		calcBreaks(out)
    emit:
    	out
}

workflow WGS_NONORM_ALL {
    take: 
    reference
	annotation
    sorted_aln
    bamfiles
    indexes
    
    main:
     	breaksize = estimateBinSize(reference, bamfiles, indexes)
		ref = doIndexWGS_UCSC_NONORM(reference, annotation, breaksize)
		out = doCNV(ref, sorted_aln)
		plotDiagram(out)
		plotScatter(out)
		calcBreaks(out)
    emit:
    	out
}


workflow GET_VERSION {
    main:
		getVersion()
    emit:
    	getVersion.out
}    
