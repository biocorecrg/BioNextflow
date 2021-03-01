/*
*  lumpy module + samtools 
*/

params.LABEL = ""
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

process doIndexWGS_UCSC_NONORM {
    tag { genomefile }
    label (params.LABEL)
    container params.CONTAINER

    input:
    path(genomefile)
    path(annotation)

    output:
    path("reference") 
 
 	script:
 	def unzip_data = unzipCmd(genomefile)
 	def cmd = unzip_data[1]
 	def fname = unzip_data[0]
    """	
    ${cmd}
    cnvkit.py batch --annotate ${annotation} --target-avg-size 100000 -m wgs -n -f ${fname} -d reference
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

    output:
    path("reference") 
 
 	script:
 	def unzip_data = unzipCmd(genomefile)
 	def cmd = unzip_data[1]
 	def fname = unzip_data[0]
 	def normals = norm_samples.join(' ')
    """	
    ${cmd}
    cnvkit.py batch --annotate ${annotation} --drop-low-coverage --segment-method hmm --target-avg-size 100000 -m wgs -n ${normals} -f ${fname} -d reference
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
    tuple pair_id, path("${pair_id}_out/${pair_id}.cns"), path("${pair_id}_out/${pair_id}.cnr") 
 
 	script:
    """	
    cnvkit.py batch --drop-low-coverage --segment-method flasso -m wgs -d ${pair_id}_out -r ${reference}/reference.cnn ${reads}
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
    tuple pair_id, path("${pair_id}.diag.pdf") 
 
 	script:
    """	
	cnvkit.py diagram -s ${cnsfile} ${cnrfile} -o ${pair_id}.diag.pdf
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
    tuple pair_id, path("${pair_id}.scatter.pdf") 
 
 	script:
    """	
	cnvkit.py scatter -s ${cnsfile} ${cnrfile} -o ${pair_id}.scatter.pdf
    """    
}



workflow CNVKIT_WGS_NORM_ALL {
    take: 
    reference
    annotation
    normals_aln
    sorted_aln
    
    main:
    	normals_aln.map{
    		it[1]
    	}.collect().set{
    		normals_files
    	}
    	
		ref = doIndexWGS_UCSC(reference, annotation, normals_files)
		out = doCNV(ref, sorted_aln)
		plotDiagram(out)
		plotScatter(out)
    emit:
    	out
}

workflow CNVKIT_WGS_NONORM_ALL {
    take: 
    reference
	annotation
    sorted_aln
    
    main:
		ref = doIndexWGS_UCSC_NONORM(reference, annotation)
		out = doCNV(ref, sorted_aln)
		plotDiagram(out)
		plotScatter(out)
    emit:
    	out
}


workflow GET_VERSION {
    main:
		getVersion()
    emit:
    	getVersion.out
}    
