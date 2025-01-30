/*
* NanoPolish
*/

params.LABEL = ""
params.LABELPREP = ""
params.LABELINF = ""
params.EXTRAPARS = ""
params.OUTPUT = ""
params.OUTPUTMODE = "copy"
params.CONTAINER = 	"biocorecrg/remora:3.0.0"

include { unzipCmd } from '../global_functions.nf'

/*
*/
process getVersion {
    container params.CONTAINER

    output:
	stdout emit: out    
    
    shell:
    """
		remora -v       
    """
}

process dataprep {
//    if (params.OUTPUT != "") { publishDir(params.OUTPUT, mode:'copy') }

    container params.CONTAINER
    label (params.LABEL)
    label (params.LABELPREP)
    tag "${id}" 
 		
    input:
    tuple val(id), path(pod5), path(bam)
    
    //output:
    //tuple val(id), path("${id}")
    
    script:

	"""
	remora \
  		dataset prepare \
		${pod5} \
  		${bam} \
  		--output-path can_chunks \
  		--refine-kmer-level-table levels.txt \
  		--refine-rough-rescale \
  		--motif CG 0 \
  		--mod-base-control
    """



}




/*
*/

process inference {
    if (params.OUTPUT != "") { publishDir(params.OUTPUT, mode:'copy') }

    container params.CONTAINER
    label (params.LABEL)
    label (params.LABELINF)
    
    tag "${sampleID}" 
 		
    input:
    tuple val(sampleID), path(input)
    
    output:
    tuple val(sampleID), path("${sampleID}_inf")
    
    script:
	"""
    m6anet inference --input_dir ./${input} --out_dir ./${sampleID}_inf ${params.EXTRAPARS}
    """



}

/*
*/

workflow INFERENCE_RNA004 {

    take: 
    input
    
    main:
    m6data = dataprep(input)
    
	out = inference(m6data)

	emit:
	m6data
 }


/*
*/

workflow GET_VERSION {
    main:
		getVersion()
    emit:
    	getVersion.out
}    
