/*
*  Genrich 
*/

params.LABEL = ""
params.EXTRAPARS = ""

params.OUTPUT = ""
params.CONTAINER = "quay.io/biocontainers/genrich:0.6.1--he4a0461_4"

include { unzipCmd } from '../global_functions.nf'

process getVersion {
    container params.CONTAINER

    output:
	stdout emit: out    
    
    shell:
    """
    Genrich --version 2>&1 | grep "version"
    """
}


process peakCallAtac {
    label (params.LABEL)
    tag { comp_id }
    container params.CONTAINER
    if (params.OUTPUT != "") { publishDir(params.OUTPUT, mode:'copy') }

    input:
    tuple val(comp_id), path(sample)

    output:
    tuple val(comp_id), path("${comp_id}_peaks.narrowPeak"), optional: true, emit: narrowPeaks
    
	script:
    """
	Genrich ${params.EXTRAPARS} -r -j -t ${sample}  -o ${comp_id}_peaks.narrowPeak
    """
}




workflow ATAC_CALL {
    take: 
    samples
    
    main:
		peakCallAtac(samples)
    emit:
    	narrowPeaks = peakCallAtac.out.narrowPeaks
}


workflow GET_VERSION {
    main:
		getVersion()
    emit:
    	getVersion.out
}    
