/*
*  Epic2 
*/

params.LABEL = ""
params.EXTRAPARS = ""

params.OUTPUT = "epic2_out"
params.CONTAINER = "quay.io/biocontainers/epic2:0.0.48--py37hd0e48df_0"

include { unzipCmd } from '../global_functions.nf'

process getVersion {
    container params.CONTAINER

    output:
	stdout emit: out    
    
    shell:
    """
    echo "epic2 "`epic2 --version`
    """
}


process peakCall {

    label (params.LABEL)
    tag { comp_id }
    container params.CONTAINER
    if (params.OUTPUT != "") { publishDir(params.OUTPUT, mode:'copy') }

    input:
    tuple val(comp_id), path(sample), path(input), path(index_sample), path(index_input), val(gfrac)

    output:
    tuple val(comp_id), path("${comp_id}_epic2_peaks.bed")
    
	script:
    """
    epic2 ${params.EXTRAPARS} -t ${sample} -c ${input} -a --effective-genome-fraction ${gfrac} --output ${comp_id}_epic2_peaks.bed
    """
}

workflow CALL {
    take: 
    comparisons
    gsize
    
    main:
		out = peakCall(comparisons.combine(gsize))
    emit:
    	out
}


workflow GET_VERSION {
    main:
		getVersion()
    emit:
    	getVersion.out
}    
