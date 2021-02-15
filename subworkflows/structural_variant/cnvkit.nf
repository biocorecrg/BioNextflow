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

process doIndexWGS {
    label (params.LABEL)
    container params.CONTAINER

    input:
    path(genomefile)

    output:
    tuple val(pair_id), path("${pair_id}.discordants.bam") 
 
 	script:
 	
    """
    cnvkit.py batch -m wgs -n -f ${genomefile} -d reference 
    """    
}



workflow CNVKIT_WGS_ALL {
    take: 
    reference
    
    main:
		doIndexWGS(reference)
    //emit:
    //	out
}


workflow GET_VERSION {
    main:
		getVersion()
    emit:
    	getVersion.out
}    
