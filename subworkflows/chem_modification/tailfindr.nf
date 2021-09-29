/*
* Epinano
*/

params.LABEL = ""
params.EXTRAPARS = ""
params.OUTPUT = ""
params.OUTPUTMODE = "copy"
params.CONTAINER = "biocorecrg/mopnanotail:0.2"

/*
*/
process getVersion {
    container params.CONTAINER

    output:
	stdout emit: out    
    
    shell:
    """
    	echo tailfindr' '`Rscript -e "library('tailfindr'); packageVersion('tailfindr')"` 2>/dev/null
    """
}


/*
*/

process estimateTailSize {
    if (params.OUTPUT != "") { publishDir(params.OUTPUT, mode:'copy', pattern: '*.csv.gz') }

    container params.CONTAINER
    label (params.LABEL)
    tag "${sampleID}" 
 	
	input:
	tuple val(sampleID), path(fast5)

	output:
	tuple val(sampleID), path("*_findr.csv") 

	script:
	"""
	R --slave -e "library(tailfindr); find_tails(fast5_dir = './' , save_dir = './', ${params.EXTRAPARS}, csv_filename = \'${sampleID}_findr.csv\', num_cores = ${task.cpus})"
	"""
}


/*
*/

workflow ESTIMATE_TAIL {

    take: 
    fast5
    
    main:
    	out = estimateTailSize(fast5)

	emit:
		out
}

/*
*/

workflow GET_VERSION {
    main:
		getVersion()
    emit:
    	getVersion.out
}    
