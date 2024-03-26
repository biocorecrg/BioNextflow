/*
* Tailfindr
*/

params.LABEL = ""
params.EXTRAPARS = ""
params.OUTPUT = ""
params.OUTPUTMODE = "copy"
params.CONTAINER = "biocorecrg/moptail:1.3"
params.MODE = "default"

def mycontainer = params.CONTAINER
if (params.MODE == 'n3ps_r9') {
	mycontainer = 'biocorecrg/moptail:nano3p_5'
} else if (params.MODE == 'n3ps_r10') {
        mycontainer = 'biocorecrg/moptail_nano3p_5_r10:0.2'
}


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

    container mycontainer
    label (params.LABEL)
    tag "${sampleID}" 
 	
	input:
	tuple val(sampleID), path(fast5)

	output:
	tuple val(sampleID), path("*_findr.csv.gz") 

	script:
	"""
	R --vanilla --slave -e "library(tailfindr); find_tails(fast5_dir = './' , save_dir = './', ${params.EXTRAPARS}, csv_filename = \'${sampleID}_findr.csv\', num_cores = ${task.cpus})"
	gzip *_findr.csv
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
