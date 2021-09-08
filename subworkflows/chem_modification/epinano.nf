/*
* Epinano
*/

params.LABEL = ""
params.EXTRAPARS = ""
params.OUTPUT = ""
params.OUTPUTMODE = "copy"
params.CONTAINER = "biocorecrg/mopepinano:0.1"

/*
*/
process getVersion {
    container params.CONTAINER

    output:
	stdout emit: out    
    
    shell:
    """
		echo "Epinano v1.1.1"      
    """
}


/*
*/

process calcVarFrequencies {
    if (params.OUTPUT != "") { publishDir(params.OUTPUT, mode:'copy', pattern: '*.csv.gz') }

    container params.CONTAINER
    label (params.LABEL)
    tag "${sampleID}" 
 	
    input:
    tuple val(sampleID), path(tsvfile)
    
    output:
    tuple val(sampleID), path("*.tsv.per.site.var.csv.gz"), emit: per_site_vars
    tuple val(sampleID), path("*.csv.gz"), emit: all_vars
    
    script:
	"""
	TSV_to_Variants_Freq.py3 -f ${tsvfile} -t ${task.cpus}
	for i in *.csv; do gzip \$i; done
	"""
}

/*
*/

workflow CALC_VAR_FREQUENCIES {

    take: 
    tsvfiles
    
    main:
    	calcVarFrequencies(tsvfiles)

	emit:
    	per_site_vars = calcVarFrequencies.out.per_site_vars
    	all_vars = calcVarFrequencies.out.all_vars
}

/*
*/

workflow GET_VERSION {
    main:
		getVersion()
    emit:
    	getVersion.out
}    
