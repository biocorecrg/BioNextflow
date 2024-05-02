/*
* Epinano
*/

params.LABEL = ""
params.EXTRAPARS = ""
params.OUTPUT = ""
params.OUTPUTMODE = "copy"
params.CONTAINER = "biocorecrg/mopepinano1.2.4"

/*
*/
process getVersion {
    container params.CONTAINER

    output:
	stdout emit: out    
    
    shell:
    """
		echo "Epinano v1.2.4"      
    """
}


/*
*/

process calcVarFrequencies {
    if (params.OUTPUT != "") { publishDir(params.OUTPUT, mode:'copy', pattern: '*.csv.gz') }

    container params.CONTAINER
    label (params.LABEL)
    tag "${alnfile} on ${sampleID}" 
 	
    input:
    tuple val(sampleID), path(alnfile), path(bais), path(reference), path(dict_index), path(faiidx) 

    
    output:
    tuple val(sampleID), path("*.per.site.csv.gz"), emit: per_site_vars
    
    script:
	"""
	Epinano_Variants.py -c ${task.cpus} -r ${reference} -b ${alnfile}
	for i in *.csv; do gzip \$i; done
	"""
}

/*
*/

workflow CALC_VAR_FREQUENCIES {

    take: 
    bams
    indexes

    main:
    	calcVarFrequencies(bams.join(bais).combine(indexes))

	emit:
    	per_site_vars = calcVarFrequencies.out.per_site_vars
//    	all_vars = calcVarFrequencies.out.all_vars
}

/*
*/

workflow GET_VERSION {
    main:
		getVersion()
    emit:
    	getVersion.out
}    
