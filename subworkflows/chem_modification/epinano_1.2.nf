/*
* Epinano
*/

params.LABEL = ""
params.EXTRAPARS = ""
params.OUTPUT = ""
params.OUTPUTMODE = "copy"
params.CONTAINER = "biocorecrg/mopepinano:0.2"

/*
*/
process getVersion {
    container params.CONTAINER

    output:
	stdout emit: out    
    
    shell:
    """
		echo "Epinano v1.2"      
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
    tuple val(sampleID), path(alnfile), path(bais), path(reference), path(dict_index), path(faiidx) 

    
    output:
    tuple val(sampleID), path("*.per.site.csv.gz"), emit: per_site_vars
    
    script:
	"""
	Epinano_Variants.py -n ${task.cpus} -R ${reference} -b ${alnfile} -s \$SAM2TSV --type t 
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
