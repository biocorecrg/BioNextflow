/*
* NanoPolish
*/

params.LABEL = ""
params.EXTRAPARS = ""
params.OUTPUT = ""
params.OUTPUTMODE = "copy"
params.CONTAINER = 	"quay.io/biocontainers/f5c:1.5--h56e2c18_1"


/*
*/
process getVersion {
    container params.CONTAINER

    output:
	stdout emit: out    
    
    shell:
    """
		f5c --version | grep version      
    """
}


/*
*/

process index_slow {

    container params.CONTAINER
    label (params.LABEL)
    tag "${idsample}" 
 		
    input:
    tuple val(idsample), path(slowfile), path(fastq)
    
    output:
    tuple val(idsample), path ("${fastq}.*"), emit: fastq_idx
    tuple val(idsample), path ("${slowfile}.idx"), optional: true, emit: slow_idx

    script:
    
    """ 
      f5c index --slow5 ${slowfile} ${fastq}
    """

}

/*
*/
process eventalign_slow {
    container params.CONTAINER
    label (params.LABEL)
    tag "${idsample}" 
    if (params.OUTPUT != "") { publishDir(params.OUTPUT, mode:params.OUTPUTMODE) }
 	
    input:
    tuple val(idsample), path(slow), path(fastq), path(bam), path(bai), path(reference), path("*")
    
    output:
    tuple val(idsample), path("${idsample}_events.tsv.gz")

    script:
    """ 
    f5c eventalign ${params.EXTRAPARS} -t ${task.cpus} -b ${bam} -g ${reference} -r ${fastq} --slow5 ${slow} | gzip -c > ${idsample}_events.tsv.gz
    """
}




/*
*/

workflow GET_VERSION {
    main:
		getVersion()
    emit:
    	getVersion.out
}    
