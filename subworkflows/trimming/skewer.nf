/*
*  skewer module
* warning is failing because of the biocontainers is not using gzip2
*/

params.LABEL = ""
params.EXTRAPARS = ""

params.OUTPUT = ""
params.CONTAINER = "biocorecrg/skewer:0.2.2"

process getVersion {
    container params.CONTAINER

    output:
	stdout emit: out    
    
    shell:
    """
	skewer --version | grep skewer | sed s/\\\t//g
	"""
}


process trimWithSkewer {
    label (params.LABEL)
    tag { pair_id }
    container params.CONTAINER
    if (params.OUTPUT != "") { publishDir(params.OUTPUT, mode:'copy') }

    input:
    tuple val(pair_id), path(reads)

    output:
    tuple val(pair_id), path("*trimmed*.fastq.gz"), emit: trimmed_reads
    path "*trimmed.log", emit: trim_log
    
    """
    skewer ${params.EXTRAPARS} -t ${task.cpus} -n -u -o ${pair_id} -z ${reads}
    """
}


workflow TRIMMING {
    take: 
    fastq
    
    main:
        out = trimWithSkewer(fastq)
    emit:
        trimmed_reads = out.trimmed_reads.map{
        	def id = it[0]
        	def reads = it.remove(1)
        	[ id, [reads] ]
        }
        trim_log = out.trim_log
}

workflow GET_VERSION {
    main:
		getVersion()
    emit:
    	getVersion.out
}    


 

