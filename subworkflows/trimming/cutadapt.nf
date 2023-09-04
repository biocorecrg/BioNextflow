/*
* cutadapt module
*/

include { separateSEandPE } from '../global_functions.nf'

params.CONTAINER = "quay.io/biocontainers/cutadapt:4.4--py39hf95cd2a_1"
params.OUTPUT = ""

params.LABEL = ""
params.EXTRAPARS = ""


process getVersion {
    container params.CONTAINER

    output:
    stdout emit: out

    shell:
    """
    cutadapt --version | awk '{print "cutadapt "\$0}'
    """
}

process filterFastaOpts {
    if (params.OUTPUT != "") {publishDir(params.OUTPUT, mode: 'copy') }
    
    tag {id }
    container params.CONTAINER
    label (params.LABEL)

    input:
    tuple val (id), path(seqs), val(extraopts)

    output:
    tuple val (id), path ("${id}_trimmed.fasta"), emit: trimmed_reads

    script:

    """
    cutadapt \
    -j ${task.cpus} \
    -o ${id}_trimmed.fasta ${extraopts} \
    ${params.EXTRAPARS} ${seqs} 

    """
}



workflow FILTER_FASTA_OPTS{
    take: 
    seqs
    
    main:
        
    out = filterFastaOpts(seqs)

    emit:
 		out
}


workflow GET_VERSION {
    main:
		getVersion()
    emit:
    	getVersion.out
} 


