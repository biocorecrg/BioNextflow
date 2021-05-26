/*
*  Epic 
*/

params.LABEL = ""
params.EXTRAPARS = ""

params.OUTPUT = "epic_out"
params.CONTAINER = "biocorecrg/epic:0.2.12"

include { unzipCmd } from '../global_functions.nf'
include { CALC_AVG_READSIZE as CALC_AVG_READSIZE } from "../misc/misc"


process getVersion {
    container params.CONTAINER

    output:
	stdout emit: out    
    
    shell:
    """
	epic --version
    """
}


process getEffectiveGenomeSize {
    label (params.LABEL)
    tag { reference }
    container params.CONTAINER

    input:
    val read_size
    path(reference)

    output:
	path("epic-effective.out")
    
	script:
    def unzip_ref = unzipCmd(reference)
    def cmd_ref = unzip_ref[1]
    def ref_name = unzip_ref[0]
    def cmd_clean = unzip_ref[2]
    """
    ${cmd_ref}
    epic-effective -t ./ --read-length=${read_size} -n ${task.cpus} ${ref_name} 2>/dev/null > epic-effective.out
    ${cmd_clean}
    """
}

workflow GET_EFFECTIVE_GENOME_SIZE {
    take: 
    reads
    genome
    
    main:
    	readsize = CALC_AVG_READSIZE(reads, "illumina")
		raw_out = getEffectiveGenomeSize(readsize, genome)
		raw_out.splitCsv(sep: ":", skip: 1 ).map{
			[ it[1].replaceAll("\\s","") ]
		}.toList().set{out}
    emit:
    	out
}


workflow GET_VERSION {
    main:
		getVersion()
    emit:
    	getVersion.out
}    
