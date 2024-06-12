/*
* ALIGNER 
*/

moduleFolder   = "../subworkflows"

params.outmode = "copy"
params.output = ""
params.label = ""
params.extrapars = ""
params.type = "bwa"

// include subworkflows
include { GET_VERSION as BWA_VER; ALL as BWA_ALL } from "${moduleFolder}/alignment/bwa" addParams(EXTRAPARS: params.progPars["bwa"], OUTPUT: params.output , LABEL: params.label) 
include { GET_VERSION as STAR_VER; ALL as STAR_ALL } from "${moduleFolder}/alignment/star" addParams(EXTRAPARS: params.extrapars, OUTPUT: params.output , LABEL: params.label) 
include { GET_VERSION as BOWTIE2_VER; ALL as BOWTIE2_ALL } from "${moduleFolder}/alignment/bowtie2" addParams(EXTRAPARS: params.extrapars, OUTPUT: params.output, LABEL: params.label) 
include { GET_VERSION as BISCUIT_VER; ALL as BISCUIT_ALL } from "${moduleFolder}/alignment/biscuit" addParams(EXTRAPARS: params.extrapars, OUTPUT: params.output, LABEL: params.label) 

workflow ALIGN {	

    take: 
	reads
	genome

    main:
	  switch( params.type ) {
      case "bwa":
		map_res = BWA_ALL(genome, reads)
      break;   
      case "bowtie2":
	    map_res = BOWTIE2_ALL(genome, reads).aln
      break;   
      case "star":
        map_res = channel.empty()
      break;   
      case "biscuit":
        map_res = BISCUIT_ALL(genome, reads)
      break;   
    }
	emit:
		map_res
	
}
