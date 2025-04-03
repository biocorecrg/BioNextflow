/*
* ALIGNER 
*/

moduleFolder   = "../subworkflows"

params.outmode = "copy"
params.output = ""
params.label = ""
params.extrapars = ""
params.extraparsind = ""
params.type = "bwa"

// include subworkflows
include { GET_VERSION as BWA_VER; ALL as BWA_ALL } from "${moduleFolder}/alignment/bwa" addParams(EXTRAPARS: params.extrapars, OUTPUT: params.output , LABEL: params.label) 
include { GET_VERSION as STAR_VER; ALL as STAR_ALL } from "${moduleFolder}/alignment/star" addParams(EXTRAPARS: params.extrapars, OUTPUT: params.output , LABEL: params.label) 
include { GET_VERSION as BOWTIE2_VER; ALL as BOWTIE2_ALL } from "${moduleFolder}/alignment/bowtie2" addParams(EXTRAPARS: params.extrapars, OUTPUT: params.output, LABEL: params.label) 
include { GET_VERSION as BISCUIT_VER; ALL as BISCUIT_ALL } from "${moduleFolder}/alignment/biscuit" addParams(EXTRAPARS: params.extrapars, OUTPUT: params.output, LABEL: params.label) 
include { GET_VERSION as GRAPHMAP2_VER; MAP as GRAPHMAP2} from "${moduleFolder}/alignment/graphmap2" addParams(EXTRAPARS: params.extrapars, OUTPUT: params.output, LABEL: params.label)
include { GET_VERSION as GRAPHMAP_VER; MAP as GRAPHMAP} from "${moduleFolder}/alignment/graphmap" addParams(EXTRAPARS: params.extrapars, OUTPUT: params.output, LABEL: params.label)
include { GET_VERSION as MINIMAP2_VER; MAP as MINIMAP2} from "${moduleFolder}/alignment/minimap2" addParams(EXTRAPARS: params.extrapars, OUTPUT: params.output, LABEL: params.label)
include { GET_VERSION as WINNOWMAP_VER; ALL as WINNOWMAP} from "${moduleFolder}/alignment/winnowmap" addParams(EXTRAPARSINDEX: params.extraparsind, EXTRAPARS: params.extrapars, OUTPUT: params.output, LABEL: params.label)

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
      case "biscuit":
        map_res = BISCUIT_ALL(genome, reads)
      break;   
      case "graphmap2":
        map_res = GRAPHMAP2(reads, genome).bams
      break;
      case "graphmap":
        map_res = GRAPHMAP(reads, genome).bams
      break;
      case "minimap2":
        map_res = MINIMAP2(reads, genome).bams
      break;     
      case "winnowmap":      	
        map_res = WINNOWMAP(reads, genome).bams
      break;
      default:
      	map_res = channel.empty() 
      	map_res.ifEmpty{exit 1, "UNSUPPORTED ALIGNER \nENDING NOW, BYE!!!"}
      break;   
       
    }
	emit:
		out = map_res
	
}


