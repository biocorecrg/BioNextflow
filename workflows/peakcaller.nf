/*
* ALIGNER 
*/

moduleFolder   = "../subworkflows"

params.outmode = "copy"
params.output = ""
params.label = ""
params.extrapars = ""
params.type = "macs2"

// include subworkflows
include { GET_VERSION as EPIC2_VER;  CALL as EPIC2_CALL_CHIP} from "${moduleFolder}/peak_calling/epic2" addParams(OUTPUT: params.output, LABEL: params.label, EXTRAPARS: params.progPars["epic2"])  
include { GET_VERSION as SEACR_VER; CALL as SEACR_CALL } from "${moduleFolder}/peak_calling/seacr" addParams(OUTPUT: params.output, LABEL: params.label, EXTRAPARS: params.progPars["seacr"]) 
include { MAKE_GRAPHS as SEACR_MAKE_GRAPHS; MAKE_NORMALIZED_GRAPHS as SEACR_MAKE_NORMALIZED_GRAPHS } from "${moduleFolder}/peak_calling/seacr" 
include { CALL as MACS2_CALL_CHIP } from "${moduleFolder}/peak_calling/macs2" addParams(OUTPUT: params.output, LABEL: params.label, EXTRAPARS: params.progPars["macs2"]) 
include { CALL as MACS3_CALL_CHIP } from "${moduleFolder}/peak_calling/macs3" addParams(OUTPUT: params.output, LABEL: params.label, EXTRAPARS: params.progPars["macs3"]) 
include { GET_VERSION as BB_VER; BEDCLIP } from "${moduleFolder}/misc/bedclip"  




workflow PEAKCALL {	

    take: 
    peakcall_params
    dedup_indexed_reads
    peak_call_data
    
    main:

 	// PREPARE PEAK CALL DATA  
  	dedup_indexed_reads.combine(dedup_indexed_reads).map{
  		["${it[0]}-${it[3]}", [it[1], it[4]]]
  	}.set{comb_reads}
  	
 	peakcall_params.join(comb_reads).map{
 		[it[0], it[4][0], it[4][1]]
 	}.set{peak_data}
 	
 	// RAISE ERROR IF NO OVERLAP BETWEEN ID COMBINATIONS AND SAMPLE NAMES
 	peak_data.ifEmpty{exit 1, "PEAK DATA IS EMPTY. PLEASE CHECK THE ID COMBINATIONS YOU INDICATED ARE COINCIDING WITH SAMPLE NAMES\n\n!!!"}
 	 
   	// IF NO INDEXES, DO NOT FAIL 
  	dedup_indexed_reads.combine(dedup_indexed_reads).map{
  		if (it[5]) {
	  		["${it[0]}-${it[3]}", [it[2], it[5]]]
	  	} else {
	  		[null]
	  	}
  	}.set{comb_indexes} 	
 
 	peakcall_params.join(comb_indexes).map{
 	 		[it[0], it[4][0], it[4][1]]
 	}.set{index_data}	

	peak_res = channel.empty()

	switch( params.type ) {
      case "macs2":
		 macs_res = MACS2_CALL_CHIP(peak_data, peak_call_data["mappable_gsize"])
		 macs_peaks = macs_res.broadPeaks.mix(macs_res.narrowPeaks)
		 peak_res = convertTo6Bed(macs_peaks)
      break;   
      case "macs3":
		 macs_res = MACS3_CALL_CHIP(peak_data, peak_call_data["mappable_gsize"])
		 macs_peaks = macs_res.broadPeaks.mix(macs_res.narrowPeaks)
		 peak_res = convertTo6Bed(macs_peaks)
      break;   
      case "epic2":
         epic2_res = EPIC2_CALL_CHIP(peak_data.join(index_data), peak_call_data["gfrac"])
         epic_clip = BEDCLIP(peak_call_data["genome_fai"], epic2_res)
         epic_6_peaks = convertEpicTo6Bed(epic_clip)      
      break;   
      case "seacr":
        peak_res = channel.empty()
      break;   
    }
    
	emit:
		peaks = peak_res
		combs = peak_data
	
}


process convertTo6Bed {
    tag { id }

    input:
    tuple val(id), path(peak)

    output:
    tuple val(id), path("${id}_peaks.bed")
    
	script:
    """
    awk -F"\t" '{OFS="\t"; print \$1,\$2,\$3,\$4,\$5,\$6}' ${peak} | awk '{OFS="\t"; score=int(\$5); \
    if (score>1000) {score = 1000}; print \$1,\$2,\$3,\$4, score, \$6}' > ${id}_peaks.bed
    """
}

process convertEpicTo6Bed {
    tag { id }

    input:
    tuple val(id), path(peak)

    output:
    tuple val(id), path("${id}_peaks.bed")
    
	script:
    """
    grep -v "#" ${peak} | awk '{OFS="\t"; num++; score=int(\$5); \
    if (score>1000) {score = 1000}; print \$1,\$2,\$3,"peak_"num, score, \$6}' > ${id}_peaks.bed
    """
}