/*
* DEMULTIPLEXING 
*/

moduleFolder   = "../subworkflows"

params.gpu = ""
params.outmode = "copy"
params.output = ""
params.label = ""
params.extrapars = ""
params.container = ""
params.models = ""

// INCLUDE MODULES / SUBWORKFLOWS
def type = params.type
if (params.type == "seqtagger-trna"){
	type = "seqtagger"
}

include { DEMULTIPLEX as SEQTAGGER_DEMULTIPLEX_TRNA } from "${moduleFolder}/demultiplexing/seq_tagger" addParams(CONTAINER: params.container, EXTRAPARS: params.extrapars, LABEL: params.label)
include { DEMULTIPLEX as SEQTAGGER_DEMULTIPLEX } from "${moduleFolder}/demultiplexing/seq_tagger" addParams(EXTRAPARS: params.extrapars, LABEL: params.label)
include { DEMULTIPLEX as DEMULTIPLEX_DEEPLEXICON } from "${moduleFolder}/demultiplexing/deeplexicon" addParams(EXTRAPARS: params.extrapars, LABEL: params.label, GPU: params.gpu)
include { DEMULTI_FASTQ } from "${moduleFolder}/misc/demulti_fastq" addParams(TYPE: type)

/*
* Wrapper for demultiplexing
*/
workflow DEMULTIPLEX {
    take: 
        fast5_4_analysis
        basecalled_fastq
       
    main:
        demufq = channel.empty()
        demux = channel.empty()
        switch(params.type) {                      
           case "deeplexicon":
			   demux_models = Channel.fromPath( "${params.models}/*.h5")
               demux = DEMULTIPLEX_DEEPLEXICON(demux_models, fast5_4_analysis)
               break;
           case "seqtagger": 
			   demux_models = Channel.fromPath( "${params.models}/*", type: 'dir').collect()
               demux = SEQTAGGER_DEMULTIPLEX(fast5_4_analysis, demux_models)
               break;
           case "seqtagger-trna": 
			   demux_models = Channel.fromPath( "${params.models}/*", type: 'dir').collect()			   
               demux = SEQTAGGER_DEMULTIPLEX_TRNA(fast5_4_analysis, demux_models)
               break;
        }        
        		
        demufq = DEMULTI_FASTQ(demux, basecalled_fastq)

    emit:
          demultiplexed_fastq =  demufq
          demultiplexed_tsv =  demux
}


