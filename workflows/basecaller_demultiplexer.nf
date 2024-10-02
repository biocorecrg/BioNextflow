/*
* DEMULTIPLEXING 
*/

moduleFolder   = "../subworkflows"

params.gpu = ""
params.outmode = "copy"
params.output = ""
params.label = ""
params.label2 = ""
params.extrapars = ""
params.models = ""
params.type = "guppy"

def cuda_cont = (params.gpu == 'cuda11' ? 'biocorecrg/mopbasecallc11:0.3' : 'biocorecrg/mopbasecall:0.3')

// INCLUDE WORKFLOWS

include { separateGuppy } from "./basecaller" 


// INCLUDE MODULES / SUBWORKFLOWS

include { BASECALL_DEMULTI as GUPPY_BASECALL_DEMULTI } from "${moduleFolder}/basecalling/guppy" addParams(EXTRAPARS_DEM: params.extrapars, LABEL: params.label, GPU: params.gpu, MOP: "YES", OUTPUT: params.output, OUTPUTMODE: params.outmode, CONTAINER: cuda_cont)
include { BASECALL_DEMULTI as GUPPY6_BASECALL_DEMULTI } from "${moduleFolder}/basecalling/guppy" addParams(VERSION:"6", EXTRAPARS_DEM: params.extrapars, LABEL: params.label, GPU: params.gpu, MOP: "YES", OUTPUT: params.output, OUTPUTMODE: params.outmode, CONTAINER: cuda_cont)
include { BASECALL_DEMULTI as GUPPY65_BASECALL_DEMULTI } from "${moduleFolder}/basecalling/guppy" addParams(VERSION:"6.4", EXTRAPARS_DEM: params.extrapars, LABEL: params.label, GPU: params.gpu, MOP: "YES", OUTPUT: params.output, OUTPUTMODE: params.outmode, CONTAINER: cuda_cont)
include { BASECALL_DEMULTI as DORADO_BASECALL_DEMULTI } from "${moduleFolder}/basecalling/dorado" addParams(EXTRAPARS: params.extrapars, LABELBC: params.label, LABELCONV: params.label2,  GPU: params.gpu, OUTPUT: params.output, OUTPUTMODE: params.outmode)
include { DEMULTIPLEX as READUCKS_DEMULTIPLEX } from "${moduleFolder}/demultiplexing/readucks" addParams(EXTRAPARS: params.extrapars, LABEL: params.label, OUTPUT: params.output, OUTPUTMODE: params.outmode)


/*
* Wrapper for demultiplexing
*/
workflow BASECALL_DEMULTIPLEX {
    take: 
        fast5_4_analysis  
       
    main:
        basecalled_fastq = channel.empty()
        basecalling_stats = channel.empty()    
        basecalled_fast5 = channel.empty()    
        
        switch(params.type) {                      
           case "guppy":
           case "readucks":
               (newer, middle, older) = separateGuppy(fast5_4_analysis) 
               
               outbc65 = GUPPY65_BASECALL_DEMULTI(newer)  
               outbc6 = GUPPY6_BASECALL_DEMULTI(middle)
               outbc = GUPPY_BASECALL_DEMULTI(older)  
 
               basecalled_fastq = outbc.basecalled_fastq.concat(outbc6.basecalled_fastq).concat(outbc65.basecalled_fastq)    
               basecalled_fast5 = outbc.basecalled_fast5.concat(outbc6.basecalled_fast5).concat(outbc65.basecalled_fast5)    
               basecalling_stats = outbc.basecalling_stats.concat(outbc6.basecalling_stats).concat(outbc65.basecalling_stats)    

  
               if (params.type == "readucks") {
                    basecalled_fastq = READUCKS_DEMULTIPLEX(basecalled_fastq)
               }
    
           break;
           case "dorado":
           	   dorado_models = Channel.fromPath("${params.models}", type: "dir", checkIfExists:true)
               //dorado_models = Channel.fromPath("${params.models}/*", type: "dir", checkIfExists:true)
               out_demulti = DORADO_BASECALL_DEMULTI(fast5_4_analysis, dorado_models)
               basecalled_fastq = out_demulti.basecalled_fastq
			   basecalling_stats = out_demulti.demulti_report
		   break;
			
        }        

    emit:
        demultiplexed_fastqs =  basecalled_fastq
        basecalling_stats = basecalling_stats
        basecalled_fast5 = basecalled_fast5
}



