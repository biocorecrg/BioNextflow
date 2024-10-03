/*
* BASECALLING 
*/

moduleFolder   = "../subworkflows"

params.gpu = ""
params.outmode = "copy"
params.output = ""
params.label = ""
params.label2 = ""
params.extrapars = ""
params.models = ""
params.type = "dorado"

def cuda_cont = (params.gpu == 'cuda11' ? 'biocorecrg/mopbasecallc11:0.3' : 'biocorecrg/mopbasecall:0.3')

include { GET_VERSION as GUPPY_VERSION } from "${moduleFolder}/basecalling/guppy" 
include { print_log_message } from "${moduleFolder}/global_functions.nf"
include { BASECALL as GUPPY_BASECALL } from "${moduleFolder}/basecalling/guppy" addParams(LABEL: params.label, EXTRAPARS_BC: params.extrapars, GPU: params.gpu, MOP: "YES", OUTPUT: params.output, CONTAINER: cuda_cont, OUTPUTMODE: params.outmode )
include { BASECALL as GUPPY6_BASECALL } from "${moduleFolder}/basecalling/guppy" addParams(VERSION:"6", EXTRAPARS_BC: params.extrapars, LABEL: params.label, GPU: params.gpu, MOP: "YES", OUTPUT: params.output, CONTAINER: cuda_cont, OUTPUTMODE: params.outmode)
include { BASECALL as GUPPY64_BASECALL } from "${moduleFolder}/basecalling/guppy" addParams(VERSION:"6.4", EXTRAPARS_BC: params.extrapars, LABEL: params.label, GPU: params.gpu, MOP: "YES", OUTPUT: params.output, CONTAINER: cuda_cont, OUTPUTMODE: params.outmode)
include { BASECALLMOD as DORADO_BASECALL } from "${moduleFolder}/basecalling/dorado" addParams(EXTRAPARS: params.extrapars, LABELBC: params.label, LABELCONV:params.label2, GPU: params.gpu, MOP: "YES", OUTPUT: params.output, OUTPUTMODE: params.outmode)

// ADD A CHECK FOR GUPPY FOR DISABLING SCORE

def separateGuppy (fast5) {

	ver = GUPPY_VERSION()
    ver.view { gver -> print_log_message("The GUPPY VERSION IS: ${gver.trim()}") }
	
	data_and_ver = ver.map{
        def vals = it.split("\\.")
    	"${vals[0]}.${vals[1]}"
    }.toBigDecimal().combine(fast5)	
	
    older = data_and_ver.map{ if (it[0] < 6 ) [ it[1], it[2] ]}
    middle = data_and_ver.map{ if (it[0] >= 6 && it[0] <= 6.3) [ it[1], it[2] ]}
    newer = data_and_ver.map{ if (it[0] >= 6.4 ) [ it[1], it[2] ]}

    return([newer, middle, older])   
}
    

/*
* Wrapper for basecalling
*/
workflow BASECALL {

    take: 
        fast5_4_analysis
        
    main:
        switch(params.type) {                      
           case "guppy":
               
               (newer, middle, older) = separateGuppy(fast5_4_analysis) 
    
               outbc65 = GUPPY64_BASECALL(newer)  
               outbc6 = GUPPY6_BASECALL(middle)
               outbc = GUPPY_BASECALL(older)  
               
               basecalled_fastq = outbc.basecalled_fastq.concat(outbc6.basecalled_fastq).concat(outbc65.basecalled_fastq)    
               basecalled_fast5 = outbc.basecalled_fast5.concat(outbc6.basecalled_fast5).concat(outbc65.basecalled_fast5)    
               basecalling_stats = outbc.basecalling_stats.concat(outbc6.basecalling_stats).concat(outbc65.basecalling_stats)    
               break; 
               
           case "dorado": 
               dorado_models = Channel.fromPath("${params.models}", type: "dir", checkIfExists: true)
               outbc = DORADO_BASECALL (fast5_4_analysis, dorado_models)
               basecalled_fastq = outbc.basecalled_fastq
               basecalled_fast5 = Channel.empty()
               basecalling_stats = Channel.empty()
               break; 
        }        
		
    emit:
       basecalled_fast5 = basecalled_fast5
       basecalled_fastq = basecalled_fastq
       basecalling_stats = basecalling_stats
}
