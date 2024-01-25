/*
* BASECALLING 
*/

moduleFolder   = "../subworkflows"

params.gpu = ""
params.outmode = "copy"
params.output = ""
params.label = ""
params.extrapars = ""
params.models = ""

def cuda_cont = (params.gpu == 'cuda11' ? 'biocorecrg/mopbasecallc11:0.3' : 'biocorecrg/mopbasecall:0.3')

include { GET_VERSION as GUPPY_VERSION } from "${moduleFolder}/basecalling/guppy" 
include { BASECALL as GUPPY_BASECALL } from "${moduleFolder}/basecalling/guppy" addParams(LABEL: params.label, EXTRAPARS_BC: params.extrapars["guppy"], GPU: params.gpu, MOP: "YES", CONTAINER: cuda_cont)
include { BASECALL as GUPPY6_BASECALL } from "${moduleFolder}/basecalling/guppy" addParams(VERSION:"6", EXTRAPARS_BC: params.extrapars["guppy"], LABEL: params.label, GPU: params.gpu, MOP: "YES", OUTPUT: params.output, CONTAINER: cuda_cont, OUTPUTMODE: params.outmode)
include { BASECALL as GUPPY65_BASECALL } from "${moduleFolder}/basecalling/guppy" addParams(VERSION:"6.5", EXTRAPARS_BC: params.extrapars["guppy"], LABEL: params.label, GPU: params.gpu, MOP: "YES", OUTPUT: params.output, CONTAINER: cuda_cont, OUTPUTMODE: params.outmode)
include { BASECALL as DORADO_BASECALL } from "${moduleFolder}/basecalling/dorado" addParams(EXTRAPARS:  params.extrapars["dorado"], LABEL: params.label, GPU: params.gpu, MOP: "YES", OUTPUT: params.output, OUTPUTMODE: params.outmode)

// ADD A CHECK FOR GUPPY FOR DISABLING SCORE

def separateGuppy (fast5) {

	data_and_ver = GUPPY_VERSION().map{
        def vals = it.split("\\.")
    	"${vals[0]}.${vals[1]}"
    }.toBigDecimal().combine(fast5)
	
	
    older = data_and_ver.map{ if (it[0] < 6 ) [ it[1], it[2] ]}
    middle = data_and_ver.map{ if (it[0] >= 6 && it[0] < 6.5) [ it[1], it[2] ]}
    newer = data_and_ver.map{ if (it[0] >= 6.5 ) [ it[1], it[2] ]}

    return([newer, middle, older])   
}
    

/*
* Wrapper for basecalling
*/
workflow BASECALL {

    take: 
        fast5_4_analysis
        
    main:
        switch(params.basecalling) {                      
           case "guppy":
           
               (newer, middle, older) = separateGuppy(fast5_4_analysis) 
    
               outbc65 = GUPPY65_BASECALL(newer)  
               outbc6 = GUPPY6_BASECALL(middle)
               outbc = GUPPY_BASECALL(older)  
               
               basecalled_fastq = outbc.basecalled_fastq.concat(outbc6.basecalled_fastq).concat(outbc65.basecalled_fastq)    
               basecalled_fast5 = outbc.basecalled_fast5.concat(outbc6.basecalled_fast5).concat(outbc65.basecalled_fast5)    
               basecalling_stats = outbc.basecalling_stats.concat(outbc6.basecalling_stats).concat(outbc65.basecalling_stats)    
               break; 
           case "dorado": 
               dorado_models = Channel.fromPath("${params.models}/*", type: "dir")
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
