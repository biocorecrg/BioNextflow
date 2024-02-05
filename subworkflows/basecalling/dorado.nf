/*
*  This workflow allows to wrap the DORADO program for basecalling on ONT fast5 data
*  when included you can specify the GPU param ON or OFF for using the GPU
*/

params.LABEL = ""
params.EXTRAPARS_BC = ""
params.EXTRAPARS_DEM = ""
params.OUTPUT = ""
params.OUTPUTMODE = "copy"
params.MOP = ""
params.CONTAINER = "ontresearch/dorado:shab1ff19616e2b8635791f17bef11f806628505a35"
params.GPU = ""

def gpu_cmd = ""
def library_export = ""


process getVersion {
    container params.CONTAINER
    label (params.LABEL)

    output:
	stdout emit: out    
    
    shell:
    """
		dorado --version 2>&1 | head -n1    
    """
}

process baseCall {
    tag { idfile }
    label (params.LABEL)
	if (params.MOP == "YES")  { if (params.OUTPUT != "") { publishDir(params.OUTPUT, pattern: '*_out/workspace/*.fast5',  mode: params.OUTPUTMODE, saveAs: { file -> "${idfile.split('---')[0]}/${file.split('\\/')[-1]}" } ) } }
	else if (params.OUTPUT != "") { publishDir(params.OUTPUT, pattern: '*_out/workspace/*.fast5',  mode: params.OUTPUTMODE, saveAs: { file -> "${idfile}/${file.split('\\/')[-1]}" } ) }

    container params.CONTAINER
             
    input:
    tuple val(idfile), path(fast5), path(models)
    
    output:
    tuple val(idfile), path("*.fastq.gz"), emit: basecalled_fastq

    script:

    """
          dorado basecaller ${params.EXTRAPARS} --emit-fastq ./ > ${idfile}.fastq
          bgzip -@ ${task.cpus} ${idfile}.fastq
    """
}

 workflow BASECALL {
    take: 
    input_fast5
    model_folders
    
    main:
        models = model_folders.collect().map{ [ it ] }
    	baseCall(input_fast5.combine(models))

	emit:
    	basecalled_fastq = baseCall.out.basecalled_fastq
 
}





workflow GET_VERSION {
    main:
		getVersion()
    emit:
    	getVersion.out
} 

workflow GET_WORKFLOWS {
    take: 
    flowcell
    kit

    main:
		getWorkflow(flowcell, kit)
    emit:
    	getWorkflow.out
} 
