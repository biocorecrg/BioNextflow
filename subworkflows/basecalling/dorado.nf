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

if (params.GPU == "OFF") {
	gpu_cmd = '-x "cpu"'
}

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
	if (params.OUTPUT != "") { publishDir(params.OUTPUT, pattern: '*.fastq.gz',  mode: params.OUTPUTMODE ) }

    container params.CONTAINER
             
    input:
    tuple val(idfile), path(fast5), path(models)
    
    output:
    tuple val(idfile), path("*.fastq.gz"), emit: basecalled_fastq

    script:

    """
          dorado basecaller ${gpu_cmd} ${params.EXTRAPARS} --emit-fastq ./ > ${idfile}.fastq
          bgzip -@ ${task.cpus} ${idfile}.fastq
    """
}

process baseCallMod {
    tag { idfile }
    label (params.LABEL)
	if (params.OUTPUT != "") { publishDir(params.OUTPUT, pattern: '*.fastq.gz',  mode: params.OUTPUTMODE ) }

    container params.CONTAINER
             
    input:
    tuple val(idfile), path(fast5), path(models)
    
    output:
    tuple val(idfile), path("*.fastq.gz"), emit: basecalled_fastq

    script:

    """
          dorado basecaller ${gpu_cmd} ${params.EXTRAPARS} ./ > ${idfile}.bam
          samtools fastq -@ ${task.cpus} -T MN,MM,ML ${idfile}.bam > ${idfile}.fastq
          bgzip -@ ${task.cpus} ${idfile}.fastq
          rm *.bam
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

 workflow BASECALLMOD {
    take: 
    input_fast5
    model_folders
    
    main:
        models = model_folders.collect().map{ [ it ] }
    	baseCallMod(input_fast5.combine(models))

	emit:
    	basecalled_fastq = baseCallMod.out.basecalled_fastq
 
}





workflow GET_VERSION {
    main:
		getVersion()
    emit:
    	getVersion.out
} 


