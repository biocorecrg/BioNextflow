/*
*  This workflow allows to wrap the GUPPY program for basecalling on ONT fast5 data
*  when included you can specify the GPU param ON or OFF for using the GPU
*/

params.LABEL = ""
params.EXTRAPARS_BC = ""
params.EXTRAPARS_DEM = ""
params.OUTPUT = ""
params.OUTPUTMODE = "copy"
params.MOP = ""
params.CONTAINER = "biocorecrg/mopbasecall:0.2"
params.GPU = ""

def gpu_cmd = ""
def library_export = ""

if (params.GPU == "ON") {
	gpu_cmd = '-x "cuda:0"'
	library_export = 'export LD_LIBRARY_PATH="/usr/local/nvidia/lib:/usr/local/nvidia/lib64:/.singularity.d/libs"'
}

process getVersion {
    container params.CONTAINER
    label (params.LABEL)

    output:
	stdout emit: out    
    
    shell:
    """
		guppy_basecaller --version       
    """
}

process getWorkflow {
    container params.CONTAINER

	input:
	val(flowcell) 
    val(kit)

    output:
	stdout emit: out    
    
    shell:
    """
		guppy_basecaller --version
		guppy_basecaller --print_workflows | grep ${flowcell} | grep ${kit}
    """
}


process baseCall {
    tag { idfile }
    label (params.LABEL)
	if (params.MOP == "YES")  { if (params.OUTPUT != "") { publishDir(params.OUTPUT, pattern: '*_out/workspace/*.fast5',  mode: params.OUTPUTMODE, saveAs: { file -> "${idfile.split('---')[0]}/${file.split('\\/')[-1]}" } ) } }
	else if (params.OUTPUT != "") { publishDir(params.OUTPUT, pattern: '*_out/workspace/*.fast5',  mode: params.OUTPUTMODE, saveAs: { file -> "${idfile}/${file.split('\\/')[-1]}" } ) }

    container params.CONTAINER
             
    input:
    tuple val(idfile), path(fast5)
    
    output:
    tuple val(idfile), path("${idfile}_out/workspace/*.fast5"), emit: basecalled_fast5
    tuple val(idfile), path("*.fastq.gz"), emit: basecalled_fastq
    tuple val(idfile), path("${idfile}_out/sequencing_summary.txt"), emit: basecalling_stats

    script:
    def infolder = "./"

    """
        ${library_export}
        guppy_basecaller ${gpu_cmd} \
        --fast5_out ${params.EXTRAPARS_BC} -i ${infolder} \
        --save_path ./${idfile}_out \
        --gpu_runners_per_device 1 \
        --cpu_threads_per_caller 1 \
	    --num_callers  ${task.cpus}
        cat ${idfile}_out/*.fastq >> ${idfile}.fastq
        rm ${idfile}_out/*.fastq
        gzip ${idfile}.fastq
    """
}

process baseCallAndDemultiPlex {
    tag { idfile }
    label (params.LABEL)
	if (params.MOP == "YES")  { if (params.OUTPUT != "") { publishDir(params.OUTPUT, pattern: '*_out/workspace/*.fast5',  mode: params.OUTPUTMODE, saveAs: { file -> "${idfile.split('---')[0]}/${file.split('\\/')[-1]}" } ) } }
	else if (params.OUTPUT != "") { publishDir(params.OUTPUT, pattern: '*_out/workspace/*.fast5',  mode: params.OUTPUTMODE, saveAs: { file -> "${idfile}/${file.split('\\/')[-1]}" } ) }
    
    container params.CONTAINER
             
    input:
    tuple val(idfile), path(fast5)
    
    output:
    tuple val(idfile), path("${idfile}_out/workspace/*.fast5"), emit: basecalled_fast5
    tuple val(idfile), path("${idfile}.*.gz"), emit: basecalled_fastq
    tuple val(idfile), path("${idfile}_out/sequencing_summary.txt"), emit: basecalling_stats

    script:
    def infolder = "./"

    """
     ${library_export}
		guppy_basecaller ${gpu_cmd} \
		${params.EXTRAPARS_BC} ${params.EXTRAPARS_DEM} \
		--num_barcode_threads ${task.cpus} \
		--trim_barcodes  \
		--fast5_out -i ${infolder} \
		--save_path ./${idfile}_out \
        --gpu_runners_per_device 1 \
		--cpu_threads_per_caller 1  \
		--num_callers ${task.cpus}
		cd ${idfile}_out; 
		if [ -d barcode01 ]; then
			for d in barcode*; do echo \$d; cat \$d/*.fastq > ../${idfile}.\$d.fastq; done;
		fi
		if [ -d unclassified ]; then
			cat unclassified/*.fastq > ../${idfile}.unclassified.fastq; cd ../
		fi
		for i in *.fastq; do gzip \$i; done	
        rm *_out/*/*.fastq
     """
}


 workflow BASECALL {
    take: 
    input_fast5
    
    main:
    	baseCall(input_fast5)

	emit:
    	basecalled_fast5 = baseCall.out.basecalled_fast5
    	basecalled_fastq = baseCall.out.basecalled_fastq
    	basecalling_stats = baseCall.out.basecalling_stats
 
}

 workflow BASECALL_DEMULTI {
    take: 
    input_fast5
    
    main:
    	baseCallAndDemultiPlex(input_fast5)

	emit:
    	basecalled_fast5 = baseCallAndDemultiPlex.out.basecalled_fast5
    	basecalled_fastq = baseCallAndDemultiPlex.out.basecalled_fastq
    	basecalling_stats = baseCallAndDemultiPlex.out.basecalling_stats
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
