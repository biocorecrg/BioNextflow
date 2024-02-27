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
params.CONTAINER = "biocorecrg/mopbasecall:0.3"
params.GPU = ""
params.VERSION = "5"

def EXTRAPARS_BC = params.EXTRAPARS_BC
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
		guppy_basecaller --version | awk -F "Version" '{print \$2}' | awk -F"," '{print \$1}'| awk -F"+" '{print \$1}'| head -n 1      
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

process baseCallNew {

   tag { idfile }
    label (params.LABEL)

    container params.CONTAINER
             
    input:
    tuple val(idfile), path(fast5)
    
    output:
    tuple val(idfile), path("*.fastq.gz"), emit: basecalled_fastq

    script:
    def infolder = "./"

    """
        ${library_export}
        guppy_basecaller ${gpu_cmd} \
        ${EXTRAPARS_BC} -i ${infolder} \
        -s ./${idfile}_out \
        --gpu_runners_per_device 1 \
        --cpu_threads_per_caller 1 \
        --disable_qscore_filtering \
        --num_callers  ${task.cpus}
        cat ${idfile}_out/*.fastq >> ${idfile}.fastq
        rm ${idfile}_out/*.fastq
        bgzip -@ ${task.cpus} ${idfile}.fastq
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
        --fast5_out ${EXTRAPARS_BC} -i ${infolder} \
        --save_path ./${idfile}_out \
        --gpu_runners_per_device 1 \
        --cpu_threads_per_caller 1 \
	    --num_callers  ${task.cpus}
        cat ${idfile}_out/*.fastq >> ${idfile}.fastq
        rm ${idfile}_out/*.fastq
        bgzip -@ ${task.cpus} ${idfile}.fastq
    """
}

process baseCall_mod {
    tag { idfile }
    label (params.LABEL)
	if (params.MOP == "YES")  { if (params.OUTPUT != "") { 
		publishDir(params.OUTPUT, pattern: '*_out/workspace/*.fast5',  mode: params.OUTPUTMODE, saveAs: { file -> "${idfile.split('---')[0]}/${file.split('\\/')[-1]}" } ) }
	}
		
	else if (params.OUTPUT != "") { publishDir(params.OUTPUT, pattern: '*_out/workspace/*.fast5',  mode: params.OUTPUTMODE, saveAs: { file -> "${idfile}/${file.split('\\/')[-1]}" } ) }

    container params.CONTAINER
             
    input:
    tuple val(idfile), path(fast5)
    
    output:
    tuple val(idfile), path("${idfile}_out/workspace/*.fast5"), emit: basecalled_fast5
    tuple val(idfile), path("*.fastq.gz"), emit: basecalled_fastq
    tuple val(idfile), path("${idfile}_out/sequencing_summary.txt"), emit: basecalling_stats
    tuple val(idfile), path("*.bam"), emit: basecalled_aln

    script:
    def infolder = "./"

    """
        ${library_export}
        guppy_basecaller ${gpu_cmd} \
        --fast5_out ${EXTRAPARS_BC} -i ${infolder} \
        --bam_out \
        --save_path ./${idfile}_out \
        --gpu_runners_per_device 1 \
        --cpu_threads_per_caller 1 \
	    --num_callers  ${task.cpus}
        cat ${idfile}_out/*.fastq >> ${idfile}.fastq
        mv ${idfile}_out/*.bam ${idfile}.bam
        rm ${idfile}_out/*.fastq
        bgzip -@ ${task.cpus} ${idfile}.fastq
    """
}

process baseCall_aln {
    tag { idfile }
    label (params.LABEL)
	if (params.MOP == "YES")  { if (params.OUTPUT != "") { 
		publishDir(params.OUTPUT, pattern: '*_out/workspace/*.fast5',  mode: params.OUTPUTMODE, saveAs: { file -> "${idfile.split('---')[0]}/${file.split('\\/')[-1]}" } ) } 
	
	}
		
	else if (params.OUTPUT != "") { publishDir(params.OUTPUT, pattern: '*_out/workspace/*.fast5',  mode: params.OUTPUTMODE, saveAs: { file -> "${idfile}/${file.split('\\/')[-1]}" } ) }

    container params.CONTAINER
             
    input:
    tuple val(idfile), path(fast5)
    path(reference)
    
    output:
    tuple val(idfile), path("${idfile}_out/workspace/*.fast5"), emit: basecalled_fast5
    tuple val(idfile), path("*.fastq.gz"), emit: basecalled_fastq
    tuple val(idfile), path("${idfile}_out/sequencing_summary.txt"), emit: basecalling_stats
    tuple val(idfile), path("*.bam"), emit: basecalled_aln

    script:
    def infolder = "./"

    """
        ${library_export}
        guppy_basecaller ${gpu_cmd} \
        --fast5_out ${EXTRAPARS_BC} -i ${infolder} \
        --align_ref ${reference} \
        --bam_out \
        --save_path ./${idfile}_out \
        --gpu_runners_per_device 1 \
        --cpu_threads_per_caller 1 \
	    --num_callers  ${task.cpus}
        cat ${idfile}_out/*.fastq >> ${idfile}.fastq
        mv ${idfile}_out/*.bam ${idfile}.bam
        rm ${idfile}_out/*.fastq
        bgzip -@ ${task.cpus} ${idfile}.fastq
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
		${EXTRAPARS_BC} ${params.EXTRAPARS_DEM} \
		--num_barcode_threads ${task.cpus} \
		--trim_barcodes  \
		--fast5_out -i ${infolder} \
		--save_path ./${idfile}_out \
        --gpu_runners_per_device 1 \
		--cpu_threads_per_caller 1  \
		--num_callers ${task.cpus}
		cd ${idfile}_out; 
		if [ "\$(find . -type d -name "barcode*" )" != "" ] 
		then
			for d in barcode*; do echo \$d; cat \$d/*.fastq > ../${idfile}.\$d.fastq; done;
		fi
		if [ -d unclassified ]; then
			cat unclassified/*.fastq > ../${idfile}.unclassified.fastq; cd ../
		fi
		for i in *.fastq; do bgzip -@ ${task.cpus} \$i; done	
        rm *_out/*/*.fastq
     """
}


 workflow BASECALL {
 
    take: 
    input_fast5
    
    main:
       switch(params.VERSION) {                      
           case "6":
                EXTRAPARS_BC = EXTRAPARS_BC + " --disable_qscore_filtering "
           case "5":
				baseCall(input_fast5)
     			basecalled_fast5 = baseCall.out.basecalled_fast5
    			basecalled_fastq = baseCall.out.basecalled_fastq
    			basecalling_stats = baseCall.out.basecalling_stats
           break;
           case "6.5":
 				baseCallNew(input_fast5)
     			basecalled_fast5 = channel.empty()
    			basecalled_fastq = baseCallNew.out.basecalled_fastq
    			basecalling_stats = channel.empty()
           break;           
    	}
    	
	emit:
    	basecalled_fast5 = basecalled_fast5
    	basecalled_fastq = basecalled_fastq
    	basecalling_stats = basecalling_stats
 
}

 workflow BASECALL_ALN {
    take: 
    input_fast5
    reference
    
    main:
    	baseCall_aln(input_fast5, reference)

	emit:
    	basecalled_fast5 = baseCall_aln.out.basecalled_fast5
    	basecalled_fastq = baseCall_aln.out.basecalled_fastq
    	basecalling_stats = baseCall_aln.out.basecalling_stats
    	basecalled_aln = baseCall_aln.out.basecalled_aln
 
}

workflow BASECALL_MOD {
    take: 
    input_fast5
    
    main:
    	baseCall_mod(input_fast5)

	emit:
    	basecalled_fast5 = baseCall_mod.out.basecalled_fast5
    	basecalled_fastq = baseCall_mod.out.basecalled_fastq
    	basecalling_stats = baseCall_mod.out.basecalling_stats
    	basecalled_aln = baseCall_mod.out.basecalled_aln
 
}


 workflow BASECALL_DEMULTI {
    take: 
    input_fast5
    
    main:
       switch(params.VERSION) {                      
           case "6":
                EXTRAPARS_BC = EXTRAPARS_BC + " --disable_qscore_filtering "
           case "5":
				baseCallAndDemultiPlex(input_fast5)
     			basecalled_fast5 = baseCallAndDemultiPlex.out.basecalled_fast5
    			basecalled_fastq = baseCallAndDemultiPlex.out.basecalled_fastq
    			basecalling_stats = baseCallAndDemultiPlex.out.basecalling_stats
           break;
           case "6.5":
 		//		baseCallNew(input_fast5)
     			basecalled_fast5 = channel.empty()
    			basecalled_fastq = channel.empty()
    			basecalling_stats = channel.empty()
           break;           
    	}

	emit:
    	basecalled_fast5 = basecalled_fast5
    	basecalled_fastq = basecalled_fastq
    	basecalling_stats = basecalling_stats
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
