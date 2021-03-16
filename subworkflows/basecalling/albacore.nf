/*
*  Basecalling
*  This workflow allows to make basecalling on ONT fast5 data
*  It needs fast5_files, fast5_type, basecaller, basecaller_opt
*  demultiplexer, demultiplexer_opt, seq_type
*  when included you can specify the GPU_OPTION param ON or OFF for using the GPU
*/

params.GPU_OPTION = ""

process baseCallWithGuppy {
    tag { idfile }

    label (params.GPU_OPTION == "ON" ? 'basecall_gpus': 'basecall_cpus')
            
    input:
    tuple idfile, path(fast5)
    val(multi5)
    val(basecaller_opt)
    val(seq_type) 
    
    output:
    path "${idfile}_out/workspace/*.fast5", emit: basecalled_fast5
    tuple val(idfile), path("${idfile}.*.gz"), emit: basecalled_fastq
    path "${idfile}_out/sequencing_summary.txt", emit: basecalling_stats

    script:
    def RNA_conv_cmd = ""
    def demulti_cmd = ""
    def infolder = "./"
    if (seq_type == "RNA") {    RNA_conv_cmd = " | awk '{if (NR%4==2) gsub(\"U\",\"T\"); print}' " }   

    def multi_cmd = ""
    def gpu_cmd = ""
    def gpu_prefix = ""
    if (params.GPU_OPTION == "ON") {
        gpu_prefix = 'export LD_LIBRARY_PATH="/usr/local/nvidia/lib:/usr/local/nvidia/lib64:/.singularity.d/libs"'
        gpu_cmd = '-x "cuda:0"'
    }
    // in case input files are single fast5 group them in multifast5 at the end
    if (multi5 == 0) {
       multi_cmd = "mkdir single_basecallings temp_multi; mv *_out/workspace/*.fast5 single_basecallings; single_to_multi_fast5 -i single_basecallings -s temp_multi -t ${task.cpus}; mv temp_multi/batch_0.fast5 ./${idfile}_out/workspace/batch_${idfile}.fast5; rm -fr temp_multi single_basecallings"
    }
    """
        ${gpu_prefix}
        guppy_basecaller ${gpu_cmd} --flowcell ${params.flowcell} --kit ${params.kit} --fast5_out ${basecaller_opt} --input ${infolder} --save_path ./${idfile}_out --cpu_threads_per_caller 1  --num_callers  ${task.cpus} 
        cat ${idfile}_out/*.fastq ${RNA_conv_cmd} >> ${idfile}.fastq
        rm ${idfile}_out/*.fastq
        gzip ${idfile}.fastq
        ${multi_cmd}
    """
}

process baseCallWithAlbacore {
    label ('basecall_cpus')
    tag { idfile }
            
    input:
    tuple idfile, path(fast5)
    val(multi5)
    val(basecaller_opt)
    val(seq_type) 
    
    output:
    path "${idfile}_out/workspace/*.fast5", emit: basecalled_fast5
    tuple val(idfile), path("${idfile}.*.gz"), emit: basecalled_fastq
    path "${idfile}_out/sequencing_summary.txt", emit: basecalling_stats

    script:
    // conversion command if input is RNA - have to check if this is really needed
    def RNA_conv_cmd = ""
    def demulti_cmd = ""
    def infolder = "./"
    if (seq_type == "RNA") {    RNA_conv_cmd = " | awk '{if (NR%4==2) gsub(\"U\",\"T\"); print}' " }   
    // in case input files are multi fast5 convert them in single fast5 since albacore is not able to deal with multi fast5
    if (multi5 == 1) {
        demulti_cmd = "mkdir demulti_tmp; mkdir demulti; multi_to_single_fast5 -i ${infolder} -s demulti_tmp -t ${task.cpus}; mv demulti_tmp/*/*.fast5 demulti; rm -fr demulti_tmp"
        infolder = "demulti"
    }
    """
    ${demulti_cmd}
    export PYTHONPATH=$baseDir/bin/albacore:\$PYTHONPATH
    read_fast5_basecaller.py ${basecaller_opt} --flowcell \"${params.flowcell}\" --kit \"${params.kit}\" --output_format fastq,fast5 \
        --worker_threads ${task.cpus} -s ./${idfile}_out --disable_filtering --input ${infolder};
    cat ${idfile}_out/workspace/*.fastq ${RNA_conv_cmd} >> ${idfile}.fastq
    rm ${idfile}_out/workspace/*.fastq
    gzip ${idfile}.fastq
    mkdir single_basecallings
    mv ${idfile}_out/workspace/*/*.fast5 single_basecallings
    mkdir temp_multi
    single_to_multi_fast5 -i single_basecallings -s temp_multi -t ${task.cpus}
    mv temp_multi/batch_0.fast5 ./${idfile}_out/workspace/batch_${idfile}.fast5
    if [-d demulti]; then rm -fr demulti; fi
    rm -fr single_basecallings temp_multi
    """
}



 workflow BASECALL_ONT {
    take: 
    fast5_files
    fast5_type
    basecaller
    basecaller_opt
    seq_type
    
    main:
    if (basecaller == "guppy") {
        baseCallWithGuppy(fast5_files, fast5_type, basecaller_opt, seq_type)
    }
    if (basecaller == "albacore") {
        baseCallWithAlbacore(fast5_files, fast5_type, basecaller_opt, seq_type)
    }


}




