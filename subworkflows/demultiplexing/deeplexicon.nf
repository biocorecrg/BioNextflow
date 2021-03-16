/*
*  This workflow allows to wrap the deeplexicon program for demuliplexing ONT fast5 data
*  when included you can specify the GPU param ON or OFF for using the GPU
*/

params.LABEL = ""
params.EXTRAPARS
params.OUTPUT = ""
params.CONTAINER = "biocorecrg/mopdem:0.2"
params.GPU = ""

def gpu_cmd = ""
//def library_export = ""

//if (params.GPU == "ON") {
//	library_export = 'export LD_LIBRARY_PATH="/usr/local/nvidia/lib:/usr/local/nvidia/lib64:/.singularity.d/libs"'
//}


process demultiplex {
    tag { idfile }
    label (params.LABEL)
    //if (params.OUTPUT != "") { publishDir(params.OUTPUT, pattern: '*_out/workspace/*.fast5',  mode: 'move', saveAs: { file -> "${idfile}/${file.split('\\/')[-1]}" } ) }

    container params.CONTAINER
             
    input:
    tuple val(idfile), path(fast5)
    path(deeplexicon_folder) 

    output:
		tuple val(idfile), path("${idfile}_demux.tsv"), emit: demux_files
 
    script:  
    """
		ln -s ${deeplexicon_folder}/* .
        deeplexicon.py -p ./ ${params.EXTRAPARS} -f multi -b 4000 -v > ${idfile}_demux.tsv
    """
}

	




 workflow DEEPLEXICON_DEMULTI {
    take: 
    input_fast5
    deep_folder
    
    main:
    	demultiplex(input_fast5, deep_folder)

	emit:
    	demux_files = demultiplex.out.demux_files
  
}

 
