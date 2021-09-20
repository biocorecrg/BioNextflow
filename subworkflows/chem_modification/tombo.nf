/*
* Tombo
*/

params.LABEL = ""
params.EXTRAPARS = ""
params.OUTPUT = ""
params.OUTPUTMODE = "copy"
params.CONTAINER = 	"quay.io/biocontainers/ont-tombo:1.5.1--py37r36h70f9b12_2"


/*
*/
process getVersion {
    container params.CONTAINER

    output:
	stdout emit: out    
    
    shell:
    """
		tombo -v   
    """
}


/*
*/

process resquiggle_rna {
    tag "${idsample}"
    label (params.LABEL)
    container params.CONTAINER
//    if (params.OUTPUT != "") { publishDir(params.OUTPUT, mode:'copy', pattern: '*.bam') }
	
    input:
    tuple val(idsample), path(fast5_folder)
    path(reference)
    
    output:
    tuple val(idsample), path ("*.resquiggle.failed.tsv"), emit: failed_resquiggles
    tuple val(idsample), path  (".*.tombo.index"), emit: tombo_indexes

    script:    
    """ 
    # resquiggling
    tombo resquiggle ${params.EXTRAPARS} ${fast5_folder} ${reference} --rna --processes ${task.cpus} --overwrite --failed-reads-filename ${idsample}.resquiggle.failed.tsv 
    """
}

/*
*/
process getModificationsWithModelSampleCompare {
     tag "${idA} vs ${idB}"

    label (params.LABEL)
	container params.CONTAINER
//  if (params.OUTPUT != "") { publishDir(params.OUTPUT, mode:'copy', pattern: '*.bam') }
            
    input:
    tuple val(idA), val(idB), path(indexA), path(fast5_dirA), path(indexB), path(fast5_dirB)
    file(reference)
    
    output:
    tuple val("${idA}---${idB}"), path("*.bedgraph"), emit: bedgraphs
	tuple val("${idA}---${idB}"), path("*.wig"), emit:  dampened_wiggles
 
    script:
	"""
	mkdir ${idA} ${idB}
	mv ${fast5_dirA} ${indexA} ${idA}
	mv ${fast5_dirB} ${indexB} ${idB}
	
    tombo detect_modifications model_sample_compare ${params.EXTRAPARS} --fast5-basedirs ${idA}/* \
       --control-fast5-basedirs ${idB}/* \
       --processes ${task.cpus} ${params.EXTRAPARS} \
       --statistics-file-basename ${idA}_vs_${idB}_model_sample_compare
      
    tombo text_output browser_files --fast5-basedirs ${idA}/* \
       --control-fast5-basedirs ${idB}/* \
       --browser-file-basename ${idA}_vs_${idB}.features \
       --statistics-filename ${idA}_vs_${idB}_model_sample_compare.tombo.stats \
       --file-types 'dampened_fraction' 'coverage'
    """
}

/*
*/
process getModificationsWithLevelSampleCompare {
    tag "${idA} vs ${idB}"

    label (params.LABEL)
	container params.CONTAINER
//  if (params.OUTPUT != "") { publishDir(params.OUTPUT, mode:'copy', pattern: '*.bam') }
            
    input:
    tuple val(idA), val(idB), path(indexA), path(fast5_dirA), path(indexB), path(fast5_dirB)
    file(reference)
    
    output:
    tuple val("${idA}---${idB}"), path("*.bedgraph"), emit: bedgraphs
	tuple val("${idA}---${idB}"), path("*.wig"), emit: dampened_wiggles
    
    script:
	"""
	mkdir ${idA} ${idB}
	mv ${fast5_dirA} ${indexA} ${idA}
	mv ${fast5_dirB} ${indexB} ${idB}
	
    tombo detect_modifications level_sample_compare ${params.EXTRAPARS} --fast5-basedirs ${idA}/* \
       --alternate-fast5-basedirs ${idB}/* \
       --processes ${task.cpus} ${params.EXTRAPARS} \
       --statistics-file-basename ${idA}_vs_${idB}_level_sample_compare
      
    tombo text_output browser_files --fast5-basedirs ${idA}/* \
       --control-fast5-basedirs ${idB}/* \
       --browser-file-basename ${idA}_vs_${idB}.features \
       --statistics-filename ${idA}_vs_${idB}_level_sample_compare.tombo.stats \
       --file-types {'coverage','statistic'}
    """
}

/*
*/

workflow RESQUIGGLE_RNA {

    take: 
    fast5_folders
    reference
    
    main:
	out = resquiggle_rna(fast5_folders, reference).tombo_indexes

	emit:
	out 

 }

workflow GET_MODIFICATION_MSC {

    take: 
    input_data
    reference
    
    main:   
    getModificationsWithModelSampleCompare(input_data, reference)
	bedgraphs = getModificationsWithModelSampleCompare.out.bedgraphs
	dampened_wiggles = getModificationsWithModelSampleCompare.out.dampened_wiggles

	emit:
		bedgraphs 
		dampened_wiggles
	
 }

workflow GET_MODIFICATION_LSC {

    take: 
    input_data
    reference
    
	main:
    	getModificationsWithLevelSampleCompare(input_data, reference)
		bedgraphs = getModificationsWithLevelSampleCompare.out.bedgraphs
		dampened_wiggles = getModificationsWithLevelSampleCompare.out.dampened_wiggles
	
	emit:
		bedgraphs 
		dampened_wiggles
	
 }

/*
*/

workflow GET_VERSION {
    main:
		getVersion()
    emit:
    	getVersion.out
}    
