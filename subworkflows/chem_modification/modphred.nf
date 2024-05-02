/*
* Epinano
*/

params.LABEL = ""
params.EXTRAPARS = ""
params.OUTPUT = ""
params.OUTPUTMODE = "copy"
params.CONTAINER = "lpryszcz/modphred-3.6.1:1.0d"

include { FAIDX as SAMTOOLS_FAIDX } from "../misc/samtools" addParams(EXTRAPARS: "", OUTPUT: "")


/*
*/
process getVersion {
    container params.CONTAINER

    output:
	stdout emit: out    
    
    shell:
    """
		/opt/modPhred/run --version      
    """
}


/*
*/

process runByChromosome {

    container params.CONTAINER
    label (params.LABEL)
    tag "Report on ${chr}" 
 	
    input:
    tuple path(fast5), path(bams, stageAs: 'modPhred/minimap2/*'), path(reads, stageAs: 'modPhred/reads/*'), path(reference), path(reffai), val(chr)
    
    output:
    path("modPhred/mod.gz.${chr}.gz")
    
    script:

	"""
	/opt/modPhred/src/mod_report.py -f ${reference} -ri ${fast5} --chr ${chr} -t ${task.cpus}
	"""
}


process mergeOutput {
    if (params.OUTPUT != "") { publishDir(params.OUTPUT, mode:'copy') }

    container params.CONTAINER
    label (params.LABEL)
    tag "Merging all" 
 	
    input:
    tuple path(chrom_res), path(reference), path(fai)
    
    output:
    path("mod.gz")
    
    script:
	"""
		/opt/modPhred/src/merge_chr.py ./mod.gz.*.gz
	"""
}

process runEncodeAndAlign {

    container params.CONTAINER
    label (params.LABEL)
    tag "${sampleID}" 
 	
    input:
    tuple val(sampleID), path(fast5, stageAs: 'input_folder/*'), path(reference)
    
    
    output:
    tuple val(sampleID), path("./modPhred/minimap2/${sampleID}.bam"), emit: bams
    tuple val(sampleID), path("./modPhred/reads/${sampleID}/"), emit: readdir
    tuple val(sampleID), path("./modPhred/reads/${sampleID}/modPhred.pkl"), emit: pkl
    tuple val(sampleID), path("./${sampleID}/"), emit: fast5dir
    
    script:
	"""
	ln -s input_folder ${sampleID}
	/opt/modPhred/src/guppy_encode.py -o ./modPhred/ -ri ./${sampleID}  -t ${task.cpus}
	/opt/modPhred/src/guppy_align.py -f ${reference} -o ./modPhred -ri ./modPhred/reads/* -t ${task.cpus}
	"""
}


/*
*/

workflow RUNBYCHROM {

    take: 
    fast5
    reference
    chr

    main:
    	out_raw = runEncodeAndAlign(fast5.combine(reference))
        fai = SAMTOOLS_FAIDX(reference)

    	all_bams = out_raw.bams.map{it[1]}.collect().map{
 			[it]
 		}
 
     	all_readdir = out_raw.readdir.map{it[1]}.collect().map{
 			[it]
 		}  
 		 	
 		all_fast5_folders = out_raw.fast5dir.map{it[1]}.collect().map{
 			[it]
 		}

        data_for_report = all_fast5_folders.combine(all_bams).combine(all_readdir).combine(reference).combine(fai).combine(chr)
        out_per_chr = runByChromosome(data_for_report)
        data_for_merge = out_per_chr.collect().toList().combine(reference).combine(fai)
        mergeOutput(data_for_merge)
        
}

/*
*/

workflow GET_VERSION {
    main:
		getVersion()
    emit:
    	getVersion.out
}    
