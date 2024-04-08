/*
* Epinano
*/

params.LABEL = ""
params.EXTRAPARS = ""
params.OUTPUT = ""
params.OUTPUTMODE = "copy"
params.CONTAINER = "lpryszcz/modphred-3.6.1"

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
    tag "${sampleID} on ${chr}" 
 	
    input:
    tuple val(sampleID), path(fast5), path(reference), val(chr)
    
    output:
    tuple val(sampleID), path("./${sampleID}_${chr}"), emit: outfolder
    tuple val(sampleID), path("./${sampleID}_${chr}/mod.gz.${chr}.gz"), emit: outfile
    
    script:
	"""
	/opt/modPhred/run -f ${reference} -o ./${sampleID}_${chr} -i ./ --chr ${chr} -t ${task.cpus}
	"""
}


process mergeOutputPerChr {
    if (params.OUTPUT != "") { publishDir(params.OUTPUT, mode:'copy', pattern: '*.csv.gz') }

    container params.CONTAINER
    label (params.LABEL)
    tag "${sampleID} on ${chr}" 
 	
    input:
    tuple val(sampleID), path("*"), path(reference), val(chr)
    
    output:
    tuple val(sampleID), val(chr), path("./${sampleID}_${chr}")
    
    script:
	"""
		/opt/modPhred/src/merge_chr.py  output_dir/mod.gz.chr*.gz
	"""
}

process runEncodeAndAlign {

    container params.CONTAINER
    label (params.LABEL)
    tag "${sampleID}" 
 	
    input:
    tuple val(sampleID), path(fast5), path(reference)
    
    output:
    tuple val(sampleID), path("./modPhred/minimap2/${sampleID}.bam"), emit: bams
    tuple val(sampleID), path("./modPhred/reads/${sampleID}/"), emit: readdir
    tuple val(sampleID), path("./modPhred/reads/${sampleID}/modPhred.pkl"), emit: pkl
    
    
    script:
	"""
	mkdir ${sampleID}; cd ${sampleID}; ln -s ../*.fast5 ./; cd ../
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
    	out = runEncodeAndAlign(fast5.combine(reference))
    	out.pkl.view()
        //data = fast5.combine(reference).combine(chr)
        //out_per_chr = runByChromosome(data).outfile
        //out_per_chr.groupTuple().view()

}

/*
*/

workflow GET_VERSION {
    main:
		getVersion()
    emit:
    	getVersion.out
}    
