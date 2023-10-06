/*
* trax module
*/


params.CONTAINER = "ucsclowelab/trax:latest"
params.OUTPUT = ""

params.LABEL = ""
params.EXTRAPARS = ""

include { checkSEorPE } from '../global_functions.nf' 

process getVersion {
    container params.CONTAINER

    output:
    stdout emit: out

    shell:
    """
    echo \"trax v1.0\"
    """
}

process quickdb {
    if (params.OUTPUT != "") {publishDir(params.OUTPUT, mode: 'copy') }
   
    tag {tagDB }
    container params.CONTAINER
    label (params.LABEL)

    input:
    val (tagDB)

    output:
	path("./db")

    script:

    """
	quickdb2.0.bash ${tagDB} ./db ${task.cpus}
    """
}

process trimAdapters {
    if (params.OUTPUT != "") {publishDir(params.OUTPUT, mode: 'copy') }
   
    tag ("All")
    
    container params.CONTAINER
    label (params.LABEL)

    input:
    path(runfile)
    path(reads)
    val(type)

    output:
	path("*_{merge,trimmed}.fastq.gz"), emit: res

    script:
    def par = ""
    if ("${type}" == "SE") {
    	par = "--singleend"
    }
    """
	trimadapters.py --cores ${params.EXTRAPARS} ${task.cpus} ${par} --runname "runname" --runfile ${runfile} 
    """
}


process processSamples {
    if (params.OUTPUT != "") {publishDir(params.OUTPUT, mode: 'copy') }
   
    tag ("All")
    
    container params.CONTAINER
    label (params.LABEL)

    input:
    path(reads)
    path(db)
    path(samplefile)
    path(samplepairs)

//    output:


    script:
    """
    	mkdir temp
    	export TMPDIR=\$PWD/temp
		processsamples.py --experimentname=expname --databasename=${db}/db \
		--ensemblgtf=${db}/genes.gtf --samplefile=${samplefile} \
		--exppairs=${samplepairs} ${params.EXTRAPARS}
    """
}


workflow QUICKDB {
    take: 
    tagDB
    
    main:
    out = quickdb(tagDB) 

    emit:
  	out
}

workflow TRIMADAPTERS {
    take: 
    reads
    
    main:
    reads.map{
		"${it[0]}" + "\t" + it[1].join("\t")
	}.set{read_list}

	runfile = read_list.collectFile(name: 'runfile.txt', newLine: true)
	files = reads.map{it[1]}.collect()
	type = checkSEorPE(reads)	
    out = trimAdapters(runfile, files, type) 

    emit:
  	out
}

workflow PROCESSSAMPLES {
    take: 
    trimmedreads
    db
    samplefile
    samplepairs
    
    
    main:
    samplecont = samplefile.splitCsv(sep: "\t")
    trim_index = trimmedreads.flatten().map{
    	["${it.getName()}".replaceAll("_merge.fastq.gz","").replaceAll("_trimmed.fastq.gz",""), it.getName()]
    }

    final_samplefile = samplecont.join(trim_index).map{
		"${it[0]}\t${it[1]}\t${it[2]}"
	}.collectFile(name: 'samplefile.txt', newLine: true)
    
    processSamples(trimmedreads, db, final_samplefile, samplepairs) 


}




workflow GET_VERSION {
    main:
		getVersion()
    emit:
    	getVersion.out
} 


