/*
*  QC module
*  This workflow allows to make QC on input data
*  It needs input fastq
*/

params.LABEL = ""
params.CONTAINER = "biocorecrg/alevinqc:0.1"
params.OUTPUT = ""

process alevinQC {
    tag "${id}"
    label (params.LABEL)
    container params.CONTAINER
    if (params.OUTPUT != "") { publishDir(params.OUTPUT, mode:'copy') }

    input:
    tuple val(id), path(folder)

    output:
    path("${folder}_QC") 

    script:
	"""
	mkdir ${folder}_QC
	cat > CMD.R << EOL

library("alevinQC")
alevinQCReport(baseDir = "./${folder}",
               sampleId = "${folder}", 
               outputFile = "${folder}.html", 
               outputFormat = "html_document",
               outputDir = "${folder}_QC", forceOverwrite = TRUE)

EOL

	Rscript CMD.R 
	"""
}


workflow QC {
    take: 
    folderin
    
    main:
    out = alevinQC(folderin)
    
    emit:
    out

}

process getVersion {
    container params.CONTAINER

    output:
	stdout emit: out    
    
    shell:
    """
    echo alevinQC' '`Rscript -e "library('alevinQC'); packageVersion('alevinQC')"` 2>/dev/null
    """
}



workflow GET_VERSION {
    main:
		getVersion()
    emit:
    	getVersion.out
}   
 
