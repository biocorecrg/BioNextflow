/*
*  Misc 
*/

params.LABEL = ""
params.EXTRAPARS = ""

params.OUTPUT = "misc_out"
params.CONTAINER = "biocorecrg/centos-perlbrew-pyenv"

include { unzipCmd } from '../global_functions.nf'

process printFileName {
   label (params.LABEL)
    
    tag { id }
    container params.CONTAINER

    input:
    tuple val(id), path(file)
	
    output:
    stdout emit: out

    script:
    """
        echo ${file}
    """
}

process renameFilename {
   label (params.LABEL)
    
    tag { id }
    container params.CONTAINER

    input:
    tuple val(id), path(file), val(newname)
	
    output:
    tuple val(id), path(newname)

    script:
    """
        ln -s ${file} ${newname}
    """
}



// Take first 100 bases
process calcIlluminaAvgReadSize {
    label (params.LABEL)
    
    tag { pair_id }
    container params.CONTAINER

    input:
    tuple val(pair_id), path(reads)
	
    output:
    stdout emit: readsize

    script:
    def first_pair = reads[0]
    """
    if [ `echo ${first_pair} | grep "gz"` ]; then cat="zcat"; else cat="cat"; fi
        \$cat ${first_pair} | awk '{num++}{if (num%4==2){line++; sum+=length(\$0)} if (line==100) {printf "%.0f", (sum/100); exit} }'
    """
}



workflow CALC_AVG_READSIZE {
    take: 
    reads
    type
    
    main:    	
		if (type == "illumina") {
			out = calcIlluminaAvgReadSize(reads.first())
		}
    emit:
    	out
}
 
workflow PRINT_FILE_NAME {
    take: 
    input
    
    main: 
        out = printFileName(input)
       	
   emit:
    	out
}

workflow RENAME_FILE_NAME {
    take: 
    input
    
    main: 
        out = renameFilename(input)
       	
   emit:
    	out
}