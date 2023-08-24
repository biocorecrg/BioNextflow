/*
*  Misc 
*/

params.LABEL = ""
params.EXTRAPARS = ""

params.OUTPUT = "misc_out"
params.CONTAINER = "biocorecrg/centos-perlbrew-pyenv"

include { unzipCmd } from '../global_functions.nf'
include { zcatOrCat } from '../global_functions.nf'

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

// Take first 100 bases
process downSamplePairs {
    
    tag { "${id}" }
    container params.CONTAINER

    input:
    tuple val(id), path(reads) 
	val(readnum)
	
    output:
    tuple val(id), path("${id}_sub_*.fq.gz")

    script:
    def readA = reads[0]
    def cmdA = zcatOrCat(readA)
    def rownum = readnum*4

    def readB = reads[1]
    def cmdB = zcatOrCat(readB)
    
	"""
		${cmdA} | head -n ${rownum} | gzip > ${id}_sub_1.fq.gz
		${cmdB} | head -n ${rownum} | gzip > ${id}_sub_2.fq.gz
    """
}

// Take first 100 bases
process PossiblyUnzipGenome {
    label (params.LABEL)
    
    tag { "${genome}" }
    container params.CONTAINER

    input:
    path(genome)
	
    output:
    path("outgenome.fa")

    script:
	"""
	if [[ ${genome} == *.gz ]] 
	then
		zcat ${genome}  > outgenome.fa
	else
		cp ${genome} outgenome.fa
	fi
    """
}

workflow DOWNSAMPLE_PAIRS {

    take: 
    reads
    value
    
    main:   	
		out = downSamplePairs(reads, value)
    emit:
    	out
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
