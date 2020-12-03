/*
*  QC module
*  This workflow allows to get some info on fastq data
*/

params.LABEL = ""
params.CONTAINER = ""

process getReadLen {
    tag { fastq }
    label (params.LABEL)

    input:
    path(fastq)

    output:
	stdout()


    script:
	"""
       if [ `echo ${fastq} | grep "gz"` ]; then cat="zcat"; else cat="cat"; fi
        \$cat ${fastq} | awk '{num++}{if (num%4==2){line++; sum+=length(\$0)} if (line==100) {printf "%.0f", sum/100; exit} }'
	"""
}

workflow FASTQ_STATS {
    take: 
    fastq
    stat
    
    main:
	if (stat == "get_read_len") {
		seq_len_raw = getReadLen(fastq)
    	out = seq_len_raw.map {  it.trim().toInteger() }
    	out.map{"The average read length is ${it}" }.view()
	}
    emit:
        out
}

 
