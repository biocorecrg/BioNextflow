/*
* NanoPolish
*/

params.LABEL = ""
params.EXTRAPARS = ""
params.OUTPUT = ""
params.OUTPUTMODE = "copy"
params.CONTAINER = 	"quay.io/biocontainers/nanocompore:1.0.0rc3.post1--py38_0"


/*
*/
process getVersion {
    container params.CONTAINER

    output:
	stdout emit: out    
    
    shell:
    """
		 echo nanocompore `nanocompore -v`      
    """
}


/*
*/

process sampleCompare {
    if (params.OUTPUT != "") { publishDir(params.OUTPUT, mode:'copy') }

    container params.CONTAINER
    label (params.LABEL)
    tag "${folder_name_A} vs ${folder_name_B}" 
 		
    input:
    tuple val(folder_name_A), val(folder_name_B), path(tsv_file_A), path(tsv_file_B)
    file(reference)
    
    output:
    file("${folder_name_A}-${folder_name_B}_nanocompore")
    
    script:
	"""
	nanocompore sampcomp --nthreads ${task.cpus}\
    --file_list1 ${tsv_file_A}/out_eventalign_collapse.tsv \
    --file_list2 ${tsv_file_B}/out_eventalign_collapse.tsv \
    --label1 ${folder_name_A} \
    --label2 ${folder_name_B} \
    --fasta ${reference} \
    --outpath ./${folder_name_A}-${folder_name_B}_nanocompore/ \
    ${params.EXTRAPARS} --pvalue_thr 1 --outprefix ${folder_name_A}_vs_${folder_name_B} --logit --comparison_methods GMM,KS,MW,TT --overwrite 
	"""


}

/*
*/

workflow SAMPLE_COMPARE {

    take: 
    input
    reference
    
    main:
	out = sampleCompare(input, reference)

	emit:
	out
 }


/*
*/

workflow GET_VERSION {
    main:
		getVersion()
    emit:
    	getVersion.out
}    
