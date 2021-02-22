/*
*  MoAIMS 
*/

params.LABEL = ""
params.EXTRAPARS = ""

params.OUTPUT = "moaims_out"
params.CONTAINER = "biocorecrg/moaims:1.0"

include { unzipCmd } from '../global_functions.nf'

process getVersion {
    container params.CONTAINER

    output:
	stdout emit: out    
    
    shell:
    """
    Rscript -e "library('moaims'); packageVersion('moaims')" 2>/dev/null
    """
}


process peakCall {
    label (params.LABEL)
    tag { comp_id }
    container params.CONTAINER
    publishDir(params.OUTPUT, mode:'copy')

    input:
    tuple val(comp_id), path(sample), path(input)
    path(annofile)

    output:
    tuple val(pair_id), path("*.bed") 
    
	script:
	def filecontent = "SampleID\tBamIP\tBamInput"
    def unzip_anno = unzipCmd(annofile)
    def cmd_anno = unzip_anno[1]
    def anno_name = unzip_anno[0]
    def cmd_clean = unzip_anno[2]

    """
	${cmd_anno}
	echo "${filecontent}" > sample_info_file
	echo -e "${comp_id}\t\$PWD/${sample}\t\$PWD/${input}" >> sample_info_file 
	echo 'library(moaims)' > R.cmd;
	echo 'moaims(sample_info_file = "'\$PWD/sample_info_file'", gtf_file ="'\$PWD/${anno_name}'", strand_specific = 1, is_paired = F, proj_name="'${comp_id}'")' >> R.cmd   
	Rscript R.cmd
	${cmd_clean}

    """
}



workflow MOAIMS_CALL {
    take: 
    comparisons
    annotation
    
    main:
    	
		out = peakCall(comparisons, annotation)
    emit:
    	out
}


workflow GET_VERSION {
    main:
		getVersion()
    emit:
    	getVersion.out
}    
