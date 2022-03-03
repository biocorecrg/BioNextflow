/*
*  MoAIMS 
*/
 
params.LABEL = ""
params.EXTRAPARS = ""

params.OUTPUT = ""
params.CONTAINER = "biocorecrg/exomepeak2:2.0"

include { unzipCmd } from '../global_functions.nf'

process getVersion {
    container params.CONTAINER

    output:
	stdout emit: out    
    
    shell:
    """
    echo exomePeak2' '`Rscript -e "library('exomePeak2'); packageVersion('exomePeak2')"` 2>/dev/null
    """
}


process peakCall {
    label (params.LABEL)
    tag { comp_id }
    container params.CONTAINER
    if (params.OUTPUT != "") { publishDir(params.OUTPUT, mode:'copy') }

    input:
    tuple val(comp_id), path(sample), path(input)
    path(annofile)

    output:
    tuple val(comp_id), path("exomePeak2_output/Mod.bed"), emit: bedPeaks
    
	script:
	def filecontent = "SampleID\tBamIP\tBamInput"
    def unzip_anno = unzipCmd(annofile)
    def cmd_anno = unzip_anno[1]
    def anno_name = unzip_anno[0]
    def cmd_clean = unzip_anno[2]
    def extrapar = ""	
    if (params.EXTRAPARS != "") {
    	extrapar = ", ${params.EXTRAPARS}"
    }
    """
	${cmd_anno}
	
	cat > CMD.R << EOL

library('exomePeak2')
set.seed(1)

GENE_ANNO_GTF = "\$PWD/${anno_name}"

library("GenomicFeatures")
txdb_obj <- makeTxDbFromGFF(GENE_ANNO_GTF, format= "gtf")


f1 = "${sample}"
IP_BAM = c(f1)

f2 = "${input}"
INPUT_BAM = c(f2)

exomePeak2(bam_ip = IP_BAM,
           bam_input = INPUT_BAM,
           txdb = txdb_obj ${extrapar})

EOL

	Rscript  CMD.R 
	${cmd_clean}
    """
}



workflow CALL {
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
