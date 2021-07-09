/*
*  Bedtools 
*/

params.LABEL = ""
params.EXTRAPARS = ""

params.OUTPUT = ""
params.CONTAINER = "quay.io/biocontainers/bioconductor-tximeta:1.10.0--r41hdfd78af_0"

process getVersion {
    container params.CONTAINER

    output:
	stdout emit: out    
    
    shell:
    """
    echo tximeta' '`Rscript -e "library('tximeta'); packageVersion('tximeta')"` 2>/dev/null
    """
}

process processSalmon {
    label (params.LABEL)
    tag { "all" }
    container params.CONTAINER
    if (params.OUTPUT != "") { publishDir(params.OUTPUT, mode:'copy') }

    input:
    path(input)
    path(index)
    path(transcriptome)
    path(annotation)
    val(anno_type)
    val(org_name)
    val(release)

    output:
	path("gene_counts.txt") , emit: gcounts
	path("transcript_counts.txt"), emit: tcounts
	path("gse.rds"), emit: gse
    
	script:
def Rscript = """
# R --slave --args index fasta gtf gtf_t organism genome release < import.R

args<-commandArgs(TRUE)
index_i<-args[1]
fasta_i<-args[2]
gtf_i<-args[3]
gtf_t<-args[4]
organism_n<-args[5]
genome_n<-args[6]
release<-args[7]

suppressPackageStartupMessages(library(tximeta))
dir<-getwd()
files <- dir(dir, recursive=TRUE, pattern="quant.sf", full.names=TRUE)
names<-basename(dirname(files))
coldata <- data.frame(files, names=names, stringsAsFactors=FALSE)

#se <- tximeta(coldata)
dir.create(file.path(dir, "tmp"), showWarnings = FALSE)

setTximetaBFC(file.path(dir, "tmp"), quiet = FALSE)

indexDir <- file.path(index_i)
fasta <- file.path(fasta_i)
gtfPath <- file.path(gtf_i)

makeLinkedTxome(indexDir=indexDir,
                source=gtf_t,
                organism=organism_n,
                genome=genome_n,
                fasta=fasta,
                release=release,
                gtf=gtfPath,
                write=FALSE)
 
se <- tximeta(coldata)
gse <- summarizeToGene(se)
saveRDS(gse, file = "gse.rds")

library("SummarizedExperiment")
gcounts<-assay(gse)
tcounts<-assay(se)

write.table(gcounts, "gene_counts.txt", col.names = TRUE, sep="\t")
write.table(tcounts, "transcript_counts.txt", col.names = TRUE, sep="\t")
"""

    """
	cat > import.R << 'EOF'
${Rscript}
EOF
	
	R --slave --args \$PWD/${index} \$PWD/${transcriptome} \$PWD/${annotation} ${anno_type} ${org_name} ${release} < import.R
    """
}

workflow PROCESS_SALMON {
    take: 
    input
    index
    transcriptome
    annotation
    anno_type
    org_name
    release
    
    main:
		processSalmon(input, index, transcriptome, annotation, anno_type, org_name, release)
    emit:
    	gcounts = processSalmon.out.gcounts
    	tcounts = processSalmon.out.tcounts
    	gse = processSalmon.out.gse
}





workflow GET_VERSION {
    main:
		getVersion()
    emit:
    	getVersion.out
}    
