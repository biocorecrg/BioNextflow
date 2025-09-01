/*
*  Seurat module
*/

params.LABEL = ""
params.CONTAINER = 'biocorecrg/sc_benchmark:0.2'

params.OUTPUT = ""


process preproc_for_cellranger {
    tag "${id}"
    label (params.LABEL)
    container params.CONTAINER


    input:
    tuple val(id), path(quants_folder)

    output:
    tuple val(id), path("${id}.rds")

    script:
	"""
cat > CMD.R << 'EOL'

library(dplyr)
library(Seurat)
library(patchwork)
library(tximport)
library("ggplot2")
library("SingleR")
library("celldex")
library("stringr")
 
data_dir <- "${quants_folder}/outs/filtered_feature_bc_matrix"
data <- Read10X(data.dir = data_dir)


if(is.list(data)) {
	olddata <- data
	data <- olddata[["Gene Expression"]]
}

seurObj = CreateSeuratObject(counts = data)
# save seurat object
saveRDS(seurObj, file = "${id}.rds")
quit("no")

EOL

	Rscript CMD.R 

	"""


}

process preproc_for_cellranger_multi {
    tag "${id}"
    label (params.LABEL)
    container params.CONTAINER


    input:
    tuple val(id), path(quants_folder)

    output:
    tuple val(id), path("${id}.rds")

    script:
	"""
cat > CMD.R << 'EOL'

library(dplyr)
library(Seurat)
library(patchwork)
library(tximport)
library("ggplot2")
library("SingleR")
library("celldex")
library("stringr")
 
data_dir <- "${quants_folder}/count/sample_filtered_feature_bc_matrix"
data <- Read10X(data.dir = data_dir)


if(is.list(data)) {
	olddata <- data
	data <- olddata[["Gene Expression"]]
}

seurObj = CreateSeuratObject(counts = data)
# save seurat object
saveRDS(seurObj, file = "${id}.rds")
quit("no")

EOL

	Rscript CMD.R 

	"""


}

process preproc_for_alevin {

    tag "${id}"
    label (params.LABEL)
    container params.CONTAINER

    input:
    tuple val(id), path(quants_folder)
    path(tx2gene)

    output:
    tuple val(id), path("${id}.rds")

    script:
	"""
cat > CMD.R << 'EOL'

library(dplyr)
library(Seurat)
library(patchwork)
library(tximport)
library("ggplot2")
library("SingleR")
library("celldex")
library("stringr")
 
tx2gene<-read.delim("${tx2gene}", header=FALSE)
colnames(tx2gene)<-c("enst_id", "ensg_id", "ens_name")

ens_gene_ids<-str_split(tx2gene\$ensg_id , "\\\\.", n = Inf, simplify = TRUE) 

tx2gene\$ensg_id <- ens_gene_ids[, 1]
tx2gene2<-unique(tx2gene[, c(2,3)])

files<- file.path("${quants_folder}/alevin/quants_mat.gz")
txi <- tximport(files, type="alevin")

ens_gene_ids2<-str_split(row.names(txi\$counts) , "\\\\.", n = Inf, simplify = TRUE) 

row.names(txi\$counts)<-ens_gene_ids2[, 1]

seurObj <- CreateSeuratObject(counts = txi\$counts , min.cells = 3, min.features = 200, project = "${id}")


# save seurat object
saveRDS(seurObj, file = "${id}.rds")
quit("no")

EOL

	Rscript CMD.R 

	"""


}

process preproc_for_parse {

    tag "${id}"
    label (params.LABEL)
    container params.CONTAINER
    

    input:
    tuple val(id), path(quants_folder)

    output:
    tuple val(id), path("${id}.rds")

    script:
	"""
cat > CMD.R << 'EOL'

library(dplyr)
library(Seurat)
library(patchwork)
library(tximport)
library("ggplot2")
library("SingleR")
library("celldex")
library("stringr")

expression_matrix <- ReadParseBio("${quants_folder}/DGE_filtered")
gene_genomes <- read.csv(paste0("${quants_folder}/DGE_filtered", "/all_genes.csv"))

# For multi species... 
#comb_gene<-paste(gene_genomes\$genome, rownames(expression_matrix), sep="-")

# if empty gene names are present, name them unknown.
rownames(expression_matrix)[rownames(expression_matrix) == ""] <- "unknown"
# Read in cell meta data
cell_meta <- read.csv(paste0("${quants_folder}/DGE_filtered", "/cell_metadata.csv"), row.names = 1)
seurObj <- CreateSeuratObject(expression_matrix, names.field = 0, meta.data = cell_meta, min.cells = 3, min.features = 200, project = "${id}")
Idents(seurObj) <- seurObj@meta.data\$orig.ident


# save seurat object
saveRDS(seurObj, file = "${id}.rds")
quit("no")

EOL

	Rscript CMD.R 

	"""


}


process seurat {
    tag "${id}"
    label (params.LABEL)
    container params.CONTAINER

    if (params.OUTPUT != "") { publishDir(params.OUTPUT, mode:'copy') }

    input:
    tuple val(id), path("input_ori.rds")
    val(genome)

    output:
    tuple val(id), path("${id}.rds"),  path("${id}_*.pdf")

    script:
    def globalsMax = (task.memory.toBytes() / 2) as long

    def script_anno = ""
    
    if (genome == "human" || genome == "mouse") {
    	if (genome == "human") {
			script_anno = "ref.data <- HumanPrimaryCellAtlasData()"
    	} else if (genome == "mouse") {
        	script_anno = "ref.data <- MouseRNAseqData()"
    	}
    	script_anno = script_anno + """    
sce <- as.SingleCellExperiment(DietSeurat(seurObj))

predictions.main <- SingleR(test=sce, assay.type.test=1, 
    ref=ref.data, labels=ref.data\$label.main)

table(predictions.main\$labels)

seurObj@meta.data\$predictions.main <- predictions.main\$labels

seurObj <- SetIdent(seurObj, value = "predictions.main")

pdf(paste("${id}", "_ann.pdf", sep=""), width=10, height=10)
DimPlot(seurObj, label = T , repel = T, label.size = 3) + NoLegend()
dev.off()
"""
	} 	


  
	"""
cat > CMD.R << 'EOL'

library(dplyr)
library(Seurat)
library(patchwork)
library(tximport)
library("ggplot2")
library("SingleR")
library("celldex")
library("stringr")
 
options(future.globals.maxSize = ${globalsMax} ) 

seurObj <-readRDS("input_ori.rds")

# Extract mitochondrial genes (case-insensitive)
mt_genes <- grep("^mt-", rownames(seurObj[["RNA"]]), value = TRUE, ignore.case = TRUE)

# Check if any mitochondrial genes were found
if (length(mt_genes) == 0) {
  stop("No mitochondrial genes (starting with 'mt-' or 'MT-'') found in the dataset.")
}

# Compute percent mitochondrial content
seurObj[["percent.mt"]] <- PercentageFeatureSet(seurObj, features = mt_genes, assay = 'RNA')
seurObj[["percent.rb"]] <- PercentageFeatureSet(seurObj, pattern = "^RP[SL]",  assay = 'RNA')

pdf(paste("${id}", "_vp.pdf", sep=""), width=10, height=20)
VlnPlot(seurObj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb") , ncol = 4)
dev.off()

# subsetting and normalize 
# We subset cells with less than 0.05 percentile MT and between 200 and 0.99 percentile of nFeatures 
cutoff.mt<-round(quantile(seurObj[["percent.mt"]][, 1], c(.95)))[[1]]
cutoff.nFeat<-round(quantile(seurObj[["nFeature_RNA"]][, 1], c(.99)))[[1]]


pdf(paste("${id}", "_fc.pdf", sep=""), width=10)
plot1 <- FeatureScatter(seurObj, feature1 = "nCount_RNA", feature2 = "percent.mt") + theme(legend.position="none") + geom_hline(yintercept = cutoff.mt) + annotate("text", x=-100, y=cutoff.mt+2, label= cutoff.mt)
plot2 <- FeatureScatter(seurObj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + theme(legend.position="none") + geom_hline(yintercept = 200) + geom_hline(yintercept = cutoff.nFeat) + annotate("text", x=c(-400,-400), y=c(300, cutoff.nFeat+100), label= c(200, cutoff.nFeat))
plot1 + plot2
dev.off()


seurObj <- subset(seurObj, subset = nFeature_RNA > 200 & nFeature_RNA < cutoff.nFeat[1] & percent.mt < cutoff.mt)

seurObj <- SCTransform(seurObj)
seurObj <- FindVariableFeatures(seurObj, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(seurObj), 10)

pdf(paste("${id}", "_vf.pdf", sep=""), width=10)
plot1 <- VariableFeaturePlot(seurObj)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
dev.off()

# PCA
seurObj <- RunPCA(seurObj, features = VariableFeatures(object = seurObj))

pdf(paste("${id}", "_ep.pdf", sep=""), width=10)
ElbowPlot(seurObj)
dev.off()

# Here we use 1:15 to PCs. So it might be wise to have a look at elbowplot 
# before trusting these results!

seurObj <- FindNeighbors(seurObj, dims = 1:15)
seurObj <- FindClusters(seurObj, resolution = 0.5)

seurObj <- RunUMAP(seurObj, dims = 1:15)

pdf(paste("${id}", "_dp.pdf", sep=""), width=10)
DimPlot(seurObj, reduction = "umap")
dev.off()

# find all markers of cluster 2
cluster2.markers <- FindMarkers(seurObj, ident.1 = 2, min.pct = 0.25)
head(cluster2.markers, n = 5)

# find markers for every cluster compared to all remaining cells, report only the positive
# ones

seurObj.markers <- FindAllMarkers(seurObj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
seurObj.markers %>%
    group_by(cluster) %>%
    slice_max(n = 2, order_by = avg_log2FC)

cluster0.markers <- FindMarkers(seurObj, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)

seurObj.markers %>%
    group_by(cluster) %>%
    top_n(n = 10, wt = avg_log2FC) -> top10
    

pdf(paste("${id}", "_hm.pdf", sep=""), width=10, height=20)
DoHeatmap(seurObj, features = top10\$gene) + NoLegend()
dev.off()

${script_anno}

# save seurat object
saveRDS(seurObj, file = "${id}.rds")
quit("no")

EOL

	Rscript CMD.R 

	"""

}


workflow PREPROCESS {
    take: 
    input
    type
    tx2gene
    genome
    
    main:
    if(type == "alevin") {
	    rds = preproc_for_alevin(input, tx2gene)
    } else if (type == "cellranger") {
    	rds = preproc_for_cellranger(input)
    } else if (type == "parse") {
		rds =  preproc_for_parse(input)
    }
    else if (type == "cellranger_multi") {
    	rds = preproc_for_cellranger_multi(input)
    }
    else {
		 error "ERROR!!! Specify a method for import\n" 
    }
    out = seurat(rds, genome)
    
    emit:
    out

}

process getVersion {
    container params.CONTAINER

    output:
	stdout emit: out    
    
    shell:
    """
    fastqc --version
    """
}




workflow GET_VERSION {
    main:
		getVersion()
    emit:
    	getVersion.out
}   
 
