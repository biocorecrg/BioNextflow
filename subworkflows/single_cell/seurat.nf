/*
*  Seurat module
*/

params.LABEL = ""
params.CONTAINER = "biocorecrg/seurat:4.1"
params.OUTPUT = ""


process preproc_alevin {
    tag "${id}"
    label (params.LABEL)
    container params.CONTAINER
    if (params.OUTPUT != "") { publishDir(params.OUTPUT, mode:'copy') }

    input:
    tuple val(id), path(quants_folder)
    path(tx2gene)

    output:
    tuple val(id), path("${id}.rds"),  path("${id}_*.pdf")

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

mt_genes<-subset(tx2gene2, grepl("^MT-", tx2gene2\$ens_name))
rb_genes<-subset(tx2gene2, grepl("^RP[SL]", tx2gene2\$ens_name))
files<- file.path("${quants_folder}/alevin/quants_mat.gz")
txi <- tximport(files, type="alevin")

ens_gene_ids2<-str_split(row.names(txi\$counts) , "\\\\.", n = Inf, simplify = TRUE) 

row.names(txi\$counts)<-ens_gene_ids2[, 1]

seurObj <- CreateSeuratObject(counts = txi\$counts , min.cells = 3, min.features = 200, project = "${id}")

counts <- seurObj@assays\$RNA@counts
mt <- mt_genes\$ensg_id[mt_genes\$ensg_id %in% rownames(counts)]
rb <- rb_genes\$ensg_id[rb_genes\$ensg_id %in% rownames(counts)]

seurObj[["percent.mt"]] <- PercentageFeatureSet(seurObj, features = mt,  assay = 'RNA')
seurObj[["percent.rb"]] <- PercentageFeatureSet(seurObj, features = rb,  assay = 'RNA')

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
top10g <- tx2gene2[match(top10, tx2gene2\$ensg_id),]

pdf(paste("${id}", "_vf.pdf", sep=""), width=10)
plot1 <- VariableFeaturePlot(seurObj)
plot2 <- LabelPoints(plot = plot1, points = top10, labels=top10g\$ens_name, repel = TRUE)
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
    
top10_cls_g <- tx2gene2[match(top10\$gene, tx2gene2\$ensg_id),]
top10_cls<-top10_cls_g\$ens_name
names(top10_cls)<-top10_cls_g\$ensg_id

pdf(paste("${id}", "_hm.pdf", sep=""), width=10, height=20)
DoHeatmap(seurObj, features = top10\$gene) + NoLegend() + scale_y_discrete(labels = top10_cls)
dev.off()

# We use HumanPrimaryCell as default annotation

ref.data <- HumanPrimaryCellAtlasData(ensembl=TRUE)
sce <- as.SingleCellExperiment(DietSeurat(seurObj))

predictions.main <- SingleR(test=sce, assay.type.test=1, 
    ref=ref.data, labels=ref.data\$label.main)

table(predictions.main\$labels)

seurObj@meta.data\$predictions.main <- predictions.main\$labels

seurObj <- SetIdent(seurObj, value = "predictions.main")

pdf(paste("${id}", "_ann.pdf", sep=""), width=10, height=10)
DimPlot(seurObj, label = T , repel = T, label.size = 3) + NoLegend()
dev.off()




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
    
    main:
    if(type == "alevin") {
	    out = preproc_alevin(input, tx2gene)
    } else {
		 error "ERROR!!! Specify a method for import\n" 
    }
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
 
