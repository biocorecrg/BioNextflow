suppressMessages(library("argparse"))
suppressMessages(library("stringr"))
suppressMessages(library("DESeq2"))
suppressMessages(library("EnhancedVolcano"))
suppressMessages(library("tximport"))
suppressMessages(library("pheatmap"))


##Argument parser:
#Create parser object

# functions
#Get common pars
getCommonPars <- function(desc_exp) {
	parser <- ArgumentParser()

	#Define desired outputs:
	#GLOBAL FEATURES:
	parser$add_argument("-desc", "--Desc_exp", default="desc.txt", type="character", help="File describing experiments, [default %(default)s]")
	parser$add_argument("-dgenes", "--Desc_genes", default="gene_desc.txt", type="character", help="Files describing genes, [default %(default)s]")
    parser$add_argument("-assay", "--Assay_field", required="true", type="character", help="Formula to be used for DE analysis, example \"~ time\"")
	parser$add_argument("-min", "--min_count", type="integer", default=1, help="Number of minimum read per row [default %(default)s]")
	return (parser)
}


makeDDSFromCounts <- function(input, desc, field) {
	# Read counts and make count table
	ff <- list.files( path = input, recursive=F, pattern = "*.counts$", full.names = TRUE )
	counts.files <- lapply( ff, read.table)
	counts <- as.data.frame( sapply( counts.files, function(x) x[ , 2 ] ) )
	fn <- basename(ff)
	fn <- gsub( ".counts", "", fn)
	colnames(counts) <- fn
	row.names(counts) <- counts.files[[1]]$V1

	#Make coldata for DESEq2
	coldata <- makeColData(desc, fn)

	#DESEQ2
	dds <- DESeqDataSetFromMatrix(countData = counts,
                    colData = coldata,
                    design = as.formula(field))

	return(list("dds" = dds, "coldata" = coldata))
}

makeDDSFromSalmon <- function(input, desc, field, desctx) {
	# List the quantification files from Salmon: one quant.sf file per sample
	files <- dir(input, recursive=TRUE, pattern="quant.sf", full.names=TRUE)
	names(files) <-basename(dirname(files))

	tx2gene <- read.table(desctx, 
		sep="\t",
		header=F)

	txi <- tximport(files,
		type = "salmon", 
		tx2gene = tx2gene)

	#Make coldata for DESEq2
	coldata <- makeColData(desc, names(files))

	#DESEQ2
	dds <- DESeqDataSetFromTximport(txi,
                        colData = coldata,
                        design = as.formula(field))
	return(list("dds" = dds, "coldata" = coldata))
}

makeDDSFromStar <- function(input, desc, field, number) {

	# Read counts and make count table
	ff <- list.files( path = input, recursive=F, pattern = "*ReadsPerGene.out.tab$", full.names = TRUE )
	counts.files <- lapply( ff, read.table, skip = 4 )
	counts <- as.data.frame( sapply( counts.files, function(x) x[ , number ] ) )
	fn <- basename(ff)
	fn <- gsub( "ReadsPerGene.out.tab", "", fn)
	colnames(counts) <- fn
	row.names(counts) <- counts.files[[1]]$V1

	#Make coldata for DESEq2
	coldata <- makeColData(desc, fn)

	#DESEQ2
	dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = coldata,
                              design = as.formula(field))
	return(list("dds" = dds, "coldata" = coldata))
}

makeDDSFromMatrix <- function(input, desc, field) {

	if (basename(input) == "annotated_counts.txt.gz") {
 	 	counts <- read.delim(gzfile(input), row.names = 4, check.names = FALSE)
	} else {
 	 	counts <- read.csv(input, row.names = 1, check.names = FALSE)
    }
    
	if ('gene.name' %in% names(counts)) {
		counts$gene.name<-NULL
		counts$gene.type<-NULL
	}
	
	else if ('Strand' %in% names(counts)) {
		counts$Chr<-NULL
		counts$Start<-NULL
		counts$End<-NULL
		counts$Score<-NULL
		counts$Strand<-NULL
		counts$gene.id<-NULL
		counts$gene.name<-NULL
		counts$gene.type<-NULL
	}


    
	#Make coldata for DESEq2
	coldata <- makeColData(desc, names(counts))

	#DESEQ2
	dds <- DESeqDataSetFromMatrix(countData = counts,
                        colData = coldata,
                        design = as.formula(field))

	return(list("dds" = dds, "coldata" = coldata))
}


#Make coldata for DESEq2
makeColData <- function(desc_exp, fn) {
	desc_data<-read.table(desc_exp, sep="\t", header=T)
	coldata<-desc_data[match(fn, desc_data$file),]
	row.names(coldata)<-fn
	if ("group" %in% names(coldata)) {
		stop("****** 'group' cannot be in the names of the description file, please change it ******")
	}
	all(rownames(coldata) == colnames(counts))
	return(coldata)
}

#READ gene desc file
makeDesc <- function(desc_gene) {
	desc<-read.table(file=desc_gene, sep="\t", header=F)
	names(desc)<-c("gene.id", "gene.name", "gene.type")
	desc <- desc[-grep("gene_id", desc$gene.id), ]
	return(desc)
}

filterDDS <- function(dds, mincount) {
	dds <- dds[ rowSums(counts(dds)) > mincount, ]
	dds <- DESeq(dds)
	return(dds)
}

# functions
ddsResult <- function(dds, desc, assay_field_raw, cond1, cond2 ) {
    assay_field <- str_remove(assay_field_raw, '~ ')
	filename<-paste(assay_field, cond1, cond2, sep="_")
	filename<-paste(filename, ".csv", sep="")
	res<-results(dds, contrast=c(assay_field,cond1,cond2))
	resOrdered <- res[order(res$padj),]
	resOrdered$ids<-row.names(resOrdered)
	resMerged<-merge(as.data.frame(resOrdered),desc, all = TRUE, by.x="ids", by.y="gene.id", sort=FALSE)
	resMerged$ID<-NULL
	write.csv(resMerged, file=filename, row.names = FALSE)          
    return (resMerged)
}

makeVolcano <- function(results, assay_field_raw, treat, ctrl, pcut, l2fc, h, w) {
    assayfield <- str_remove(assay_field_raw, '~ ')

	res_for_volc<-results[, c(8,3,7)]
	volcano_filename<-paste(assayfield, treat, ctrl, sep="_")
	volcano_file<-paste(volcano_filename, "_volc.pdf", sep="")

	subtitle<-paste(treat, ctrl, sep=" vs ")

	myplot<-EnhancedVolcano(res_for_volc,
    title = "Volcano plot",
    subtitle = subtitle,
    lab = res_for_volc$gene.name,
    x = 'log2FoldChange',
    pCutoff = pcut,
    FCcutoff = l2fc,
    y = 'padj')

	pdf(volcano_file, width = w, height = h)
    print(myplot)
	dev.off()
}


makeVST <- function(dds, blindval = FALSE) {
    vsd <- tryCatch( 
    {
      vst(dds, blind=blindval)
    },
      error = function(e) {
          varianceStabilizingTransformation(dds, blind=blindval)
      }
   )
   return(vsd)
}

makeHeatmapDE <- function(results, dds, assay_field_raw, treat, ctrl, pcut, l2fc, h, w, maxnum=100000, sample_list=NULL, genes_file=NULL) {

    assayfield <- str_remove(assay_field_raw, '~ ')
    sig_res <- subset(results, padj <= pcut)
    if (maxnum < 100000) {
        sig_res<-sig_res[1:maxnum, ]
    }

    sig_ids <- sig_res$ids 

    # varianceStabilizingTransformation or vst?
    vsd<-assay(makeVST(dds, FALSE))
    sig_vsd_ids <- vsd[sig_ids, ]

    conv<-data.frame("ids"=results["ids"], "gene.name"=results["gene.name"])
    sig_vsd<-merge(as.data.frame(sig_vsd_ids),conv, all = FALSE, by.y="ids", by.x="row.names", sort=FALSE)
    sig_vsd["Row.names"]<-NULL
    row.names(sig_vsd)<-deduplicateIDs(sig_vsd$gene.name)
    sig_vsd["gene.name"]<-NULL

    if (!is.null(sample_list)) {
        samples<-str_split(sample_list, ",", simplify=TRUE)
        sig_vsd <- sig_vsd[ , samples]
    }

    if (!is.null(genes_file)) {
        mygenes<-trimws(readLines(genes_file))
	    mygenes <- mygenes[mygenes != ""]		
        sig_vsd <- sig_vsd[intersect(mygenes, rownames(sig_vsd)), ]
	    sig_size <- dim(sig_vsd)
	    print(paste0("we will analyze only ", sig_size[1], " genes"))
    } 

    heatmap_filename<-paste(assayfield, treat, ctrl, sep="_")
    heatmap_file<-paste(heatmap_filename, "_heat.pdf", sep="")

    subtitle<-paste("Heatmap of", treat, "vs",  ctrl, "( padj <=", pcut, ")",  sep=" ")
    pheatmap(sig_vsd, main=subtitle, cluster_rows=TRUE, filename = heatmap_file, show_rownames=TRUE, cluster_cols=TRUE, height=h, width=w)
}



deduplicateIDs <- function(vec) {
    dedup = ave(vec, vec, FUN =  \(x) if(length(x) > 1) paste(x, seq_along(x), sep = "__") else x)
    return(dedup)
}


plotPCA_plus <- function (object, intgroup = "condition", ntop = 500, returnData = FALSE, pcnum=6) {
    rv <- rowVars(assay(object), useNames=TRUE)
    select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
    pca <- prcomp(t(assay(object)[select, ]))
    percentVar <- pca$sdev^2/sum(pca$sdev^2)
    if (!all(intgroup %in% names(colData(object)))) {
        stop("the argument 'intgroup' should specify columns of colData(dds)")
    }
    intgroup.df <- as.data.frame(colData(object)[, intgroup,  drop = FALSE])
    group <- if (length(intgroup) > 1) {
        factor(apply(intgroup.df, 1, paste, collapse = " : "))
    }
    else { colData(object)[[intgroup]] }
    d <- data.frame(group = group, intgroup.df, name = colnames(object))
	for (i in 1:pcnum) {
	d[[paste0("PC", i)]] <- pca$x[, i]
	}
	attr(d, "percentVar") <- percentVar[1:pcnum]
	pv <- round(100*attr(d, "percentVar"))
	PC_names <- c(paste("PC", 1:pcnum, sep=""))
	pv_df<-data.frame(PC_names, pv)

    if (returnData) {
 	    rots<-pca$rotation
 	    mylist<-list("data"=d,"rotation"=rots,"variance"=pv_df)
    	return(mylist)
	}
	

}

makePCA <- function(vsd, pca_fields, labels="", all_fields, PCA=1, PCB=2) {
    
    filename<-paste("PCA", pca_fields, sep="_")    
    filename<-paste(filename, "pdf", sep=".")    

	pdf(filename)
    PCA_lab<-paste("PC", PCA, sep="")
    PCB_lab<-paste("PC", PCB, sep="")
	cond=sym(pca_fields)
	PCA_labs=sym(PCA_lab)
	PCB_labs=sym(PCB_lab)
    pcaData_plus <- plotPCA_plus(vsd, intgroup = all_fields, returnData=TRUE)
    pcaData <- pcaData_plus[["data"]]
    rot.df<-as.data.frame(pcaData_plus[["rotation"]])
    rot.sel<-rot.df[, c(PCA_lab, PCB_lab)]
	percentVar <- round(100 * attr(pcaData, "percentVar"))
	myplot <- ggplot(pcaData, aes(!!PCA_labs, !!PCB_labs, color=!!cond))
	if (labels!="") {
		lab=sym(labels)
	    myplot <- ggplot(pcaData, aes(!!PCA_labs, !!PCB_labs, color=!!cond, label=!!lab)) +
	    geom_text(alpha = 0.8, size = 4, show.legend = FALSE)
	}
	myplot <- myplot + geom_point(size=2, alpha = 0.5) +
	xlab(paste0(PCA_lab, ": ",percentVar[PCA],"% variance")) +
	ylab(paste0(PCB_lab, ": ",percentVar[PCB],"% variance")) +
	scale_shape_identity() +
	scale_x_reverse() + theme_classic()
	print(myplot)
	dev.off()
	return(rot.sel)
}

printCounts <- function(dds, desc, ofile=NULL) {
    if (!is.null(ofile)) {
       vst_genes_f<-paste(ofile, "vst.genes", sep="_")
       norm_genes_f<-paste(ofile, "norm_counts.genes", sep="_")
       raw_genes_f<-paste(ofile, "raw_counts.genes", sep="_")
    } else {
       vst_genes_f<-"vst.genes"
       norm_genes_f<-"norm_counts.genes"
       raw_genes_f<-"raw_counts.genes"
    }

    vsd = makeVST(dds) 

    resMerged<-merge(as.data.frame(assay(vsd)), desc, all.x=F, by.x=0, by.y="gene.id", sort=FALSE)
    resMerged$ID<-NULL
    write.csv(resMerged, file=vst_genes_f, row.names = FALSE)          
    
    norcounts<-counts(dds, normalized=TRUE)
    rawcounts<-counts(dds, normalized=FALSE)
    resMerged2<-merge(as.data.frame(norcounts) ,desc, all.x=F, by.x=0, by.y="gene.id", sort=FALSE)

    resMerged2$ID<-NULL
    write.csv(resMerged2, file=norm_genes_f, row.names = FALSE)

    resMerged3<-merge(as.data.frame(rawcounts) ,desc, all.x=F, by.x=0, by.y="gene.id", sort=FALSE)
    resMerged3$ID<-NULL
 
    write.csv(resMerged3, file=raw_genes_f, row.names = FALSE)

}

printContrGenes <- function(rot_sel, desc) {

	resMerged<-merge(rot_sel, desc, all.x=F, by.x=0, by.y="gene.id", sort=FALSE)
	resOrdered <- resMerged[order(resMerged[, 2], decreasing = TRUE), ]
	write.csv(resOrdered, file="pca.genes", row.names = FALSE)          

}
