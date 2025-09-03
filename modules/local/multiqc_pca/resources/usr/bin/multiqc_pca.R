#!/usr/local/bin/Rscript --vanilla

source(Sys.which("functions.R"))



# Load common pars
parser <- getCommonPars()

#Define desired outputs:
#SPECIFIC FEATURES:
parser$add_argument("-input", "--Input_f", required="true", type="character", help="Path or matrix file to read counts")
parser$add_argument("-type", "--Input_type", type="character", required="true", help="Type of input: counts, salmon, star, matrix")
parser$add_argument("-label", "--Add_labels", type="character", default="", help="Column to be used for labels")
parser$add_argument("-dtrans", "--Desc_Tx", type="character", default="tx2gene.csv", help="Path to read tx2gene file (ONLY FOR SALMON) [default %(default)s]")
parser$add_argument("-strand", "--Strand_id", type="integer", default=4, help="Strand to analyze, 2: unstranded, 3: forward, 4 reverse (ONLY FOR STAR) [default %(default)s]")
parser$add_argument("-batch", "--add_batch", action='store_true', help="Remove batch effect using the column batch in desc file")
parser$add_argument("-pcnum", "--Number_principal_components", type="integer", default=4, help="Number of principal components (one for each column) data to extract to the data table, [default %(default)s]")

#Get command line options, if help option encountered - print help and exit:
args <- parser$parse_args()

if (args$Input_type == "counts") {
	out <- makeDDSFromCounts(args$Input_f, args$Desc_exp, args$Assay_field)
} else if (args$Input_type == "salmon") {
	out <- makeDDSFromSalmon(args$Input_f, args$Desc_exp, args$Assay_field, args$Desc_Tx)
} else if (args$Input_type == "star") {
	out <- makeDDSFromStar(args$Input_f, args$Desc_exp, args$Assay_field, args$Strand_id)
} else if (args$Input_type == "matrix") {
	out <- makeDDSFromMatrix(args$Input_f, args$Desc_exp, args$Assay_field)
} else {
	print("please define one correct input type!")
    stop()
}

dds <- filterDDS(out$dds, args$min_count)

vsd<-makeVST(dds, FALSE)

if (args$add_batch) {
    if (grepl("batch", args$Assay_field)) {
        stop("ERROR DO NOT PUT THE BATCH AS CONTROLLING FACTOR!!!")
        
    }
    mat<-assay(vsd)
    mm <- model.matrix(as.formula(args$Assay_field), colData(vsd))
    mat <- limma::removeBatchEffect(mat, batch=vsd$batch, design=mm)
    assay(vsd) <- mat
}

# Removing ~ for the condition field
condition<-gsub("~","",args$Assay_field) 
print(args$Number_principal_components)

p<-plotPCA_plus(vsd, intgroup=condition, returnData = TRUE,pcnum=args$Number_principal_components)

# Defining PCA data and variance from the list 
pca_data<-p[["data"]]
pca_variance<-p[["variance"]]

# Generate distinct colors dynamically for the factor levels of the condition column 
unique_vals<- unique(pca_data[[condition]])
num_unique <- length(unique_vals)

# Generate distinct colors dynamically
color_palette <- c("#f5d905","#db14ce","#3d98bf","#2114db","#db1432","#14db1e","#db9c14","#3dbfac","#5c5952","#4f0e29") # Creates distinct colors

# Create a lookup table
color_map <- setNames(color_palette[1:num_unique], unique_vals)

# Add a new column 'color' to pca_data
pca_data$color <- color_map[pca_data[[condition]]]

# Save PCA data and variance tables
write.table(pca_data, file="PCA_data.tsv", sep="\t", quote=FALSE, row.names=FALSE)
write.table(pca_variance, file="PCA_variance.tsv", sep="\t", quote=FALSE, row.names=FALSE)

if (args$Desc_genes != "") {
    print("found desc gene file")

    # Get results
    desc <- makeDesc(args$Desc_genes)
    printCounts(dds, desc)


}
