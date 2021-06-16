install.packages(c('usethis', 'covr', 'httr', 'roxygen2', 'rversions', 'devtools'), dependencies=TRUE, repos='http://cran.us.r-project.org')

install.packages("curl")
install.packages("openssl")

install.packages("httr")
install.packages("xml2")
install.packages("lintr")
install.packages("usethis")

install.packages("roxygen2")
install.packages("rversions")
install.packages("devtools")

install.packages("BiocManager")

BiocManager::install("SummarizedExperiment", update = TRUE, ask = FALSE)
BiocManager::install("GenomeInfoDb", update = TRUE, ask = FALSE)
BiocManager::install("mosaics", update = TRUE, ask = FALSE)
BiocManager::install("Rsubread", update = TRUE, ask = FALSE)
BiocManager::install("rtracklayer", update = TRUE, ask = FALSE)
BiocManager::install("GenomicFeatures", update = TRUE, ask = FALSE)
BiocManager::install("GenomicAlignments", update = TRUE, ask = FALSE)
BiocManager::install("GenomicRanges", update = TRUE, ask = FALSE)
  
library(devtools)
install_github("ADelgadoT/MoAIMS/codes")

install.packages("eulerr")
install.packages("data.table")
