library(devtools)
INSTALL_opts=c("--no-help", "--no-html")
library(devtools)
devtools::install_github("ZW-xjtlu/exomePeak")
install.packages("https://cran.r-project.org/src/contrib/Archive/QNB/QNB_1.1.11.tar.gz", repos = NULL, type="source")
BiocManager::install("Rsubread")
BiocManager::install("DESeq2", update = TRUE, ask = FALSE)

devtools::install_github("al-mcintyre/DEQ")
