options(repos = c(CRAN = "https://cloud.r-project.org"))
.libPaths(c("~/R/library", .libPaths()))
cat("Installing to:", .libPaths()[1], "\n")
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", lib = "~/R/library")
}
library(BiocManager)
BiocManager::install("QDNAseq", update = FALSE, ask = FALSE, lib = "~/R/library")
BiocManager::install("QDNAseq.hg38", update = FALSE, ask = FALSE, lib = "~/R/library")
cat("Verifying...\n")
library(QDNAseq)
cat("QDNAseq:", as.character(packageVersion("QDNAseq")), "\n")
