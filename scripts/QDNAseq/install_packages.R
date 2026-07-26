# One-off environment bootstrap: installs the R packages the analysis scripts need into a
# user-writable library on ALICE, since the R/4.3.1 module's system library is read-only.
#
# Run once per account before any pipeline script. Safe to re-run; already-installed
# packages are skipped.

options(repos = c(CRAN = "https://cloud.r-project.org"))

# The library directory must exist before it can be used. .libPaths() silently discards
# entries that do not exist, so on a clean account this path was dropped, the install
# below went to the system library instead (or failed for lack of permission), and the
# "Installing to:" line reported the wrong destination. Create it first.
LIB <- path.expand("~/R/library")
if (!dir.exists(LIB)) dir.create(LIB, recursive = TRUE)

.libPaths(c(LIB, .libPaths()))
cat("Installing to:", .libPaths()[1], "\n")

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", lib = LIB)
}
library(BiocManager)

# Bioconductor packages. QDNAseq and Biobase drive the copy-number pipeline; biomaRt is
# used for CNA-to-gene mapping; GenVisR draws the cohort-level CNA figures.
bioc_pkgs <- c("QDNAseq", "Biobase", "biomaRt", "GenVisR")

# CRAN packages used across the analysis and plotting scripts.
cran_pkgs <- c("dplyr", "tidyr", "readr", "readxl", "ggplot2", "ggrepel", "patchwork")

for (p in c(bioc_pkgs, cran_pkgs)) {
  if (!requireNamespace(p, quietly = TRUE)) {
    BiocManager::install(p, update = FALSE, ask = FALSE, lib = LIB)
  } else {
    cat("Already installed:", p, "\n")
  }
}

# NOTE: the RASCAL scripts also need the `rascal` package, which is not on CRAN or
# Bioconductor and must be installed from source separately. It is deliberately not
# installed here because the correct source repository has not been confirmed - add it
# once verified, rather than guessing at a URL.
if (!requireNamespace("rascal", quietly = TRUE)) {
  cat("\nWARNING: 'rascal' is not installed and is not handled by this script.\n")
  cat("         Install it manually before running anything in scripts/RASCAL/.\n")
}

# QDNAseq.hg38 was previously installed here but is never loaded by any script; the
# pipeline uses the custom hg38_bins_15kb_annotated.rds bin set instead. Dropped.

cat("\nVerifying...\n")
library(QDNAseq)
cat("QDNAseq:", as.character(packageVersion("QDNAseq")), "\n")

missing <- c(bioc_pkgs, cran_pkgs)[
  !vapply(c(bioc_pkgs, cran_pkgs), requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0) {
  cat("FAILED to install:", paste(missing, collapse = ", "), "\n")
} else {
  cat("All required packages present.\n")
}
