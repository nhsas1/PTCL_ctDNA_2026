# ichorCNA installation on ALICE 
# R/4.3.1, Bioconductor 3.17, Rocky Linux 9
# Run interactively on a login node:
#   module load R/4.3.1
#   R
#   source("install_ichorCNA.R")
#
# ichorCNA requires two components:
#   1. The R package (installed below)
#   2. HMMcopy readCounter binary — built separately from source,
#      see note at the bottom of this file
#
# Broadinstitute/ichorCNA: https://github.com/broadinstitute/ichorCNA

# Personal R library on ALICE
.libPaths(c("/home/n/nhsas1/R/library", .libPaths()))

# BiocManager is the package manager for Bioconductor packages
install.packages("BiocManager",
                 lib   = "/home/n/nhsas1/R/library",
                 repos = "https://cloud.r-project.org")

# HMMcopy provides the core read-depth normalisation routines that ichorCNA
# wraps. Must be installed before ichorCNA or the GitHub install will fail.
BiocManager::install("HMMcopy",
                     lib    = "/home/n/nhsas1/R/library",
                     ask    = FALSE,
                     update = FALSE)

# GenomicRanges and GenomeInfoDb handle chromosome coordinate operations
# throughout the ichorCNA pipeline
BiocManager::install(c("GenomicRanges", "GenomeInfoDb"),
                     lib    = "/home/n/nhsas1/R/library",
                     ask    = FALSE,
                     update = FALSE)

# plyr is used internally by ichorCNA for data reshaping
install.packages("plyr",
                 lib   = "/home/n/nhsas1/R/library",
                 repos = "https://cloud.r-project.org")

# devtools provides install_github()
install.packages("devtools",
                 lib   = "/home/n/nhsas1/R/library",
                 repos = "https://cloud.r-project.org")

# Install ichorCNA from GitHub. upgrade = "never" avoids silently updating
# already-installed dependencies, which can break other tools sharing this
# R library.
library(devtools)
install_github("broadinstitute/ichorCNA",
               lib     = "/home/n/nhsas1/R/library",
               upgrade = "never")

# Quick sanity check — if this loads without error, installation succeeded
library(ichorCNA)
cat("ichorCNA loaded successfully\n")
cat(paste("Version:", packageVersion("ichorCNA"), "\n"))

# ---- readCounter (HMMcopy C binary) ----------------------------------------
#
# ichorCNA requires readCounter to bin BAM reads into fixed-width windows.
# This is a C binary that must be compiled from source; it is not part of the
# R package.
#
# How it was built on ALICE:
#   cd /scratch/alice/n/nhsas1/PTCL/
#   wget https://github.com/shahcompbio/hmmcopy_utils/archive/refs/heads/master.tar.gz
#   tar -xzf master.tar.gz
#   cmake-3.27.7-linux-x86_64.sh --prefix=/scratch/alice/n/nhsas1/PTCL/cmake --skip-license
#   export PATH=/scratch/alice/n/nhsas1/PTCL/cmake/bin:$PATH
#   cd hmmcopy_utils-master && mkdir build && cd build
#   cmake .. && make
#
# Binary location after build:
#   /scratch/alice/n/nhsas1/PTCL/hmmcopy_utils-master/build/bin/readCounter
#
# ---- Reference WIG files (hg38, 1 Mb bins) ----------------------------------
#
# GC content and mappability WIG files for hg38 at 1 Mb resolution are
# bundled with the ichorCNA package under inst/extdata/. After installation,
# locate them with:
#   system.file("extdata", package = "ichorCNA")
