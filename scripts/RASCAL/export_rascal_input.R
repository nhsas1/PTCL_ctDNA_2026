# Export bin-level relative copy number for Batches 1-2 as plain TSV.
#
# NAMING WARNING: despite the filename and directory, RASCAL does not read these files.
# Both RASCAL scripts read the QDNAseq .seg and .igv files directly and re-linearise them
# themselves. The only consumers of *_rascal_input.txt in this repository are the two
# plot_qdnaseq_ichorcna_concordance scripts, so what this produces is the QDNAseq side of
# the QDNAseq-vs-ichorCNA concordance analysis. Names kept for now because renaming touches
# paths in several scripts; see the repository cleanup proposal.
#
# Batch 3 equivalent is scripts/QDNAseq/export_rascal_input_batch3.R, which sits in a
# different directory from this one.

.libPaths(c("/home/n/nhsas1/R/library", .libPaths()))
library(QDNAseq)
library(Biobase)

pdf(NULL)

BINS_FILE  <- "/scratch/alice/n/nhsas1/PTCL/scripts/hg38_bins_15kb_annotated.rds"
OUTPUT_DIR <- "/scratch/alice/n/nhsas1/PTCL/RASCAL/input"
METADATA   <- "/scratch/alice/n/nhsas1/PTCL/sample_metadata.csv"

bins     <- readRDS(BINS_FILE)
metadata <- read.csv(METADATA, stringsAsFactors=FALSE)

if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive=TRUE)

cat("Exporting RASCAL input files for", nrow(metadata), "samples\n\n")

for (i in seq_len(nrow(metadata))) {

  sample_id <- metadata$sample_id[i]
  bam_path  <- metadata$bam_path[i]

  cat("Processing", sample_id, "(", i, "of", nrow(metadata), ")\n")

  tryCatch({

    readCounts           <- binReadCounts(bins, bamfiles=bam_path,
                                          bamnames=sample_id)
    readCountsFiltered   <- applyFilters(readCounts,
                                         residual=FALSE, blacklist=FALSE)
    readCountsFiltered   <- estimateCorrection(readCountsFiltered,
                                               variables=c("gc"))
    copyNumbers          <- correctBins(readCountsFiltered)
    copyNumbersNormalized <- normalizeBins(copyNumbers)
    copyNumbersSmooth    <- smoothOutlierBins(copyNumbersNormalized)

    # Export linear relative copy number rather than log2, so a diploid bin sits at 1.0
    # instead of 0. The concordance script compares these against ichorCNA's copy-number
    # estimates, which are also linear.
    exportBins(copyNumbersSmooth,
               file=file.path(OUTPUT_DIR, paste0(sample_id, "_rascal_input.txt")),
               format="tsv",
               logTransform=FALSE)

    cat("  SUCCESS:", sample_id, "\n")

  }, error=function(e) {
    try(dev.off(), silent=TRUE)
    cat("  ERROR:", sample_id, "-", conditionMessage(e), "\n")
  })
}

cat("\n=== All RASCAL input files exported ===\n")
cat("Location:", OUTPUT_DIR, "\n")
