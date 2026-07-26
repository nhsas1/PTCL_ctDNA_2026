# Export bin-level relative copy number for Batch 3 as plain TSV.
#
# NAMING WARNING: despite the filename, the output directory and the job name, RASCAL does
# not read these files. Both RASCAL scripts read the QDNAseq .seg and .igv files directly
# and re-linearise them with 2^log2_ratio. The only consumers of *_rascal_input.txt in this
# repository are the two plot_qdnaseq_ichorcna_concordance scripts, so what this script
# actually produces is the QDNAseq side of the QDNAseq-vs-ichorCNA concordance analysis.
# The names are being kept for now because renaming them touches paths in several scripts.
#
# This is a third independent full pass over every Batch 3 BAM, after the main pipeline and
# the plot script. It stops after smoothing and does not segment, because the concordance
# comparison is made at bin level.

.libPaths(c("/home/n/nhsas1/R/library", .libPaths()))
library(QDNAseq)
library(Biobase)
pdf(NULL)

# Only OUTPUT_DIR and METADATA differ from export_rascal_input.R; BINS_FILE is identical.
BINS_FILE  <- "/scratch/alice/n/nhsas1/PTCL/scripts/hg38_bins_15kb_annotated.rds"
OUTPUT_DIR <- "/scratch/alice/n/nhsas1/PTCL/RASCAL/input_batch3"
METADATA   <- "/scratch/alice/n/nhsas1/PTCL/scripts/sample_metadata_batch3.csv"

bins     <- readRDS(BINS_FILE)
metadata <- read.csv(METADATA, stringsAsFactors=FALSE)
if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive=TRUE)
cat("Exporting RASCAL input files for", nrow(metadata), "samples\n\n")

for (i in seq_len(nrow(metadata))) {
  sample_id <- metadata$sample_id[i]
  bam_path  <- metadata$bam_path[i]
  cat("Processing", sample_id, "(", i, "of", nrow(metadata), ")\n")
  tryCatch({
    readCounts            <- binReadCounts(bins, bamfiles=bam_path,
                                           bamnames=sample_id)
    readCountsFiltered    <- applyFilters(readCounts,
                                          residual=FALSE, blacklist=FALSE)
    readCountsFiltered    <- estimateCorrection(readCountsFiltered,
                                                variables=c("gc"))
    copyNumbers           <- correctBins(readCountsFiltered)
    copyNumbersNormalized <- normalizeBins(copyNumbers)
    copyNumbersSmooth     <- smoothOutlierBins(copyNumbersNormalized)
    # Export linear relative copy number rather than log2, so a diploid bin sits at 1.0
    # instead of 0. The concordance script compares these values against ichorCNA's
    # copy-number estimates, which are also on a linear scale.
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
