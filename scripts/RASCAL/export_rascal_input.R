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

    # Export in non-log2 format for RASCAL
    # This gives relative copy numbers centred on 1.0 (diploid = 1.0)
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
