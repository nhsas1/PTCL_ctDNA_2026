#!/usr/bin/env Rscript
.libPaths(c("/home/n/nhsas1/R/library", .libPaths()))

BINS_FILE  <- "/scratch/alice/n/nhsas1/PTCL/scripts/hg38_bins_15kb_annotated.rds"
OUTPUT_DIR <- "/scratch/alice/n/nhsas1/PTCL/QDNAseq_output"
METADATA   <- "/scratch/alice/n/nhsas1/PTCL/sample_metadata.csv"

suppressPackageStartupMessages({
  library(QDNAseq)
  library(Biobase)
})

pdf(NULL)

cat("=== QDNAseq Pipeline Starting ===\n")
cat("Time:", as.character(Sys.time()), "\n")
cat("QDNAseq version:", as.character(packageVersion("QDNAseq")), "\n\n")

cat("Loading hg38 bin annotations...\n")
bins <- readRDS(BINS_FILE)
cat("Bins loaded:", nrow(bins@data), "total,",
    sum(bins@data$use), "flagged for use\n\n")

metadata <- read.csv(METADATA, stringsAsFactors=FALSE)
cat("Samples to process:", nrow(metadata), "\n")

missing <- metadata$bam_path[!file.exists(metadata$bam_path)]
if (length(missing) > 0) stop("Missing BAM files:\n", paste(missing, collapse="\n"))
cat("All BAM files confirmed present.\n\n")

if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive=TRUE)

for (i in seq_len(nrow(metadata))) {

  sample_id <- metadata$sample_id[i]
  bam_path  <- metadata$bam_path[i]

  cat("─────────────────────────────────────────\n")
  cat("Sample", i, "of", nrow(metadata), ":", sample_id, "\n")
  cat("Time:", as.character(Sys.time()), "\n")

  sample_dir <- file.path(OUTPUT_DIR, sample_id)
  if (!dir.exists(sample_dir)) dir.create(sample_dir, recursive=TRUE)

  tryCatch({

    cat("  [1/6] Counting reads...\n")
    readCounts <- binReadCounts(bins, bamfiles=bam_path,
                                bamnames=sample_id)

    cat("  [2/6] Applying filters...\n")
    readCountsFiltered <- applyFilters(readCounts,
                                       residual=FALSE,
                                       blacklist=FALSE)

    cat("  [3/6] Estimating GC correction...\n")
    readCountsFiltered <- estimateCorrection(readCountsFiltered,
                                             variables=c("gc"))

    cat("  [4/6] Correcting and normalising...\n")
    copyNumbers           <- correctBins(readCountsFiltered)
    copyNumbersNormalized <- normalizeBins(copyNumbers)
    copyNumbersSmooth     <- smoothOutlierBins(copyNumbersNormalized)

    exportBins(copyNumbersSmooth,
               file=file.path(sample_dir, paste0(sample_id, ".igv")),
               format="igv")

    cat("  [5/6] Segmenting...\n")
    copyNumbersSegmented <- segmentBins(copyNumbersSmooth,
                                        transformFun="sqrt")
    copyNumbersSegmented <- normalizeSegmentedBins(copyNumbersSegmented)

    cat("  [6/6] Calling CNVs and exporting seg file...\n")
    copyNumbersCalled <- callBins(copyNumbersSegmented,
                                  method="cutoff")

    sampleNames(copyNumbersCalled) <- sample_id
    exportBins(copyNumbersCalled,
               file=file.path(sample_dir, paste0(sample_id, ".seg")),
               format="seg")

    cat("  SUCCESS:", sample_id, "\n")

  }, error=function(e) {
    try(dev.off(), silent=TRUE)
    cat("  ERROR in", sample_id, ":", conditionMessage(e), "\n")
    cat("  Skipping and continuing...\n")
  })
}

cat("\n=== Combining seg files ===\n")
seg_files <- list.files(OUTPUT_DIR, pattern="\\.seg$",
                         recursive=TRUE, full.names=TRUE)
seg_files <- seg_files[!grepl("ALL_SAMPLES", seg_files)]

if (length(seg_files) > 0) {
  all_segs <- do.call(rbind, lapply(seg_files, function(f) {
    tryCatch(read.table(f, header=TRUE, sep="\t", stringsAsFactors=FALSE),
             error=function(e) NULL)
  }))
  write.table(all_segs,
              file.path(OUTPUT_DIR, "ALL_SAMPLES_combined.seg"),
              sep="\t", quote=FALSE, row.names=FALSE)
  cat("Master seg file written with", nrow(all_segs), "segments\n")
  cat("Samples completed:", length(seg_files), "of", nrow(metadata), "\n")
} else {
  cat("WARNING: No seg files found\n")
}

sink(file.path(OUTPUT_DIR, "sessionInfo.txt"))
sessionInfo()
sink()

cat("\n=== Pipeline complete ===\n")
cat("End time:", as.character(Sys.time()), "\n")
