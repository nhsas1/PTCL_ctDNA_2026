#!/usr/bin/env Rscript

# QDNAseq copy-number pipeline for Batches 1 and 2 (n = 21).
# Produces per-sample relative copy-number profiles from sWGS BAMs, exported as
# .igv (bin-level) and .seg (segment-level) files. These feed the QDNAseq-vs-ichorCNA
# concordance check; ichorCNA, not QDNAseq, is the primary CNA caller for this thesis
# because ichorCNA jointly estimates tumour fraction and ploidy, which QDNAseq does not.
# Batch 3 is handled by the separate QDNAseq_batch3.R because its BAMs are not
# deduplicated at source.

.libPaths(c("/home/n/nhsas1/R/library", .libPaths()))

# 15kb fixed-width bins with GC content and blacklist annotation, pre-built for hg38.
# NOTE: the script that generated this .rds is not in the repository. It is required by
# every QDNAseq and RASCAL export script here and cannot currently be regenerated.
BINS_FILE  <- "/scratch/alice/n/nhsas1/PTCL/scripts/hg38_bins_15kb_annotated.rds"
OUTPUT_DIR <- "/scratch/alice/n/nhsas1/PTCL/QDNAseq_output"
METADATA   <- "/scratch/alice/n/nhsas1/PTCL/sample_metadata.csv"

suppressPackageStartupMessages({
  library(QDNAseq)
  library(Biobase)
})

# Open a null graphics device. QDNAseq functions plot as a side effect; on a headless
# SLURM node with no device open this would otherwise scatter Rplots.pdf into the CWD.
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

# Fail fast if any BAM is missing, before spending hours on the samples that do exist.
# This is deliberately stricter than the per-sample tryCatch inside the loop below: a
# missing input file is a setup error worth aborting on, whereas a mid-run failure on one
# sample should not cost the whole batch.
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

    # Tally reads falling in each 15kb bin. Note two defaults that are left unset and
    # matter for interpretation:
    #   - minMapq defaults to 37. The fragmentomics arm filters at -q 20, so the two
    #     analysis arms of this thesis use different mapping-quality thresholds.
    #   - pairedEnds is not set, so reads are counted individually rather than as
    #     fragments. Acceptable for relative copy number, where only the ratio between
    #     bins matters, but it is not a fragment count.
    cat("  [1/6] Counting reads...\n")
    readCounts <- binReadCounts(bins, bamfiles=bam_path,
                                bamnames=sample_id)

    # Drop bins unsuitable for copy-number calling. Both of QDNAseq's artefact filters are
    # switched off here: residual=FALSE keeps bins whose coverage deviates from the
    # reference expectation, blacklist=FALSE keeps known problem regions. Mappability is
    # also never supplied, so no mappability correction is applied anywhere in this
    # pipeline. This is the concrete form of the mappability limitation stated in the
    # thesis, and it means low-mappability bins carry artificially low counts.
    #
    # chrX and chrY are excluded here, but only IMPLICITLY: applyFilters defaults to
    # chromosomes=c("X","Y") and that argument is never overridden. Sex-chromosome
    # exclusion is appropriate for a mixed-sex cohort (X copy number differs by sex and
    # would otherwise read as a CNA), but passing chromosomes= explicitly would silently
    # turn it off. Do not add that argument without also naming X and Y.
    cat("  [2/6] Applying filters...\n")
    readCountsFiltered <- applyFilters(readCounts,
                                       residual=FALSE,
                                       blacklist=FALSE)

    # Fit the GC-content bias curve. cfDNA library prep and PCR amplify GC-rich and
    # GC-poor fragments unevenly, which produces coverage waves that mimic broad CNAs if
    # left uncorrected.
    cat("  [3/6] Estimating GC correction...\n")
    readCountsFiltered <- estimateCorrection(readCountsFiltered,
                                             variables=c("gc"))

    # Apply the GC fit, median-centre so a diploid bin sits at log2 ratio 0, then damp
    # single-bin spikes. Smoothing matters at ~1x coverage, where per-bin counts are low
    # enough that Poisson noise alone produces isolated extreme bins.
    cat("  [4/6] Correcting and normalising...\n")
    copyNumbers           <- correctBins(readCountsFiltered)
    copyNumbersNormalized <- normalizeBins(copyNumbers)
    copyNumbersSmooth     <- smoothOutlierBins(copyNumbersNormalized)

    # Bin-level profile, exported before segmentation. This is the file the
    # QDNAseq-vs-ichorCNA concordance analysis reads.
    exportBins(copyNumbersSmooth,
               file=file.path(sample_dir, paste0(sample_id, ".igv")),
               format="igv")

    # Join adjacent bins of similar value into segments (CBS). The sqrt transform
    # stabilises variance: raw counts are Poisson, so bin variance scales with the mean,
    # and sqrt makes the spread roughly constant across the coverage range so CBS does not
    # over-segment high-coverage regions.
    cat("  [5/6] Segmenting...\n")
    copyNumbersSegmented <- segmentBins(copyNumbersSmooth,
                                        transformFun="sqrt")
    copyNumbersSegmented <- normalizeSegmentedBins(copyNumbersSegmented)

    # Assign discrete gain/loss/neutral calls by fixed log2 cutoff. QDNAseq calls are
    # relative and purity-naive; ichorCNA supplies the absolute copy number and tumour
    # fraction used downstream.
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

# Concatenate the per-sample seg files into one cohort-level table.
cat("\n=== Combining seg files ===\n")
seg_files <- list.files(OUTPUT_DIR, pattern="\\.seg$",
                         recursive=TRUE, full.names=TRUE)
# Exclude the master file itself so re-runs do not concatenate the previous output.
seg_files <- seg_files[!grepl("ALL_SAMPLES", seg_files)]

if (length(seg_files) > 0) {
  # A seg file that fails to parse yields NULL and is dropped by rbind. Previously the
  # completion count below reported length(seg_files) regardless, so a truncated or
  # malformed file was silently excluded from the master table while still being counted
  # as a completed sample. Track which files actually parsed and report that instead.
  parsed <- lapply(seg_files, function(f) {
    tryCatch(read.table(f, header=TRUE, sep="\t", stringsAsFactors=FALSE),
             error=function(e) NULL)
  })
  failed <- seg_files[vapply(parsed, is.null, logical(1))]
  if (length(failed) > 0) {
    cat("WARNING:", length(failed), "seg file(s) failed to parse and are NOT in the",
        "master file:\n")
    cat(paste0("  ", failed, collapse="\n"), "\n")
  }

  all_segs <- do.call(rbind, parsed)
  write.table(all_segs,
              file.path(OUTPUT_DIR, "ALL_SAMPLES_combined.seg"),
              sep="\t", quote=FALSE, row.names=FALSE)
  cat("Master seg file written with", nrow(all_segs), "segments\n")
  cat("Samples completed:", length(seg_files) - length(failed), "of", nrow(metadata), "\n")
} else {
  cat("WARNING: No seg files found\n")
}

sink(file.path(OUTPUT_DIR, "sessionInfo.txt"))
sessionInfo()
sink()

cat("\n=== Pipeline complete ===\n")
cat("End time:", as.character(Sys.time()), "\n")
