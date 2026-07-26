# Diagnostic QDNAseq plots for Batches 1 and 2, four per sample, showing the profile at
# each stage of processing. These are QC figures used to eyeball whether GC correction
# and segmentation behaved sensibly; they are not the thesis CNA figures, which come from
# ichorCNA and GenVisR.
#
# This script re-runs the whole QDNAseq pipeline from the BAM rather than reloading the
# objects saved by QDNAseq_ALICE_automated.R. That makes it a second independent pass over
# every BAM, so the plots are only guaranteed to match the committed .seg files as long as
# both scripts keep identical parameters. Any parameter change must be made in both.

.libPaths(c('/home/n/nhsas1/R/library', .libPaths()))
library(QDNAseq)
library(Biobase)

pdf(NULL)

BINS_FILE  <- '/scratch/alice/n/nhsas1/PTCL/scripts/hg38_bins_15kb_annotated.rds'
OUTPUT_DIR <- '/scratch/alice/n/nhsas1/PTCL/QDNAseq_output'
METADATA   <- '/scratch/alice/n/nhsas1/PTCL/sample_metadata.csv'

bins     <- readRDS(BINS_FILE)
metadata <- read.csv(METADATA, stringsAsFactors=FALSE)

cat('Generating plots for', nrow(metadata), 'samples\n')

for (i in seq_len(nrow(metadata))) {

  sample_id  <- metadata$sample_id[i]
  bam_path   <- metadata$bam_path[i]
  sample_dir <- file.path(OUTPUT_DIR, sample_id)
  # The main pipeline normally creates this directory, but it will be absent for any
  # sample that failed there. Without it pdf() errors, the error is swallowed by the
  # tryCatch below, and the sample ends up with no diagnostic plots at all - which is
  # exactly the sample whose plots are most worth looking at.
  if (!dir.exists(sample_dir)) dir.create(sample_dir, recursive=TRUE)

  cat('Plotting', sample_id, '(', i, 'of', nrow(metadata), ')\n')

  tryCatch({

    readCounts <- binReadCounts(bins, bamfiles=bam_path, bamnames=sample_id)

    # Plot 1 — Raw counts, before any correction. Used to spot gross coverage problems:
    # a sample with heavy duplication or a failed library shows up here as an unusually
    # wide or skewed spread of per-bin counts.
    pdf(file.path(sample_dir, paste0(sample_id, '_01_rawCounts.pdf')), width=14, height=5)
    plot(readCounts, logTransform=FALSE, ylim=c(-50, 200),
         main=paste(sample_id, '- Raw Read Counts'))
    dev.off()

    readCountsFiltered <- applyFilters(readCounts, residual=FALSE, blacklist=FALSE)
    readCountsFiltered <- estimateCorrection(readCountsFiltered, variables=c('gc'))

    copyNumbers           <- correctBins(readCountsFiltered)
    copyNumbersNormalized <- normalizeBins(copyNumbers)
    copyNumbersSmooth     <- smoothOutlierBins(copyNumbersNormalized)

    # Plot 2 — After GC correction, normalisation and outlier smoothing. Comparing this
    # against plot 1 shows whether the GC wave was removed.
    pdf(file.path(sample_dir, paste0(sample_id, '_02_smoothed.pdf')), width=14, height=5)
    plot(copyNumbersSmooth, main=paste(sample_id, '- Smoothed Copy Number'))
    dev.off()

    copyNumbersSegmented <- segmentBins(copyNumbersSmooth, transformFun='sqrt')
    copyNumbersSegmented <- normalizeSegmentedBins(copyNumbersSegmented)

    # Plot 3 — After CBS segmentation. Checked for over-segmentation, which at ~1x
    # coverage indicates the variance-stabilising transform is not holding.
    pdf(file.path(sample_dir, paste0(sample_id, '_03_segmented.pdf')), width=14, height=5)
    plot(copyNumbersSegmented, main=paste(sample_id, '- Segmented'))
    dev.off()

    copyNumbersCalled <- callBins(copyNumbersSegmented, method='cutoff')
    sampleNames(copyNumbersCalled) <- sample_id

    # Plot 4 — Discrete gain/loss calls. The most useful of the four for review, but note
    # these calls are relative and purity-naive: a low tumour fraction sample shows few
    # calls because the signal is diluted, not because the genome is quiet.
    pdf(file.path(sample_dir, paste0(sample_id, '_04_called.pdf')), width=14, height=5)
    plot(copyNumbersCalled, main=paste(sample_id, '- Called CNVs'))
    dev.off()

    cat('  Done:', sample_id, '\n')

  }, error=function(e) {
    try(dev.off(), silent=TRUE)
    cat('  ERROR:', sample_id, '-', conditionMessage(e), '\n')
  })
}
cat('All plots complete\n')
