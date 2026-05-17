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

  cat('Plotting', sample_id, '(', i, 'of', nrow(metadata), ')\n')

  tryCatch({

    readCounts <- binReadCounts(bins, bamfiles=bam_path, bamnames=sample_id)

    # Plot 1 — Raw counts
    pdf(file.path(sample_dir, paste0(sample_id, '_01_rawCounts.pdf')), width=14, height=5)
    plot(readCounts, logTransform=FALSE, ylim=c(-50, 200),
         main=paste(sample_id, '- Raw Read Counts'))
    dev.off()

    readCountsFiltered <- applyFilters(readCounts, residual=FALSE, blacklist=FALSE)
    readCountsFiltered <- estimateCorrection(readCountsFiltered, variables=c('gc'))

    copyNumbers           <- correctBins(readCountsFiltered)
    copyNumbersNormalized <- normalizeBins(copyNumbers)
    copyNumbersSmooth     <- smoothOutlierBins(copyNumbersNormalized)

    # Plot 2 — Smoothed
    pdf(file.path(sample_dir, paste0(sample_id, '_02_smoothed.pdf')), width=14, height=5)
    plot(copyNumbersSmooth, main=paste(sample_id, '- Smoothed Copy Number'))
    dev.off()

    copyNumbersSegmented <- segmentBins(copyNumbersSmooth, transformFun='sqrt')
    copyNumbersSegmented <- normalizeSegmentedBins(copyNumbersSegmented)

    # Plot 3 — Segmented
    pdf(file.path(sample_dir, paste0(sample_id, '_03_segmented.pdf')), width=14, height=5)
    plot(copyNumbersSegmented, main=paste(sample_id, '- Segmented'))
    dev.off()

    copyNumbersCalled <- callBins(copyNumbersSegmented, method='cutoff')
    sampleNames(copyNumbersCalled) <- sample_id

    # Plot 4 — Called CNVs — most important for thesis
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
