# Step 5f: check whether the near-significant (but not BH-significant) bins
# from Step 5e overlap telomeric regions, using the UCSC hg38 gap track.
#
# THIS TEST DOES NOT MEASURE WHAT IT APPEARS TO. Read before using its output.
#
# In hg38 the gap track's telomere records are exactly 10kb at each chromosome end. A 5Mb
# bin "overlapping" a 10kb telomere contains 0.2% telomeric sequence. Worse, those regions
# are N-masked in the reference, so they carry no mappable reads and contribute no fragment
# signal to the bin at all. The overlap is with a region that is invisible to the assay by
# construction.
#
# What overlaps_telomere actually identifies is therefore "the first or last bin of a
# chromosome". The committed output confirms this exactly: 39 of 534 bins are flagged,
# comprising 17 p-arm terminal bins (five chromosomes lost theirs to the centromere
# exclusion) and 22 q-arm terminal bins - two per chromosome.
#
# Consequently no conclusion of the form "k of the top bins overlap telomeres, therefore
# the signal is a telomeric artefact" is supportable from this script. The flag carries
# essentially no information about telomeric sequence content.
#
# If the intended hypothesis is a SUBTELOMERIC effect - which is biologically reasonable,
# since subtelomeric regions differ in chromatin state and replication timing - the test
# needs a real subtelomeric window, for example the terminal 2-5Mb, or a mappability and GC
# covariate. That is a different analysis, so it is flagged rather than rewritten here.
#
# Note also the more likely explanation for the top bins lies elsewhere: chr19 is
# significantly over-represented among them (see the comment in step5_genomewide_profile.R),
# which points at residual GC bias rather than position along the chromosome.

out_dir <- "/scratch/alice/n/nhsas1/PTCL/fragmentomics/metrics"
gap_path <- "/scratch/alice/n/nhsas1/PTCL/reference/gap.txt"

gap_raw <- read.delim(gap_path, header = FALSE, stringsAsFactors = FALSE)
colnames(gap_raw)[2:8] <- c("chr", "chromStart", "chromEnd", "ix", "n_or_bridged", "size", "type")

telomeres <- gap_raw[gap_raw$type == "telomere", c("chr", "chromStart", "chromEnd")]
cat("=== Telomere entries (hg38, chr1-22) ===\n")
print(telomeres[telomeres$chr %in% paste0("chr", 1:22), ])

results <- read.csv(file.path(out_dir, "genomewide_perbin_significance.csv"), stringsAsFactors = FALSE)
results$bin_end <- results$bin_start + 5000000

results <- merge(results, telomeres, by = "chr", all.x = TRUE)
results$overlaps_telomere <- with(results,
    !is.na(chromStart) & bin_start < chromEnd & bin_end > chromStart)

overlap_summary <- aggregate(overlaps_telomere ~ chr + bin_start + p_value + p_adj,
                              data = results, FUN = any)
overlap_summary <- overlap_summary[order(overlap_summary$p_value), ]

cat("\n=== Top 15 bins by raw p-value, now flagged for telomere overlap ===\n")
print(head(overlap_summary, 15))

# The identity of the significant bin is hardcoded here. If step5e is re-run - which it must
# be, since its inputs changed when the z-scores were corrected - this block will report on
# whichever bin used to be significant, or print an empty result with no error. Report the
# actual top bin instead of a remembered one.
top_bin <- overlap_summary[which.min(overlap_summary$p_adj), ]
cat("\n=== Most significant bin, and whether it is a terminal bin ===\n")
print(top_bin)
cat("(For the committed results this was chr4:5-10Mb, which could not have overlapped a\n")
cat(" telomere in any case: the chr4 p-telomere spans 0-10kb and this bin starts at 5Mb.)\n")

cat("\nTerminal bins among top 20 raw p-values:",
    sum(head(overlap_summary, 20)$overlaps_telomere), "of 20\n")
cat("Terminal bins genome-wide:", sum(overlap_summary$overlaps_telomere),
    "of", nrow(overlap_summary), "\n")
cat("(Compare the two proportions. Reported as terminal-bin enrichment, which is what the\n")
cat(" flag measures - not telomeric sequence content. See the header.)\n")

write.csv(overlap_summary, file.path(out_dir, "genomewide_perbin_significance_telomere_checked.csv"), row.names = FALSE)
cat("\nDone.\n")
