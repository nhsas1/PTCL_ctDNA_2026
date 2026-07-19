# Step 5f: check whether the near-significant (but not BH-significant) bins
# from Step 5e overlap telomeric regions, using the UCSC hg38 gap track.

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

cat("\n=== Does the significant bin (chr4:5-10Mb) overlap a telomere? ===\n")
print(overlap_summary[overlap_summary$chr == "chr4" & overlap_summary$bin_start == 5000000, ])

cat("\nTotal bins overlapping telomeres among top 20 raw p-values:",
    sum(head(overlap_summary, 20)$overlaps_telomere), "of 20\n")

write.csv(overlap_summary, file.path(out_dir, "genomewide_perbin_significance_telomere_checked.csv"), row.names = FALSE)
cat("\nDone.\n")
