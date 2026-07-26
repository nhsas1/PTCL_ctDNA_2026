# Step 5c: diagnostic report on the most extreme bins.
# Reads the already-computed profile (no re-fitting needed).
#
# Reports mean total_count, GC fraction and raw ratio for the bins with the most
# extreme z-scores, so they can be assessed as mappability or reference artefacts
# versus genuine signal before any decision to exclude them. This is the evidence
# behind the centromere exclusion applied in step 5d.
#
# NO LONGER PLOTS. This script previously also re-drew the step 5 figure to correct
# three problems, writing the same output filenames as step 5 - which meant re-running
# step 5 afterwards silently restored the uncorrected version. All three fixes now live
# in step5_genomewide_profile.R itself, so that script's output is correct on its own:
#   1. Chromosome ordering, which was lexicographic (chr1, chr10, chr11 ... chr2)
#      because the bin index was built from an unfactored character column.
#   2. Below_floor colour darkened for legibility (#A0B6C0 -> #6B8494).
#   3. B2_S11 highlight changed from black, which collided with High_TF's dark navy,
#      to orange (#D95F02).

suppressMessages(library(ggplot2))
suppressMessages(library(dplyr))

out_dir <- "/scratch/alice/n/nhsas1/PTCL/fragmentomics/metrics"
profile <- read.csv(file.path(out_dir, "genomewide_fragmentation_profile.csv"), stringsAsFactors = FALSE)

chr_order <- paste0("chr", 1:22)
profile$chr <- factor(profile$chr, levels = chr_order)
profile$group <- factor(profile$group, levels = c("Below_floor", "Low_TF", "High_TF"))

bin_summary <- profile %>%
    group_by(chr, bin_start) %>%
    summarise(mean_z = mean(z_ratio, na.rm = TRUE),
              mean_total_count = mean(total_count, na.rm = TRUE),
              mean_gc = mean(gc_fraction, na.rm = TRUE),
              mean_raw_ratio = mean(raw_ratio, na.rm = TRUE),
              n_samples = n(),
              .groups = "drop") %>%
    arrange(desc(mean_z))

cat("=== Top 10 bins by mean z-scored short:long ratio (across all 34 samples) ===\n")
top10 <- head(bin_summary, 10)
for (i in seq_len(nrow(top10))) {
    cat(sprintf("%s:%d-%d  mean_z=%.2f  mean_total_count=%.0f  mean_gc=%.3f  mean_raw_ratio=%.3f  n=%d\n",
                top10$chr[i], top10$bin_start[i], top10$bin_start[i] + 5000000,
                top10$mean_z[i], top10$mean_total_count[i], top10$mean_gc[i],
                top10$mean_raw_ratio[i], top10$n_samples[i]))
}

cat("\n=== Cohort-wide bin count summary (context for the values above) ===\n")
cat("Median total_count across all bins/samples:", round(median(profile$total_count, na.rm = TRUE)), "\n")
cat("5th percentile total_count:", round(quantile(profile$total_count, 0.05, na.rm = TRUE)), "\n")
cat("Median GC fraction across all bins:", round(median(profile$gc_fraction, na.rm = TRUE), 3), "\n")
cat("GC fraction range across all bins: [", round(min(profile$gc_fraction, na.rm = TRUE), 3), ",",
    round(max(profile$gc_fraction, na.rm = TRUE), 3), "]\n")

cat("\nDiagnostic report complete. Figures are produced by step 5 and step 5d.\n")
