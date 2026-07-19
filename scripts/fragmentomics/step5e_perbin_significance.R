# Step 5e: formal per-bin significance testing, requested to replace the
# earlier unverified claim that this would be underpowered, with an actual
# tested result.

suppressMessages(library(dplyr))

out_dir <- "/scratch/alice/n/nhsas1/PTCL/fragmentomics/metrics"
profile <- read.csv(file.path(out_dir, "genomewide_fragmentation_profile_final.csv"), stringsAsFactors = FALSE)
profile$group <- factor(profile$group, levels = c("Below_floor", "Low_TF", "High_TF"))

bins <- unique(profile[, c("chr", "bin_start")])
cat("Testing", nrow(bins), "bins across 3 TF groups...\n")

results <- data.frame(chr = character(), bin_start = integer(),
                       p_value = numeric(), stringsAsFactors = FALSE)

for (i in seq_len(nrow(bins))) {
    sub <- profile[profile$chr == bins$chr[i] & profile$bin_start == bins$bin_start[i], ]
    if (length(unique(sub$group)) < 3 || any(table(sub$group) < 2)) next
    p <- tryCatch(kruskal.test(sub$z_ratio ~ sub$group)$p.value, error = function(e) NA)
    results <- rbind(results, data.frame(chr = bins$chr[i], bin_start = bins$bin_start[i], p_value = p))
}

results$p_adj <- p.adjust(results$p_value, method = "BH")
results <- results[order(results$p_value), ]

write.csv(results, file.path(out_dir, "genomewide_perbin_significance.csv"), row.names = FALSE)

cat("\nBins tested:", nrow(results), "\n")
cat("Bins with raw p < 0.05:", sum(results$p_value < 0.05, na.rm = TRUE), "\n")
cat("Bins surviving BH correction (p_adj < 0.05):", sum(results$p_adj < 0.05, na.rm = TRUE), "\n")
cat("\nTop 10 bins by raw p-value:\n")
print(head(results, 10))
