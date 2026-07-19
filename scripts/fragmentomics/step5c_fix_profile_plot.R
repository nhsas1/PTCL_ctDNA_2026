# Step 5c: diagnose extreme bins + fix chromosome ordering + fix color palette
# Reads the already-computed profile (no re-fitting needed).
#
# Fixes applied:
# 1. Chromosome order: gc_bins$chr was sorted lexicographically ("chr1","chr10",
#    "chr11"...) instead of numerically. Now explicitly factor-leveled chr1-chr22.
# 2. Below_floor color darkened for legibility (#A0B6C0 -> #6B8494).
# 3. B2_S11 highlight changed from black (visually collides with High_TF's
#    dark navy) to a distinct orange (#D95F02), unambiguous against every
#    other color in the palette.
# 4. Diagnostic report on the most extreme bins (mean total_count, GC%, raw
#    ratio) printed before plotting, to assess whether they are mappability/
#    reference artifacts vs. genuine signal, before any decision to exclude.

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

bin_order <- profile %>%
    distinct(chr, bin_start) %>%
    arrange(chr, bin_start) %>%
    mutate(bin_index = row_number())

profile <- merge(profile, bin_order, by = c("chr", "bin_start"))
profile <- profile[order(profile$bin_index), ]

chr_labels <- profile %>%
    group_by(chr) %>%
    summarise(mid_index = mean(bin_index), .groups = "drop")

group_summary <- profile %>%
    group_by(bin_index, chr, group) %>%
    summarise(mean_z = mean(z_ratio, na.rm = TRUE),
              se_z = sd(z_ratio, na.rm = TRUE) / sqrt(n()),
              .groups = "drop")

group_colors <- c("Below_floor" = "#6B8494", "Low_TF" = "#0063A3", "High_TF" = "#2E4057")

fig <- ggplot(group_summary, aes(x = bin_index, y = mean_z, color = group, fill = group)) +
    geom_ribbon(aes(ymin = mean_z - se_z, ymax = mean_z + se_z), alpha = 0.15, color = NA) +
    geom_line(linewidth = 0.7) +
    geom_vline(data = chr_labels[seq(1, 22, 2), ], aes(xintercept = mid_index),
               linetype = "dotted", color = "grey85", linewidth = 0.3) +
    scale_x_continuous(breaks = chr_labels$mid_index, labels = gsub("chr", "", chr_labels$chr)) +
    scale_color_manual(values = group_colors, name = "TF group") +
    scale_fill_manual(values = group_colors, name = "TF group") +
    labs(
        title = "Genome-wide fragmentation profile (DELFI-style, GC-corrected)",
        subtitle = "Short:long ratio (z-scored within sample) across 5Mb bins, chr1-22; mean +/- SE by TF group",
        x = "Chromosome", y = "Mean z-scored short:long ratio"
    ) +
    theme_classic(base_size = 12) +
    theme(axis.text.x = element_text(size = 8),
          plot.title = element_text(size = 13, face = "bold"),
          plot.subtitle = element_text(size = 9, color = "grey30"))

ggsave(file.path(out_dir, "fig_genomewide_fragmentation_profile.pdf"), fig, width = 14, height = 6, dpi = 300)
ggsave(file.path(out_dir, "fig_genomewide_fragmentation_profile.png"), fig, width = 14, height = 6, dpi = 300)

b2s11_data <- profile[profile$sample == "B2_S11", ]
if (nrow(b2s11_data) > 0) {
    fig_flagged <- fig +
        geom_line(data = b2s11_data, aes(x = bin_index, y = z_ratio),
                  inherit.aes = FALSE, color = "#D95F02", linewidth = 0.6, alpha = 0.9) +
        labs(subtitle = "Short:long ratio (z-scored within sample), chr1-22; mean +/- SE by TF group; B2_S11 shown in orange (flagged, pending final decision)")

    ggsave(file.path(out_dir, "fig_genomewide_fragmentation_profile_B2S11flagged.pdf"), fig_flagged, width = 14, height = 6, dpi = 300)
    ggsave(file.path(out_dir, "fig_genomewide_fragmentation_profile_B2S11flagged.png"), fig_flagged, width = 14, height = 6, dpi = 300)
}

cat("\nDone. Corrected figures written to:", out_dir, "\n")
