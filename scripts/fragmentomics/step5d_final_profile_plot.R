# Step 5d: exclude centromere-overlapping bins, produce final genome-wide profile
#
# Exclusion basis: UCSC hg38 centromere track (centromeres.txt), downloaded
# 2026-07-17 from https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/centromeres.txt.gz
# A 5Mb bin is excluded if it overlaps the centromeric span (min start to max
# end across all contig records) for its chromosome. This substitutes for the
# ENCODE blacklist used in the original DELFI pipeline (Cristiano et al. 2019);
# it is a narrower, coordinate-verified approximation covering centromeric
# regions specifically, not the full blacklist (segmental duplications,
# additional low-mappability regions are not covered). Four bins (chr3, chr4,
# chr5, chr10) were independently flagged as extreme (z>3), group-invariant
# outliers prior to this exclusion, then confirmed to overlap centromeres.

suppressMessages(library(ggplot2))
suppressMessages(library(dplyr))

out_dir <- "/scratch/alice/n/nhsas1/PTCL/fragmentomics/metrics"
cen_path <- "/scratch/alice/n/nhsas1/PTCL/reference/centromeres.txt"

profile <- read.csv(file.path(out_dir, "genomewide_fragmentation_profile.csv"), stringsAsFactors = FALSE)

cen_raw <- read.delim(cen_path, header = FALSE, stringsAsFactors = FALSE,
                       col.names = c("bin", "chr", "chromStart", "chromEnd", "contig"))
centromeres <- cen_raw %>%
    group_by(chr) %>%
    summarise(cen_start = min(chromStart), cen_end = max(chromEnd), .groups = "drop")

cat("=== Centromeric spans (chr1-22) ===\n")
print(centromeres)

profile <- merge(profile, centromeres, by = "chr", all.x = TRUE)
profile$bin_end <- profile$bin_start + 5000000
profile$overlaps_centromere <- with(profile,
    !is.na(cen_start) & bin_start < cen_end & bin_end > cen_start)

n_excluded_bins <- length(unique(paste(profile$chr[profile$overlaps_centromere],
                                        profile$bin_start[profile$overlaps_centromere])))
cat("\nBins excluded as centromere-overlapping:", n_excluded_bins, "\n")
excluded_summary <- profile %>%
    filter(overlaps_centromere) %>%
    distinct(chr, bin_start) %>%
    arrange(chr, bin_start)
print(excluded_summary)

profile_clean <- profile[!profile$overlaps_centromere, ]
cat("\nRows retained after exclusion:", nrow(profile_clean), "of", nrow(profile), "\n")

chr_order <- paste0("chr", 1:22)
profile_clean$chr <- factor(profile_clean$chr, levels = chr_order)
profile_clean$group <- factor(profile_clean$group, levels = c("Below_floor", "Low_TF", "High_TF"))

bin_order <- profile_clean %>%
    distinct(chr, bin_start) %>%
    arrange(chr, bin_start) %>%
    mutate(bin_index = row_number())

profile_clean <- merge(profile_clean, bin_order, by = c("chr", "bin_start"))
profile_clean <- profile_clean[order(profile_clean$bin_index), ]

chr_labels <- profile_clean %>%
    group_by(chr) %>%
    summarise(mid_index = mean(bin_index), .groups = "drop")

group_summary <- profile_clean %>%
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
        subtitle = paste0("Short:long ratio (z-scored within sample) across 5Mb bins, chr1-22, ",
                           "excl. centromeric bins (n=", n_excluded_bins, "); mean +/- SE by TF group"),
        x = "Chromosome", y = "Mean z-scored short:long ratio"
    ) +
    theme_classic(base_size = 12) +
    theme(axis.text.x = element_text(size = 8),
          plot.title = element_text(size = 13, face = "bold"),
          plot.subtitle = element_text(size = 9, color = "grey30"))

ggsave(file.path(out_dir, "fig_genomewide_fragmentation_profile_final.pdf"), fig, width = 14, height = 6, dpi = 300)
ggsave(file.path(out_dir, "fig_genomewide_fragmentation_profile_final.png"), fig, width = 14, height = 6, dpi = 300)

b2s11_data <- profile_clean[profile_clean$sample == "B2_S11", ]
if (nrow(b2s11_data) > 0) {
    fig_flagged <- fig +
        geom_line(data = b2s11_data, aes(x = bin_index, y = z_ratio),
                  inherit.aes = FALSE, color = "#D95F02", linewidth = 0.6, alpha = 0.9) +
        labs(subtitle = paste0("Short:long ratio (z-scored within sample), excl. centromeric bins (n=",
                                n_excluded_bins, "); mean +/- SE by TF group; B2_S11 in orange (flagged, pending decision)"))

    ggsave(file.path(out_dir, "fig_genomewide_fragmentation_profile_final_B2S11flagged.pdf"), fig_flagged, width = 14, height = 6, dpi = 300)
    ggsave(file.path(out_dir, "fig_genomewide_fragmentation_profile_final_B2S11flagged.png"), fig_flagged, width = 14, height = 6, dpi = 300)
}

write.csv(profile_clean, file.path(out_dir, "genomewide_fragmentation_profile_final.csv"), row.names = FALSE)

cat("\nDone. Final figures and cleaned profile table written to:", out_dir, "\n")
