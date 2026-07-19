# Step 5: DELFI-style genome-wide fragmentation profile
# GC-correction methodology follows the standard LOESS bin-GC approach used
# across ichorCNA, DELFI, and WisecondorX (Benjamini & Speed 2012 concept),
# as implemented and validated in FinaleToolkit (Bioinformatics Advances,
# https://github.com/epifluidlab/FinaleToolkit) against the original DELFI output.
#
# Deviations from Cristiano et al. 2019, documented for Methods:
# - No ENCODE blacklist file; low-validity (N-masked) bins from the reference
#   itself (valid_base_fraction < 0.5) are excluded as a substitute, covering
#   centromere/telomere/assembly-gap regions.
# - Z-scoring is applied within-sample (across bins), not across a healthy
#   reference cohort (the original DELFI method), since no reference cohort
#   is available here. This is appropriate for visualising genome-wide shape
#   independent of each sample's overall short:long baseline (already
#   characterised in Steps 3-4), but is not equivalent to DELFI's
#   case-vs-healthy z-score classifier feature.
# - Descriptive/exploratory analysis only; no per-bin statistical testing
#   (587 bins x 3 groups would require heavy multiple-testing correction
#   not supported by this cohort size).

suppressMessages(library(ggplot2))
suppressMessages(library(dplyr))

gc_path <- "/scratch/alice/n/nhsas1/PTCL/fragmentomics/genomewide_bins/gc_content_5mb_bins.tsv"
bins_dir <- "/scratch/alice/n/nhsas1/PTCL/fragmentomics/genomewide_bins"
ichor_summary_path <- "/scratch/alice/n/nhsas1/PTCL/ichorCNA/ichorCNA_summary_integrated.csv"
out_dir <- "/scratch/alice/n/nhsas1/PTCL/fragmentomics/metrics"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

gc_bins <- read.delim(gc_path, stringsAsFactors = FALSE)
gc_bins <- gc_bins[gc_bins$valid_base_fraction >= 0.5 & !is.na(gc_bins$gc_fraction), ]
cat("GC reference bins retained after excluding low-validity regions:", nrow(gc_bins), "\n")

ichor <- read.csv(ichor_summary_path, stringsAsFactors = FALSE)
ichor$group <- ifelse(ichor$detection_category == "Below_floor", "Below_floor",
                ifelse(ichor$detection_category == "Detectable", "Low_TF", "High_TF"))
ichor$group <- factor(ichor$group, levels = c("Below_floor", "Low_TF", "High_TF"))

MIN_BIN_FRAGMENTS <- 100

process_sample <- function(sample) {
    bin_path <- file.path(bins_dir, paste0(sample, "_5Mb_bins.tsv"))
    if (!file.exists(bin_path)) return(NULL)

    b <- read.delim(bin_path, stringsAsFactors = FALSE)
    b <- merge(b, gc_bins, by = c("chr", "bin_start"))
    b <- b[b$total_count >= MIN_BIN_FRAGMENTS & b$long_count > 0, ]
    if (nrow(b) < 50) return(NULL)

    b$raw_ratio <- b$short_count / b$long_count

    loess_fit <- tryCatch(
        loess(raw_ratio ~ gc_fraction, data = b, span = 0.75),
        error = function(e) NULL
    )
    if (is.null(loess_fit)) return(NULL)

    b$fitted <- predict(loess_fit)
    b$corrected_ratio <- b$raw_ratio - b$fitted + mean(b$raw_ratio, na.rm = TRUE)
    b$z_ratio <- (b$corrected_ratio - mean(b$corrected_ratio, na.rm = TRUE)) / sd(b$corrected_ratio, na.rm = TRUE)
    b$sample <- sample
    b
}

all_results <- list()
for (i in seq_len(nrow(ichor))) {
    sample <- ichor$sample[i]
    res <- process_sample(sample)
    if (!is.null(res)) {
        res$batch <- ichor$batch[i]
        res$TF_1Mb <- ichor$TF_1Mb[i]
        res$group <- ichor$group[i]
        all_results[[sample]] <- res
    } else {
        cat("WARNING: sample", sample, "skipped (insufficient bins or fit failure)\n")
    }
}

profile_table <- do.call(rbind, all_results)
write.csv(profile_table, file.path(out_dir, "genomewide_fragmentation_profile.csv"), row.names = FALSE)
cat("Genome-wide profile table written:", nrow(profile_table), "rows across",
    length(unique(profile_table$sample)), "samples\n")

bin_order <- gc_bins[order(gc_bins$chr, gc_bins$bin_start), c("chr", "bin_start")]
bin_order$bin_index <- seq_len(nrow(bin_order))
profile_table <- merge(profile_table, bin_order, by = c("chr", "bin_start"))

chr_order <- paste0("chr", 1:22)
profile_table$chr <- factor(profile_table$chr, levels = chr_order)
profile_table <- profile_table[order(profile_table$bin_index), ]

chr_labels <- profile_table %>%
    group_by(chr) %>%
    summarise(mid_index = mean(bin_index), .groups = "drop")

group_summary <- profile_table %>%
    group_by(bin_index, chr, group) %>%
    summarise(mean_z = mean(z_ratio, na.rm = TRUE),
              se_z = sd(z_ratio, na.rm = TRUE) / sqrt(n()),
              .groups = "drop")

group_colors <- c("Below_floor" = "#A0B6C0", "Low_TF" = "#0063A3", "High_TF" = "#2E4057")

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

b2s11_data <- profile_table[profile_table$sample == "B2_S11", ]
if (nrow(b2s11_data) > 0) {
    fig_flagged <- fig +
        geom_line(data = b2s11_data, aes(x = bin_index, y = z_ratio),
                  inherit.aes = FALSE, color = "black", linewidth = 0.5, alpha = 0.6) +
        labs(subtitle = "Short:long ratio (z-scored within sample), chr1-22; mean +/- SE by TF group; B2_S11 shown in black (flagged, pending final decision)")

    ggsave(file.path(out_dir, "fig_genomewide_fragmentation_profile_B2S11flagged.pdf"), fig_flagged, width = 14, height = 6, dpi = 300)
    ggsave(file.path(out_dir, "fig_genomewide_fragmentation_profile_B2S11flagged.png"), fig_flagged, width = 14, height = 6, dpi = 300)
}

cat("Done. Outputs in:", out_dir, "\n")
