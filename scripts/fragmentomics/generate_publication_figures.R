# Publication-ready fragmentomics diagnostic figures
# Reads from the verified Step 2 outputs (regenerated 2026-07-16 after the
# 50bp lower-bound fix). Does not recompute metrics; only visualises them.

suppressMessages(library(ggplot2))

metrics_path <- "/scratch/alice/n/nhsas1/PTCL/fragmentomics/metrics/fragmentomics_metrics_summary.csv"
comparison_path <- "/scratch/alice/n/nhsas1/PTCL/fragmentomics/metrics/batch3_dedup_vs_original_comparison.csv"
out_dir <- "/scratch/alice/n/nhsas1/PTCL/fragmentomics/metrics"

metrics <- read.csv(metrics_path, stringsAsFactors = FALSE)
comparison <- read.csv(comparison_path, stringsAsFactors = FALSE)

metrics$batch <- factor(metrics$batch, levels = c("Batch1", "Batch2", "Batch3"))
batch_colors <- c("Batch1" = "#2E4057", "Batch2" = "#0063A3", "Batch3" = "#327A95")
batch_n <- table(metrics$batch)
batch_labels <- setNames(
    paste0(names(batch_n), " (n=", as.integer(batch_n), ")"),
    names(batch_n)
)

wrap_text <- function(txt, width) paste(strwrap(txt, width = width), collapse = "\n")

# Agreement between two measurements of the same sample is not the same thing as
# correlation between them. Pearson r measures association with a straight line of any
# intercept, so two variables can differ by a large constant offset and still return
# r close to 1 - which is exactly the failure these sensitivity figures exist to exclude.
#
# Report mean bias with 95% limits of agreement (Bland-Altman) instead, which is the same
# approach already used in the QDNAseq-vs-ichorCNA concordance analysis. The delta columns
# were already computed in calculate_fragment_metrics.R and previously went unused.
agreement <- function(delta) {
    b <- mean(delta, na.rm = TRUE)
    s <- sd(delta, na.rm = TRUE)
    list(bias = b, lower = b - 1.96 * s, upper = b + 1.96 * s,
         n_pos = sum(delta > 0, na.rm = TRUE), n = sum(!is.na(delta)))
}

# Figure 1: Batch3 median fragment length, dedup sensitivity check
r1 <- cor(comparison$median_original, comparison$median_dedup, method = "pearson")
ag1 <- agreement(comparison$delta_median)
axis_range1 <- range(c(comparison$median_original, comparison$median_dedup))

fig1 <- ggplot(comparison, aes(x = median_original, y = median_dedup)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey55", linewidth = 0.6) +
    geom_point(color = "#0063A3", size = 3, alpha = 0.85) +
    coord_equal(xlim = axis_range1, ylim = axis_range1) +
    annotate("text", x = axis_range1[1], y = axis_range1[2],
             label = paste0("mean bias = ", sprintf("%+.3f", ag1$bias),
                            "\n95% LoA ", sprintf("%.3f", ag1$lower), " to ", sprintf("%.3f", ag1$upper),
                            "\nn = ", nrow(comparison)),
             hjust = 0, vjust = 1, size = 3.8, color = "grey20") +
    labs(
        title = wrap_text("Duplicate marking does not distort median fragment length", 45),
        subtitle = wrap_text("Batch 3 (n = 13): original (unmarked) vs. Picard-deduplicated BAMs", 55),
        x = "Median fragment length, original (bp)",
        y = "Median fragment length, deduplicated (bp)"
    ) +
    theme_classic(base_size = 13) +
    theme(plot.title = element_text(size = 13, face = "bold"),
          plot.subtitle = element_text(size = 10, color = "grey30"),
          plot.margin = margin(10, 15, 10, 10))

ggsave(file.path(out_dir, "fig_batch3_median_dedup_sensitivity.pdf"), fig1, width = 7.5, height = 6.5, dpi = 300)
ggsave(file.path(out_dir, "fig_batch3_median_dedup_sensitivity.png"), fig1, width = 7.5, height = 6.5, dpi = 300)

# Figure 2: Batch3 short:long ratio, dedup sensitivity check
r2 <- cor(comparison$sl_ratio_original, comparison$sl_ratio_dedup, method = "pearson")
# This is the figure where the distinction matters. delta_sl_ratio is positive in all 13
# Batch 3 samples (sign test p = 0.0002), a perfectly systematic upward shift after
# duplicate marking, yet Pearson r is 0.9998 and reports nothing about it. The bias and
# limits of agreement below are what actually support the figure's claim of a small but
# consistent effect.
ag2 <- agreement(comparison$delta_sl_ratio)
axis_range2 <- range(c(comparison$sl_ratio_original, comparison$sl_ratio_dedup))

fig2 <- ggplot(comparison, aes(x = sl_ratio_original, y = sl_ratio_dedup)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey55", linewidth = 0.6) +
    geom_point(color = "#327A95", size = 3, alpha = 0.85) +
    coord_equal(xlim = axis_range2, ylim = axis_range2) +
    annotate("text", x = axis_range2[1], y = axis_range2[2],
             label = paste0("mean bias = ", sprintf("%+.4f", ag2$bias),
                            "\n95% LoA ", sprintf("%.4f", ag2$lower), " to ", sprintf("%.4f", ag2$upper),
                            "\nhigher after dedup in ", ag2$n_pos, "/", ag2$n, " samples",
                            "\nn = ", nrow(comparison)),
             hjust = 0, vjust = 1, size = 3.8, color = "grey20") +
    labs(
        title = wrap_text("Duplicate marking has a small, consistent effect on short:long ratio", 45),
        subtitle = wrap_text("Short:long defined as 100-150bp : 151-220bp fragment counts", 55),
        x = "Short:long ratio, original",
        y = "Short:long ratio, deduplicated"
    ) +
    theme_classic(base_size = 13) +
    theme(plot.title = element_text(size = 13, face = "bold"),
          plot.subtitle = element_text(size = 10, color = "grey30"),
          plot.margin = margin(10, 15, 10, 10))

ggsave(file.path(out_dir, "fig_batch3_shorttolong_dedup_sensitivity.pdf"), fig2, width = 7.5, height = 6.5, dpi = 300)
ggsave(file.path(out_dir, "fig_batch3_shorttolong_dedup_sensitivity.png"), fig2, width = 7.5, height = 6.5, dpi = 300)

# Figure 3: Median fragment length by batch, with B2_S11 explicitly annotated
metrics$flag_label <- ifelse(metrics$sample == "B2_S11", "B2_S11", NA)

fig3 <- ggplot(metrics, aes(x = batch, y = median_fragment_length, color = batch)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.15, width = 0.5, linewidth = 0.6) +
    geom_jitter(width = 0.12, size = 2.6, alpha = 0.85) +
    geom_text(aes(label = flag_label), na.rm = TRUE, hjust = -0.35, vjust = 0.3,
              size = 3.4, color = "grey20", fontface = "italic") +
    scale_color_manual(values = batch_colors, labels = batch_labels, name = "Batch") +
    scale_x_discrete(labels = batch_labels) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.12))) +
    labs(
        title = wrap_text("Median fragment length is consistent across sequencing batches", 48),
        subtitle = wrap_text("Batch 3 shown after Picard-deduplication; B2_S11 flagged for elevated duplication and atypical fragment size profile (see Methods)", 62),
        x = NULL,
        y = "Median fragment length (bp)"
    ) +
    theme_classic(base_size = 13) +
    theme(legend.position = "none",
          plot.title = element_text(size = 13, face = "bold"),
          plot.subtitle = element_text(size = 9.5, color = "grey30"),
          plot.margin = margin(10, 15, 10, 10))

ggsave(file.path(out_dir, "fig_median_fragment_length_by_batch_annotated.pdf"), fig3, width = 8, height = 6.5, dpi = 300)
ggsave(file.path(out_dir, "fig_median_fragment_length_by_batch_annotated.png"), fig3, width = 8, height = 6.5, dpi = 300)

cat("Publication figures written to:", out_dir, "\n")
cat("Batch3 median dedup sensitivity: bias =", sprintf("%+.4f", ag1$bias),
    "| LoA", sprintf("%.4f", ag1$lower), "to", sprintf("%.4f", ag1$upper),
    "| Pearson r =", sprintf("%.4f", r1), "\n")
cat("Batch3 S:L ratio dedup sensitivity: bias =", sprintf("%+.4f", ag2$bias),
    "| LoA", sprintf("%.4f", ag2$lower), "to", sprintf("%.4f", ag2$upper),
    "| higher after dedup in", ag2$n_pos, "of", ag2$n, "samples",
    "| Pearson r =", sprintf("%.4f", r2), "\n")
