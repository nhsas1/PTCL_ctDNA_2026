# Step 2: per-sample fragmentomics metrics
# Primary cohort (34 samples): Batch1/Batch2 use standard extraction, Batch3 uses
# the deduplicated extraction (fragment_histograms_batch3_dedup/), matching the
# duplicate-handling used for Batch1/2.
# A secondary comparison table is also produced using the original (non-dedup)
# Batch3 histograms, to quantify the effect of duplicate marking on Batch3 metrics.
#
# long_fragment_fraction: proportion of fragments >250bp, used as a genomic DNA
# contamination QC indicator, independent of the short:long ctDNA ratio. Threshold
# fixed at 250bp a priori (above the existing 100-220bp short:long windows), not
# tuned to any specific sample.

suppressMessages(library(ggplot2))

hist_dir_main <- "/scratch/alice/n/nhsas1/PTCL/fragmentomics/fragment_histograms"
hist_dir_b3_dedup <- "/scratch/alice/n/nhsas1/PTCL/fragmentomics/fragment_histograms_batch3_dedup"
ichor_summary_path <- "/scratch/alice/n/nhsas1/PTCL/ichorCNA/ichorCNA_summary_integrated.csv"
out_dir <- "/scratch/alice/n/nhsas1/PTCL/fragmentomics/metrics"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

ichor <- read.csv(ichor_summary_path, stringsAsFactors = FALSE)

LONG_FRAGMENT_THRESHOLD <- 250

calc_metrics <- function(hist_path) {
    if (!file.exists(hist_path)) {
        return(list(n_fragments = NA, median_fragment_length = NA,
                    peak_fragment_length = NA, short_to_long_ratio = NA,
                    long_fragment_fraction = NA))
    }
    h <- read.delim(hist_path, stringsAsFactors = FALSE)
    if (nrow(h) == 0) {
        return(list(n_fragments = NA, median_fragment_length = NA,
                    peak_fragment_length = NA, short_to_long_ratio = NA,
                    long_fragment_fraction = NA))
    }

    h <- h[order(h$fragment_length), ]
    total <- sum(h$count)

    # Completeness check on the histogram.
    #
    # Extraction writes fragment lengths 50-1000bp, so a histogram whose maximum stops well
    # short of 1000 has been truncated - by an interrupted write, a sort spilling to a full
    # temp directory, or a chunk that failed without the submitter noticing. Truncation near
    # 250bp is invisible in n_fragments and median_fragment_length, because the affected
    # tail is a tiny fraction of reads, but it silently zeroes long_fragment_fraction.
    #
    # This matters: six samples in the committed metrics table report a long fragment
    # fraction of exactly 0 over 7.5-11.3 million fragments, which plasma cfDNA cannot
    # produce - it always carries a di-nucleosomal population near 330bp. See
    # RERUN_REQUIRED.md. Warn rather than stop, so a partial run still produces a table,
    # but make the problem visible instead of letting it pass as a real measurement.
    max_len <- max(h$fragment_length)
    if (max_len < 500) {
        warning(sprintf(
            "%s: histogram stops at %dbp, expected up to 1000bp. Likely truncated - long_fragment_fraction will be understated.",
            basename(hist_path), max_len), call. = FALSE)
        cat(sprintf("  WARNING: %s truncated at %dbp (expected 1000bp)\n",
                    basename(hist_path), max_len))
    }

    cum <- cumsum(h$count)
    median_idx <- which(cum >= total / 2)[1]
    median_len <- h$fragment_length[median_idx]

    peak_len <- h$fragment_length[which.max(h$count)]

    short_count <- sum(h$count[h$fragment_length >= 100 & h$fragment_length <= 150])
    long_count <- sum(h$count[h$fragment_length >= 151 & h$fragment_length <= 220])
    sl_ratio <- if (long_count > 0) short_count / long_count else NA

    long_frag_count <- sum(h$count[h$fragment_length > LONG_FRAGMENT_THRESHOLD])
    long_frag_fraction <- if (total > 0) long_frag_count / total else NA

    list(n_fragments = total, median_fragment_length = median_len,
         peak_fragment_length = peak_len, short_to_long_ratio = sl_ratio,
         long_fragment_fraction = long_frag_fraction)
}

primary_results <- list()
for (i in seq_len(nrow(ichor))) {
    sample <- ichor$sample[i]
    batch <- ichor$batch[i]

    if (batch == "Batch3") {
        hist_path <- file.path(hist_dir_b3_dedup, paste0(sample, "_dedup_fragment_length_hist.tsv"))
        dedup_status <- "dedup"
    } else {
        hist_path <- file.path(hist_dir_main, paste0(sample, "_fragment_length_hist.tsv"))
        dedup_status <- "source_deduplicated"
    }

    m <- calc_metrics(hist_path)
    primary_results[[sample]] <- data.frame(
        sample = sample,
        batch = batch,
        original_filename = ichor$original_filename[i],
        TF_1Mb = ichor$TF_1Mb[i],
        detection_category = ichor$detection_category[i],
        included_in_final_cohort = ichor$included_in_final_cohort[i],
        dedup_status = dedup_status,
        n_fragments = m$n_fragments,
        median_fragment_length = m$median_fragment_length,
        peak_fragment_length = m$peak_fragment_length,
        short_to_long_ratio = m$short_to_long_ratio,
        long_fragment_fraction = m$long_fragment_fraction,
        stringsAsFactors = FALSE
    )
}
primary_table <- do.call(rbind, primary_results)
write.csv(primary_table, file.path(out_dir, "fragmentomics_metrics_summary.csv"), row.names = FALSE)
cat("Primary metrics table written:", nrow(primary_table), "samples\n")
cat("Long-fragment fraction (>", LONG_FRAGMENT_THRESHOLD, "bp) summary:\n")
print(summary(primary_table$long_fragment_fraction))

b3_samples <- ichor$sample[ichor$batch == "Batch3"]
comparison_results <- list()
for (sample in b3_samples) {
    orig_path <- file.path(hist_dir_main, paste0(sample, "_fragment_length_hist.tsv"))
    dedup_path <- file.path(hist_dir_b3_dedup, paste0(sample, "_dedup_fragment_length_hist.tsv"))

    m_orig <- calc_metrics(orig_path)
    m_dedup <- calc_metrics(dedup_path)

    comparison_results[[sample]] <- data.frame(
        sample = sample,
        n_fragments_original = m_orig$n_fragments,
        n_fragments_dedup = m_dedup$n_fragments,
        median_original = m_orig$median_fragment_length,
        median_dedup = m_dedup$median_fragment_length,
        peak_original = m_orig$peak_fragment_length,
        peak_dedup = m_dedup$peak_fragment_length,
        sl_ratio_original = m_orig$short_to_long_ratio,
        sl_ratio_dedup = m_dedup$short_to_long_ratio,
        long_frag_fraction_original = m_orig$long_fragment_fraction,
        long_frag_fraction_dedup = m_dedup$long_fragment_fraction,
        delta_median = m_dedup$median_fragment_length - m_orig$median_fragment_length,
        delta_sl_ratio = m_dedup$short_to_long_ratio - m_orig$short_to_long_ratio,
        delta_long_frag_fraction = m_dedup$long_fragment_fraction - m_orig$long_fragment_fraction,
        stringsAsFactors = FALSE
    )
}
comparison_table <- do.call(rbind, comparison_results)
write.csv(comparison_table, file.path(out_dir, "batch3_dedup_vs_original_comparison.csv"), row.names = FALSE)
cat("Batch3 comparison table written:", nrow(comparison_table), "samples\n")

p1 <- ggplot(comparison_table, aes(x = median_original, y = median_dedup)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey60") +
    geom_point(color = "#0063A3", size = 3) +
    labs(title = "Batch 3: median fragment length, original vs deduplicated",
         x = "Median fragment length (original, bp)",
         y = "Median fragment length (deduplicated, bp)") +
    theme_minimal(base_size = 12)

p2 <- ggplot(comparison_table, aes(x = sl_ratio_original, y = sl_ratio_dedup)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey60") +
    geom_point(color = "#327A95", size = 3) +
    labs(title = "Batch 3: short-to-long ratio, original vs deduplicated",
         x = "Short:long ratio (original)",
         y = "Short:long ratio (deduplicated)") +
    theme_minimal(base_size = 12)

p3b <- ggplot(comparison_table, aes(x = long_frag_fraction_original, y = long_frag_fraction_dedup)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey60") +
    geom_point(color = "#2E4057", size = 3) +
    labs(title = "Batch 3: long-fragment fraction (>250bp), original vs deduplicated",
         x = "Long-fragment fraction (original)",
         y = "Long-fragment fraction (deduplicated)") +
    theme_minimal(base_size = 12)

ggsave(file.path(out_dir, "batch3_median_dedup_comparison.pdf"), p1, width = 6, height = 5, dpi = 300)
ggsave(file.path(out_dir, "batch3_median_dedup_comparison.png"), p1, width = 6, height = 5, dpi = 300)
ggsave(file.path(out_dir, "batch3_sl_ratio_dedup_comparison.pdf"), p2, width = 6, height = 5, dpi = 300)
ggsave(file.path(out_dir, "batch3_sl_ratio_dedup_comparison.png"), p2, width = 6, height = 5, dpi = 300)
ggsave(file.path(out_dir, "batch3_long_frag_fraction_dedup_comparison.pdf"), p3b, width = 6, height = 5, dpi = 300)
ggsave(file.path(out_dir, "batch3_long_frag_fraction_dedup_comparison.png"), p3b, width = 6, height = 5, dpi = 300)

p3 <- ggplot(primary_table, aes(x = batch, y = median_fragment_length, color = batch)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.3) +
    geom_jitter(width = 0.15, size = 2) +
    scale_color_manual(values = c("Batch1" = "#2E4057", "Batch2" = "#0063A3", "Batch3" = "#327A95")) +
    labs(title = "Median fragment length by batch (primary dataset)",
         x = "Batch", y = "Median fragment length (bp)") +
    theme_minimal(base_size = 12) +
    theme(legend.position = "none")

ggsave(file.path(out_dir, "median_fragment_length_by_batch.pdf"), p3, width = 6, height = 5, dpi = 300)
ggsave(file.path(out_dir, "median_fragment_length_by_batch.png"), p3, width = 6, height = 5, dpi = 300)

cat("Done. Outputs in:", out_dir, "\n")
