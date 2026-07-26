# Concordance between the two independent copy-number callers used in this thesis:
# QDNAseq (relative, purity-naive) and ichorCNA (purity- and ploidy-aware).
#
# Agreement between two methods on the same BAMs is evidence that the copy-number signal is
# real rather than an artefact of either tool's assumptions. It is not independent
# validation - both read the same reads - but disagreement would be a red flag.
#
# The two callers work at different resolutions, so QDNAseq's 15kb bins are averaged up to
# ichorCNA's 1Mb bins before comparison. Both are placed on a log2 scale first: QDNAseq
# exports linear relative copy number centred on 1, ichorCNA exports corrected depth
# already in log2.
#
# This is version 2, which adds the Bland-Altman analysis. It writes the same output
# filenames as plot_qdnaseq_ichorcna_concordance.R, so only one of the two should be run.
# This is the version that produced the committed results.

library(ggplot2)
library(patchwork)

# NOTE: despite the directory name, these are not RASCAL inputs. See the header of
# scripts/QDNAseq/export_rascal_input_batch3.R.
qdnaseq_dir_b12 <- "/scratch/alice/n/nhsas1/PTCL/RASCAL/input"
qdnaseq_dir_b3  <- "/scratch/alice/n/nhsas1/PTCL/RASCAL/input_batch3"
ichor_dir_b12   <- "/scratch/alice/n/nhsas1/PTCL/ichorCNA/output"
ichor_dir_b3    <- "/scratch/alice/n/nhsas1/PTCL/ichorCNA/output_batch3"
output_dir      <- "/scratch/alice/n/nhsas1/ptcl_ctdna_thesis/results/ichorCNA/qdnaseq_concordance"

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# A 1Mb ichorCNA bin spans roughly 66 QDNAseq 15kb bins. Requiring at least 10 usable ones
# before averaging discards windows where most QDNAseq bins were filtered out, which would
# otherwise contribute a noisy mean from a handful of bins.
min_valid_bins <- 10

samples_b12 <- c(sprintf("B1_S%02d", 1:5), sprintf("B2_S%02d", 1:16))
samples_b3  <- sprintf("B3_S%02d", 1:13)
all_samples <- c(samples_b12, samples_b3)

get_paths <- function(sample) {
  if (startsWith(sample, "B3")) {
    q <- file.path(qdnaseq_dir_b3, paste0(sample, "_rascal_input.txt"))
    i <- file.path(ichor_dir_b3, sample, paste0(sample, ".correctedDepth.txt"))
  } else {
    q <- file.path(qdnaseq_dir_b12, paste0(sample, "_rascal_input.txt"))
    i <- file.path(ichor_dir_b12, sample, paste0(sample, ".correctedDepth.txt"))
  }
  list(qdnaseq = q, ichor = i)
}

aggregate_sample <- function(sample) {
  paths <- get_paths(sample)
  if (!file.exists(paths$qdnaseq) || !file.exists(paths$ichor)) return(NULL)

  qd <- read.delim(paths$qdnaseq, stringsAsFactors = FALSE)
  colnames(qd)[5] <- "value"
  qd$chromosome <- as.character(qd$chromosome)
  # QDNAseq values are linear relative copy number, diploid at 1.0. Zeros are set to NA
  # before the log2 so they become missing rather than -Inf, which would otherwise poison
  # the bin mean below.
  qd$value[qd$value == 0] <- NA
  qd$log2 <- log2(qd$value)

  ic <- read.delim(paths$ichor, stringsAsFactors = FALSE)
  colnames(ic) <- c("chr", "start", "end", "ichor_log2")
  ic$chr <- as.character(ic$chr)
  ic$ichor_log2 <- suppressWarnings(as.numeric(ic$ichor_log2))
  ic <- ic[!is.na(ic$ichor_log2), ]

  ic$qdnaseq_log2 <- NA_real_

  # Assign each QDNAseq bin to the ichorCNA window containing its midpoint, then average.
  # Using the midpoint rather than any overlap means a QDNAseq bin straddling a boundary is
  # counted once, in a single window, so no bin is double-weighted.
  for (chrom in unique(ic$chr)) {
    ic_rows <- which(ic$chr == chrom)
    qd_chr <- qd[qd$chromosome == chrom, ]
    if (nrow(qd_chr) == 0) next
    qd_mid <- (qd_chr$start + qd_chr$end) / 2
    for (r in ic_rows) {
      in_window <- qd_mid >= ic$start[r] & qd_mid <= ic$end[r]
      vals <- qd_chr$log2[in_window]
      vals <- vals[!is.na(vals) & is.finite(vals)]
      if (length(vals) >= min_valid_bins) ic$qdnaseq_log2[r] <- mean(vals)
    }
  }

  merged <- ic[!is.na(ic$qdnaseq_log2), c("chr", "start", "end", "ichor_log2", "qdnaseq_log2")]
  if (nrow(merged) < 50) return(NULL)
  merged$sample <- sample
  merged
}

results_list <- list()
summary_rows <- list()

for (s in all_samples) {
  agg <- aggregate_sample(s)
  if (is.null(agg)) { cat("Skipped:", s, "\n"); next }
  r_val <- cor(agg$ichor_log2, agg$qdnaseq_log2, method = "pearson")
  results_list[[s]] <- agg
  summary_rows[[s]] <- data.frame(sample = s, n_bins = nrow(agg), pearson_r = round(r_val, 4))
  cat("Processed:", s, "| bins:", nrow(agg), "| r:", round(r_val, 3), "\n")
}

per_sample <- do.call(rbind, summary_rows)
pooled <- do.call(rbind, results_list)

pooled$mean_log2 <- (pooled$ichor_log2 + pooled$qdnaseq_log2) / 2
pooled$diff_log2 <- pooled$qdnaseq_log2 - pooled$ichor_log2

# Bland-Altman on the pooled bins. Bias near zero means neither caller systematically reads
# higher; the limits of agreement bound how far an individual bin can differ.
#
# CAVEAT for the write-up: bins are pooled across all samples and treated as independent,
# but they are clustered within sample. The limits of agreement are therefore narrower than
# a mixed-effects treatment would give, and should be read as descriptive rather than as a
# formal interval. The per-sample correlations in panel C are the honest per-sample view.
bias <- mean(pooled$diff_log2)
sd_diff <- sd(pooled$diff_log2)
upper_loa <- bias + 1.96 * sd_diff
lower_loa <- bias - 1.96 * sd_diff

# Regressing the difference on the mean tests for proportional bias, i.e. whether the two
# callers diverge more at high or low copy number. A slope near zero means the agreement
# holds across the range rather than only near diploid.
prop_bias_fit <- lm(diff_log2 ~ mean_log2, data = pooled)
prop_slope <- coef(prop_bias_fit)[2]
prop_p <- summary(prop_bias_fit)$coefficients[2, 4]

pct_within <- mean(pooled$diff_log2 >= lower_loa & pooled$diff_log2 <= upper_loa) * 100

pooled_r <- cor(pooled$ichor_log2, pooled$qdnaseq_log2, method = "pearson")
median_r <- median(per_sample$pearson_r)

write.csv(per_sample, file.path(output_dir, "concordance_per_sample.csv"), row.names = FALSE)

stats_out <- data.frame(
  metric = c("pooled_pearson_r", "median_per_sample_r", "bland_altman_bias",
             "sd_of_differences", "upper_LoA", "lower_LoA",
             "proportional_bias_slope", "proportional_bias_pvalue", "pct_within_LoA", "n_samples", "n_bins_pooled"),
  value = c(round(pooled_r,4), round(median_r,4), round(bias,4),
            round(sd_diff,4), round(upper_loa,4), round(lower_loa,4),
            round(prop_slope,5), signif(prop_p,3), round(pct_within,2), nrow(per_sample), nrow(pooled))
)
write.csv(stats_out, file.path(output_dir, "concordance_statistics.csv"), row.names = FALSE)

scatter_plot <- ggplot(pooled, aes(x = ichor_log2, y = qdnaseq_log2)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40") +
  geom_point(size = 0.4, alpha = 0.12, color = "#1A2E4A") +
  annotate("text", x = min(pooled$ichor_log2), y = max(pooled$qdnaseq_log2),
           hjust = 0, vjust = 1,
           label = paste0("pooled r = ", round(pooled_r, 3),
                          "\nmedian per-sample r = ", round(median_r, 3),
                          "\nn = ", nrow(per_sample), " samples"),
           size = 3.4) +
  labs(x = "ichorCNA corrected log2 (1Mb bins)",
       y = "QDNAseq aggregated log2 (15kb to 1Mb)",
       title = "A. Signal correlation") +
  coord_equal() +
  theme_bw(base_size = 11)

ba_plot <- ggplot(pooled, aes(x = mean_log2, y = diff_log2)) +
  geom_point(size = 0.4, alpha = 0.10, color = "#1A2E4A") +
  geom_hline(yintercept = 0, color = "grey70") +
  geom_hline(yintercept = bias, linetype = "solid", color = "#DC2626") +
  geom_hline(yintercept = upper_loa, linetype = "dashed", color = "#D97706") +
  geom_hline(yintercept = lower_loa, linetype = "dashed", color = "#D97706") +
  annotate("text", x = max(pooled$mean_log2), y = upper_loa, vjust = -0.4, hjust = 1,
           label = paste0("+1.96 SD = ", round(upper_loa, 3)), size = 3.0, color = "#D97706") +
  annotate("text", x = max(pooled$mean_log2), y = lower_loa, vjust = 1.3, hjust = 1,
           label = paste0("-1.96 SD = ", round(lower_loa, 3)), size = 3.0, color = "#D97706") +
  annotate("text", x = max(pooled$mean_log2), y = bias, vjust = -0.4, hjust = 1,
           label = paste0("bias = ", round(bias, 3)), size = 3.0, color = "#DC2626") +
  labs(x = "Mean log2 of the two methods",
       y = "Difference (QDNAseq minus ichorCNA)",
       title = "B. Bland-Altman agreement") +
  theme_bw(base_size = 11)

per_sample$sample <- factor(per_sample$sample, levels = per_sample$sample)
bar_plot <- ggplot(per_sample, aes(x = sample, y = pearson_r)) +
  geom_col(fill = "#1A2E4A") +
  geom_hline(yintercept = median_r, linetype = "dashed", color = "#DC2626") +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = "Sample", y = "Pearson r", title = "C. Per-sample concordance") +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 7))

combined_plot <- (scatter_plot | ba_plot) / bar_plot + plot_layout(heights = c(2, 1))

ggsave(file.path(output_dir, "qdnaseq_ichorcna_concordance.pdf"), combined_plot, width = 12, height = 10, dpi = 300)
ggsave(file.path(output_dir, "qdnaseq_ichorcna_concordance.png"), combined_plot, width = 12, height = 10, dpi = 300)

cat("\n=== Bland-Altman summary ===\n")
cat("Bias:", round(bias,4), "log2 units\n")
cat("Limits of agreement:", round(lower_loa,4), "to", round(upper_loa,4), "\n")
cat("Proportional bias slope:", round(prop_slope,5), "| p =", signif(prop_p,3), "\n")
cat("Percent within LoA:", round(pct_within,2), "\n")
cat("Pooled r:", round(pooled_r,4), "| Median per-sample r:", round(median_r,4), "\n")
cat("Outputs saved to:", output_dir, "\n")
