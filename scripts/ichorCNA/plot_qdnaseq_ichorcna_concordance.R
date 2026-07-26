# QDNAseq versus ichorCNA concordance - SUPERSEDED by
# plot_qdnaseq_ichorcna_concordance_v2.R, which adds the Bland-Altman agreement analysis
# and the concordance_statistics.csv output.
#
# WARNING: this script writes the same output filenames as v2
# (qdnaseq_ichorcna_concordance.pdf/png and concordance_per_sample.csv) into the same
# directory. Running it will silently overwrite the committed v2 figures with the older
# two-panel version, which reports correlation only and has no agreement analysis.
# Run v2 instead. This file is retained for provenance.
#
# See the v2 header for how the comparison works.

library(ggplot2)
library(patchwork)

qdnaseq_dir_b12 <- "/scratch/alice/n/nhsas1/PTCL/RASCAL/input"
qdnaseq_dir_b3  <- "/scratch/alice/n/nhsas1/PTCL/RASCAL/input_batch3"
ichor_dir_b12   <- "/scratch/alice/n/nhsas1/PTCL/ichorCNA/output"
ichor_dir_b3    <- "/scratch/alice/n/nhsas1/PTCL/ichorCNA/output_batch3"
output_dir      <- "/scratch/alice/n/nhsas1/ptcl_ctdna_thesis/results/ichorCNA/qdnaseq_concordance"

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

min_valid_bins <- 10

samples_b12 <- sprintf("B1_S%02d", 1:5)
samples_b12 <- c(samples_b12, sprintf("B2_S%02d", 1:16))
samples_b3  <- sprintf("B3_S%02d", 1:13)

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
  if (!file.exists(paths$qdnaseq) || !file.exists(paths$ichor)) {
    return(NULL)
  }

  qd <- read.delim(paths$qdnaseq, stringsAsFactors = FALSE)
  colnames(qd)[5] <- "value"
  qd$chromosome <- as.character(qd$chromosome)
  qd$value[qd$value == 0] <- NA
  qd$log2 <- log2(qd$value)

  ic <- read.delim(paths$ichor, stringsAsFactors = FALSE)
  colnames(ic) <- c("chr", "start", "end", "ichor_log2")
  ic$chr <- as.character(ic$chr)
  ic$ichor_log2 <- suppressWarnings(as.numeric(ic$ichor_log2))
  ic <- ic[!is.na(ic$ichor_log2), ]

  ic$qdnaseq_log2 <- NA_real_
  ic$n_bins <- 0L

  for (chrom in unique(ic$chr)) {
    ic_rows <- which(ic$chr == chrom)
    qd_chr <- qd[qd$chromosome == chrom, ]
    if (nrow(qd_chr) == 0) next
    qd_mid <- (qd_chr$start + qd_chr$end) / 2
    for (r in ic_rows) {
      in_window <- qd_mid >= ic$start[r] & qd_mid <= ic$end[r]
      vals <- qd_chr$log2[in_window]
      vals <- vals[!is.na(vals) & is.finite(vals)]
      if (length(vals) >= min_valid_bins) {
        ic$qdnaseq_log2[r] <- mean(vals)
        ic$n_bins[r] <- length(vals)
      }
    }
  }

  merged <- ic[!is.na(ic$qdnaseq_log2), c("chr", "start", "end", "ichor_log2", "qdnaseq_log2")]
  if (nrow(merged) < 50) return(NULL)
  merged$sample <- sample
  merged
}

all_samples <- c(samples_b12, samples_b3)
results_list <- list()
summary_rows <- list()

for (s in all_samples) {
  agg <- aggregate_sample(s)
  if (is.null(agg)) {
    cat("Skipped (insufficient data):", s, "\n")
    next
  }
  r_val <- cor(agg$ichor_log2, agg$qdnaseq_log2, method = "pearson")
  results_list[[s]] <- agg
  summary_rows[[s]] <- data.frame(sample = s, n_bins = nrow(agg), pearson_r = round(r_val, 4))
  cat("Processed:", s, "| bins:", nrow(agg), "| r:", round(r_val, 3), "\n")
}

per_sample <- do.call(rbind, summary_rows)
write.csv(per_sample, file.path(output_dir, "concordance_per_sample.csv"), row.names = FALSE)

pooled <- do.call(rbind, results_list)
pooled_r <- cor(pooled$ichor_log2, pooled$qdnaseq_log2, method = "pearson")
median_r <- median(per_sample$pearson_r)

scatter_plot <- ggplot(pooled, aes(x = ichor_log2, y = qdnaseq_log2)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40") +
  geom_point(size = 0.5, alpha = 0.15, color = "#1A2E4A") +
  annotate("text", x = min(pooled$ichor_log2), y = max(pooled$qdnaseq_log2),
           hjust = 0, vjust = 1,
           label = paste0("pooled r = ", round(pooled_r, 3),
                          "\nmedian per-sample r = ", round(median_r, 3),
                          "\nn = ", nrow(per_sample), " samples"),
           size = 3.6) +
  labs(
    x = "ichorCNA corrected log2 (1Mb bins)",
    y = "QDNAseq aggregated log2 (15kb to 1Mb)",
    title = "Copy number signal concordance: QDNAseq vs ichorCNA"
  ) +
  coord_equal() +
  theme_bw(base_size = 12)

per_sample$sample <- factor(per_sample$sample, levels = per_sample$sample)
bar_plot <- ggplot(per_sample, aes(x = sample, y = pearson_r)) +
  geom_col(fill = "#1A2E4A") +
  geom_hline(yintercept = median_r, linetype = "dashed", color = "#DC2626") +
  coord_cartesian(ylim = c(0, 1)) +
  labs(
    x = "Sample",
    y = "Pearson r (per sample)",
    title = "Per-sample concordance"
  ) +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 7))

combined_plot <- scatter_plot / bar_plot + plot_layout(heights = c(2, 1))

ggsave(file.path(output_dir, "qdnaseq_ichorcna_concordance.pdf"), combined_plot, width = 9, height = 10, dpi = 300)
ggsave(file.path(output_dir, "qdnaseq_ichorcna_concordance.png"), combined_plot, width = 9, height = 10, dpi = 300)

cat("\nPooled Pearson r:", round(pooled_r, 4), "\n")
cat("Median per-sample r:", round(median_r, 4), "\n")
cat("Samples processed:", nrow(per_sample), "\n")
cat("Outputs saved to:", output_dir, "\n")
