# Step 4: correlate fragmentomic metrics with ichorCNA tumour fraction
# Spearman correlation (rank-based), consistent with Step 3's non-parametric
# approach. Every correlation run twice: all 34 samples (primary) and
# excluding B2_S11 (sensitivity), per agreed methodology.

suppressMessages(library(ggplot2))

metrics_path <- "/scratch/alice/n/nhsas1/PTCL/fragmentomics/metrics/fragmentomics_metrics_summary.csv"
out_dir <- "/scratch/alice/n/nhsas1/PTCL/fragmentomics/metrics"

metrics <- read.csv(metrics_path, stringsAsFactors = FALSE)

metrics$group <- ifelse(metrics$detection_category == "Below_floor", "Below_floor",
                  ifelse(metrics$detection_category == "Detectable", "Low_TF", "High_TF"))
metrics$group <- factor(metrics$group, levels = c("Below_floor", "Low_TF", "High_TF"))

metric_cols <- c("median_fragment_length", "peak_fragment_length",
                  "short_to_long_ratio", "long_fragment_fraction")
metric_labels <- c(
    median_fragment_length = "Median fragment length (bp)",
    peak_fragment_length = "Mono-nucleosome peak position (bp)",
    short_to_long_ratio = "Short:long ratio",
    long_fragment_fraction = "Long-fragment fraction (>250bp)"
)

datasets <- list(
    All_samples = metrics,
    Excluding_B2_S11 = metrics[metrics$sample != "B2_S11", ]
)

corr_results <- data.frame(
    metric = character(), dataset = character(), n = integer(),
    spearman_rho = numeric(), p_value = numeric(),
    stringsAsFactors = FALSE
)

for (metric in metric_cols) {
    for (ds_name in names(datasets)) {
        df <- datasets[[ds_name]]
        df <- df[!is.na(df[[metric]]) & !is.na(df$TF_1Mb), ]

        ct <- suppressWarnings(cor.test(df$TF_1Mb, df[[metric]], method = "spearman"))
        corr_results <- rbind(corr_results, data.frame(
            metric = metric, dataset = ds_name, n = nrow(df),
            spearman_rho = unname(ct$estimate), p_value = ct$p.value
        ))
    }
}

# Apply BH across the four metrics within each dataset.
#
# Step 3 corrects its pairwise family but this step previously reported all eight
# correlations raw. The asymmetry between the two steps is the more visible problem: a
# marker comparing the Methods sections would find one described as BH-corrected and the
# other not, on the same cohort and the same four metrics.
#
# The family is the four metrics within one dataset, matching step 3, so All_samples and
# Excluding_B2_S11 are corrected separately rather than pooled.
#
# The three significant correlations survive correction in both datasets, so this changes
# no conclusion; it makes the reported values defensible.
#
# Ties caveat, which belongs in the Methods rather than in code: cor.test cannot compute an
# exact Spearman p-value with ties, so these are normal approximations. That matters most
# for peak_fragment_length, which takes only about six distinct values across 34 samples,
# and for long_fragment_fraction, where six samples are tied at exactly zero - see the data
# integrity note in RERUN_REQUIRED.md, since those zeros are themselves suspect.
# Initialise the column before the loop: assigning into a subset of a column that does not
# yet exist sizes it by the highest index touched, leaving a short vector.
corr_results$p_adjusted <- NA_real_
for (ds in unique(corr_results$dataset)) {
    idx <- which(corr_results$dataset == ds)
    corr_results$p_adjusted[idx] <- p.adjust(corr_results$p_value[idx], method = "BH")
}

write.csv(corr_results, file.path(out_dir, "tf_correlation_results.csv"), row.names = FALSE)
cat("TF correlation results written:", nrow(corr_results), "rows\n")
print(corr_results)

group_colors <- c("Below_floor" = "#A0B6C0", "Low_TF" = "#0063A3", "High_TF" = "#2E4057")

for (metric in metric_cols) {
    rho_all <- corr_results$spearman_rho[corr_results$metric == metric & corr_results$dataset == "All_samples"]
    p_all <- corr_results$p_value[corr_results$metric == metric & corr_results$dataset == "All_samples"]
    rho_excl <- corr_results$spearman_rho[corr_results$metric == metric & corr_results$dataset == "Excluding_B2_S11"]
    p_excl <- corr_results$p_value[corr_results$metric == metric & corr_results$dataset == "Excluding_B2_S11"]

    plot_data_all <- metrics[!is.na(metrics[[metric]]) & !is.na(metrics$TF_1Mb), ]
    plot_data_all$dataset_label <- paste0("All samples (n=", nrow(plot_data_all), ")")
    plot_data_excl <- metrics[metrics$sample != "B2_S11" & !is.na(metrics[[metric]]) & !is.na(metrics$TF_1Mb), ]
    plot_data_excl$dataset_label <- paste0("Excluding B2_S11 (n=", nrow(plot_data_excl), ")")
    plot_data <- rbind(plot_data_all, plot_data_excl)
    plot_data$dataset_label <- factor(plot_data$dataset_label, levels = c(
        unique(plot_data_all$dataset_label), unique(plot_data_excl$dataset_label)
    ))
    plot_data$flag_label <- ifelse(plot_data$sample == "B2_S11", "B2_S11", NA)

    stat_labels <- data.frame(
        dataset_label = factor(levels(plot_data$dataset_label), levels = levels(plot_data$dataset_label)),
        label = c(
            paste0("Spearman rho = ", sprintf("%.3f", rho_all), "\np = ", sprintf("%.4f", p_all)),
            paste0("Spearman rho = ", sprintf("%.3f", rho_excl), "\np = ", sprintf("%.4f", p_excl))
        )
    )

    fig <- ggplot(plot_data, aes(x = TF_1Mb, y = .data[[metric]], color = group)) +
        geom_point(size = 2.4, alpha = 0.85) +
        geom_smooth(aes(group = dataset_label), method = "loess", se = FALSE,
                     color = "grey40", linewidth = 0.6, span = 1.2) +
        geom_text(aes(label = flag_label), na.rm = TRUE, hjust = -0.25, vjust = 0.3,
                  size = 3, color = "grey20", fontface = "italic") +
        geom_text(data = stat_labels, aes(x = 0.01, y = Inf, label = label),
                  inherit.aes = FALSE, hjust = 0, vjust = 1.3, size = 3.4, color = "grey20") +
        facet_wrap(~dataset_label) +
        scale_x_log10() +
        scale_color_manual(values = group_colors, name = "TF group") +
        labs(
            title = paste0(metric_labels[metric], " vs. tumour fraction"),
            subtitle = "TF (log10 scale); B2_S11 flagged, pending final decision",
            x = "ichorCNA tumour fraction (TF_1Mb, log scale)",
            y = metric_labels[metric]
        ) +
        theme_classic(base_size = 12) +
        theme(plot.title = element_text(size = 13, face = "bold"),
              plot.subtitle = element_text(size = 9, color = "grey30"))

    fname <- paste0("fig_tfcorrelation_", metric)
    ggsave(file.path(out_dir, paste0(fname, ".pdf")), fig, width = 10, height = 6, dpi = 300)
    ggsave(file.path(out_dir, paste0(fname, ".png")), fig, width = 10, height = 6, dpi = 300)
}

cat("Done. Outputs in:", out_dir, "\n")
