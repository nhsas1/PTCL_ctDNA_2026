# Step 3: fragmentomics group comparison (Below_floor / Low_TF / High_TF)
# Non-parametric testing (Kruskal-Wallis + pairwise Wilcoxon, BH-corrected).
# Every test is run twice: all 34 samples (primary) and excluding B2_S11
# (sensitivity check), per agreed methodology. B2_S11 is flagged pending
# supervisor decision, not silently excluded.

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

test_results <- data.frame(
    metric = character(), dataset = character(), test = character(),
    comparison = character(), n1 = integer(), n2 = integer(),
    statistic = numeric(), df = numeric(),
    p_value = numeric(), p_adjusted = numeric(),
    stringsAsFactors = FALSE
)

summary_stats <- data.frame(
    metric = character(), dataset = character(), group = character(),
    n = integer(), median = numeric(), q1 = numeric(), q3 = numeric(),
    stringsAsFactors = FALSE
)

for (metric in metric_cols) {
    for (ds_name in names(datasets)) {
        df <- datasets[[ds_name]]
        df <- df[!is.na(df[[metric]]), ]
        df$group <- droplevels(df$group)

        for (g in levels(df$group)) {
            vals <- df[[metric]][df$group == g]
            summary_stats <- rbind(summary_stats, data.frame(
                metric = metric, dataset = ds_name, group = g, n = length(vals),
                median = median(vals), q1 = quantile(vals, 0.25), q3 = quantile(vals, 0.75)
            ))
        }

        if (nlevels(df$group) < 3 || any(table(df$group) < 2)) {
            next
        }

        kw <- kruskal.test(df[[metric]] ~ df$group)
        test_results <- rbind(test_results, data.frame(
            metric = metric, dataset = ds_name, test = "Kruskal-Wallis",
            comparison = "omnibus (3 groups)", n1 = NA, n2 = NA,
            statistic = unname(kw$statistic), df = unname(kw$parameter),
            p_value = kw$p.value, p_adjusted = NA
        ))

        group_levels <- levels(df$group)
        pairs <- combn(group_levels, 2, simplify = FALSE)
        raw_p <- sapply(pairs, function(p) {
            x <- df[[metric]][df$group == p[1]]
            y <- df[[metric]][df$group == p[2]]
            suppressWarnings(wilcox.test(x, y)$p.value)
        })
        adj_p <- p.adjust(raw_p, method = "BH")

        for (i in seq_along(pairs)) {
            n1 <- sum(df$group == pairs[[i]][1])
            n2 <- sum(df$group == pairs[[i]][2])
            test_results <- rbind(test_results, data.frame(
                metric = metric, dataset = ds_name, test = "Wilcoxon pairwise",
                comparison = paste(pairs[[i]][1], "vs", pairs[[i]][2]),
                n1 = n1, n2 = n2, statistic = NA, df = NA,
                p_value = raw_p[i], p_adjusted = adj_p[i]
            ))
        }
    }
}

write.csv(test_results, file.path(out_dir, "group_comparison_test_results.csv"), row.names = FALSE)
write.csv(summary_stats, file.path(out_dir, "group_comparison_summary_stats.csv"), row.names = FALSE)
cat("Test results written:", nrow(test_results), "rows\n")
cat("Summary stats written:", nrow(summary_stats), "rows\n")

group_colors <- c("Below_floor" = "#A0B6C0", "Low_TF" = "#0063A3", "High_TF" = "#2E4057")

for (metric in metric_cols) {
    kw_all <- test_results$p_value[test_results$metric == metric &
                                     test_results$dataset == "All_samples" &
                                     test_results$test == "Kruskal-Wallis"]
    kw_excl <- test_results$p_value[test_results$metric == metric &
                                      test_results$dataset == "Excluding_B2_S11" &
                                      test_results$test == "Kruskal-Wallis"]

    plot_data_all <- metrics[!is.na(metrics[[metric]]), ]
    plot_data_all$dataset_label <- paste0("All samples (n=", nrow(plot_data_all), ")")
    plot_data_excl <- metrics[metrics$sample != "B2_S11" & !is.na(metrics[[metric]]), ]
    plot_data_excl$dataset_label <- paste0("Excluding B2_S11 (n=", nrow(plot_data_excl), ")")
    plot_data <- rbind(plot_data_all, plot_data_excl)
    plot_data$dataset_label <- factor(plot_data$dataset_label, levels = c(
        unique(plot_data_all$dataset_label), unique(plot_data_excl$dataset_label)
    ))
    plot_data$flag_label <- ifelse(plot_data$sample == "B2_S11", "B2_S11", NA)

    kw_labels <- data.frame(
        dataset_label = factor(levels(plot_data$dataset_label), levels = levels(plot_data$dataset_label)),
        label = c(paste0("Kruskal-Wallis p = ", sprintf("%.3f", kw_all)),
                  paste0("Kruskal-Wallis p = ", sprintf("%.3f", kw_excl)))
    )

    y_range <- range(plot_data[[metric]], na.rm = TRUE)
    y_pos <- y_range[2] + 0.08 * diff(y_range)

    fig <- ggplot(plot_data, aes(x = group, y = .data[[metric]], color = group)) +
        geom_boxplot(outlier.shape = NA, alpha = 0.15, width = 0.5) +
        geom_jitter(width = 0.12, size = 2.4, alpha = 0.85) +
        geom_text(aes(label = flag_label), na.rm = TRUE, hjust = -0.3, vjust = 0.3,
                  size = 3, color = "grey20", fontface = "italic") +
        geom_text(data = kw_labels, aes(x = 2, y = y_pos, label = label),
                  inherit.aes = FALSE, size = 3.6, color = "grey20") +
        facet_wrap(~dataset_label) +
        scale_color_manual(values = group_colors) +
        labs(
            title = paste0(metric_labels[metric], " by TF group"),
            subtitle = "Below_floor / Low-TF (Detectable) / High-TF (Moderate+High); B2_S11 flagged, pending final decision",
            x = NULL, y = metric_labels[metric]
        ) +
        theme_classic(base_size = 12) +
        theme(legend.position = "none",
              plot.title = element_text(size = 13, face = "bold"),
              plot.subtitle = element_text(size = 9, color = "grey30"))

    fname <- paste0("fig_groupcomparison_", metric)
    ggsave(file.path(out_dir, paste0(fname, ".pdf")), fig, width = 10, height = 6, dpi = 300)
    ggsave(file.path(out_dir, paste0(fname, ".png")), fig, width = 10, height = 6, dpi = 300)
}

cat("Done. Outputs in:", out_dir, "\n")
