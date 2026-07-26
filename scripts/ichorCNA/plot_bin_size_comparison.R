# Validation that the 1Mb bin size used throughout the CNA arm does not drive the tumour
# fraction estimates. ichorCNA was re-run at 500kb on all 34 samples and the two sets of
# estimates are compared here.
#
# This matters because bin size is a free parameter: smaller bins give finer resolution but
# fewer reads per bin, so at ~1x coverage they are noisier. If tumour fraction moved
# substantially with bin size, the 3% detection floor used to define the analysable cohort
# would be arbitrary. Two views are produced - a scatter against the identity line, and a
# Bland-Altman plot, which is the correct way to assess agreement rather than association.
#
# NOTE: these paths use a different project root (ptcl_ctdna_thesis) from every other script
# in this pipeline, which use PTCL. Worth reconciling.

library(ggplot2)
library(patchwork)

data_path <- "/scratch/alice/n/nhsas1/ptcl_ctdna_thesis/results/ichorCNA/bin_size_validation/ichorCNA_bin_size_comparison.csv"
output_dir <- "/scratch/alice/n/nhsas1/ptcl_ctdna_thesis/results/ichorCNA/bin_size_validation"

df <- read.csv(data_path, stringsAsFactors = FALSE)

# Unflagged samples come through as empty strings; label them explicitly so they appear in
# the legend rather than as a blank category.
df$flag[df$flag == ""] <- "Consistent"
df$flag <- factor(df$flag, levels = c("Consistent", "THRESHOLD_FLIP", "MAJOR_OUTLIER"))

flag_colors <- c(
  "Consistent" = "#64748B",
  "THRESHOLD_FLIP" = "#D97706",
  "MAJOR_OUTLIER" = "#DC2626"
)

# Pearson is reported alongside the Bland-Altman panel rather than instead of it. On its
# own r would be misleading here: two estimates can correlate almost perfectly and still
# differ by a constant offset, which is exactly the failure this validation must exclude.
# The bias and limits of agreement below are what answer that.
r_value <- cor(df$TF_1Mb, df$TF_500kb, method = "pearson")
n_samples <- nrow(df)

scatter_plot <- ggplot(df, aes(x = TF_1Mb, y = TF_500kb, color = flag)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40") +
  # The dotted lines mark the 3% detection floor on both axes. A sample falling in an
  # off-diagonal quadrant crossed the floor when bin size changed, which is the
  # THRESHOLD_FLIP category and the practical risk this validation is testing for.
  geom_hline(yintercept = 0.03, linetype = "dotted", color = "grey60") +
  geom_vline(xintercept = 0.03, linetype = "dotted", color = "grey60") +
  geom_point(size = 2.5, alpha = 0.85) +
  scale_color_manual(values = flag_colors, name = "Category") +
  annotate("text", x = 0.02, y = 0.62, hjust = 0,
           label = paste0("r = ", round(r_value, 3), "\nn = ", n_samples),
           size = 3.6) +
  labs(
    x = "Tumour fraction, 1Mb bins",
    y = "Tumour fraction, 500kb bins",
    title = "ichorCNA tumour fraction: 1Mb vs 500kb bin size"
  ) +
  coord_equal(xlim = c(0, 0.65), ylim = c(0, 0.65)) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")

# Bland-Altman: plot the difference between the two estimates against their mean. Neither
# bin size is a gold standard, so the mean stands in for the unknown true value. A bias near
# zero means the two agree on average; the limits of agreement bound how far an individual
# sample can differ, which is the number that matters for a per-sample cohort decision.
df$mean_tf <- (df$TF_1Mb + df$TF_500kb) / 2
df$diff_tf <- df$TF_500kb - df$TF_1Mb

bias <- mean(df$diff_tf)
sd_diff <- sd(df$diff_tf)
upper_loa <- bias + 1.96 * sd_diff
lower_loa <- bias - 1.96 * sd_diff

bland_altman_plot <- ggplot(df, aes(x = mean_tf, y = diff_tf, color = flag)) +
  geom_hline(yintercept = 0, color = "grey40") +
  geom_hline(yintercept = bias, linetype = "dashed", color = "black") +
  geom_hline(yintercept = upper_loa, linetype = "dotted", color = "grey40") +
  geom_hline(yintercept = lower_loa, linetype = "dotted", color = "grey40") +
  geom_point(size = 2.5, alpha = 0.85) +
  scale_color_manual(values = flag_colors, name = "Category") +
  annotate("text", x = 0.5, y = upper_loa, vjust = -0.5, hjust = 1,
           label = paste0("+1.96 SD = ", round(upper_loa, 3)), size = 3.2) +
  annotate("text", x = 0.5, y = lower_loa, vjust = 1.5, hjust = 1,
           label = paste0("-1.96 SD = ", round(lower_loa, 3)), size = 3.2) +
  annotate("text", x = 0.5, y = bias, vjust = -0.5, hjust = 1,
           label = paste0("bias = ", round(bias, 3)), size = 3.2) +
  labs(
    x = "Mean tumour fraction (1Mb, 500kb)",
    y = "Difference (500kb minus 1Mb)",
    title = "Bland-Altman agreement plot"
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")

combined_plot <- scatter_plot + bland_altman_plot +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

ggsave(file.path(output_dir, "bin_size_comparison.pdf"), combined_plot, width = 11, height = 5.5, dpi = 300)
ggsave(file.path(output_dir, "bin_size_comparison.png"), combined_plot, width = 11, height = 5.5, dpi = 300)

cat("Pearson r:", round(r_value, 4), "\n")
cat("Mean bias:", round(bias, 4), "\n")
cat("Limits of agreement:", round(lower_loa, 4), "to", round(upper_loa, 4), "\n")
cat("Plots saved to:", output_dir, "\n")
