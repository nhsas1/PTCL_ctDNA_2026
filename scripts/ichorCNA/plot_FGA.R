# Tumour fraction against FGA for the full 34-sample cohort - SUPERSEDED by
# plot_FGA_corrected.R, which uses the curated n=20 analysable cohort.
#
# This version reads the earlier pre-curation FGA table and deliberately keeps the
# below-floor samples in the plot, greying out the region under 3% tumour fraction. Its
# purpose is to justify the detection floor visually: samples below it show FGA values that
# reflect noise rather than biology, which is the argument for excluding them. The linear
# fit is deliberately restricted to detectable samples so the below-floor points cannot
# influence the slope they are being used to argue against.
#
# Retained for that argument. It writes into FGA_results/, a different directory from the
# corrected version, so the two do not overwrite each other despite sharing a filename.

library(dplyr)
library(ggplot2)
library(ggrepel)
library(readr)

RESULTS_DIR <- "/scratch/alice/n/nhsas1/PTCL/ichorCNA/FGA_results/"

results <- read_csv(paste0(RESULTS_DIR, "FGA_results_all34.csv"),
                    show_col_types = FALSE)

# This first assignment is immediately overwritten by the recomputation below and has no
# effect. It is left in place only because the recomputation depends on nothing from it;
# the grouping that actually reaches the plot is derived from tumor_fraction directly.
results$tf_group <- factor(results$tf_group,
  levels = c("High (>=10%)", "Low-detectable (3-10%)", "Below floor (<3%)"))

# Recompute the groups from the numeric tumour fraction rather than trusting the strings in
# the CSV, which can be mangled by the >= characters on a round trip through read/write.
# The 3% and 10% cut points are the detection floor and the high-burden threshold.
# Fix factor levels from CSV (may have been saved with special chars)
results <- results %>%
  mutate(tf_group = case_when(
    tumor_fraction >= 0.10 ~ "High (>=10%)",
    tumor_fraction >= 0.03 ~ "Low-detectable (3-10%)",
    TRUE                   ~ "Below floor (<3%)"
  ),
  tf_group = factor(tf_group,
    levels = c("High (>=10%)", "Low-detectable (3-10%)", "Below floor (<3%)")))

p <- ggplot(results, aes(x = tumor_fraction * 100, y = FGA)) +
  annotate("rect",
    xmin = -Inf, xmax = 3,
    ymin = -Inf, ymax = Inf,
    fill = "grey85", alpha = 0.6) +
  annotate("text",
    x = 1.5, y = 0.95,
    label = "Below\ndetection\nfloor",
    size = 3, colour = "grey40", hjust = 0.5) +
  geom_point(aes(colour = tf_group, shape = batch),
    size = 3.5, alpha = 0.85) +
  geom_smooth(
    data = results %>% filter(tumor_fraction >= 0.03),
    method = "lm",
    colour = "black",
    linewidth = 0.7,
    linetype = "dashed",
    se = TRUE,
    alpha = 0.15) +
  geom_text_repel(
    data = results %>% filter(tumor_fraction >= 0.03),
    aes(label = sample),
    size = 2.8,
    max.overlaps = 20) +
  scale_colour_manual(values = c(
    "High (>=10%)"            = "#2166AC",
    "Low-detectable (3-10%)"  = "#F4A582",
    "Below floor (<3%)"       = "#CCCCCC")) +
  scale_shape_manual(values = c(
    "Batch1" = 16, "Batch2" = 17, "Batch3" = 15)) +
  labs(
    title    = "Tumour Fraction vs Fraction of Genome Altered (FGA)",
    subtitle = "PTCL ctDNA cohort (n=34) | ichorCNA Corrected_Copy_Number",
    x        = "Tumour Fraction (%)",
    y        = "FGA (Fraction of Genome Altered)",
    colour   = "TF Group",
    shape    = "Batch",
    caption  = "Grey shading = below 3% detection floor | Dashed line = linear fit (detectable samples only)") +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "right",
    plot.title      = element_text(face = "bold", size = 13),
    plot.subtitle   = element_text(size = 10, colour = "grey40"),
    plot.caption    = element_text(size = 8,  colour = "grey50"))

ggsave(
  filename = paste0(RESULTS_DIR, "FGA_TF_correlation.pdf"),
  plot = p, width = 10, height = 7)

cat("Plot saved successfully\n")
