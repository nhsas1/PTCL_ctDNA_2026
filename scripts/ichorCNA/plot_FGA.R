library(dplyr)
library(ggplot2)
library(ggrepel)
library(readr)

RESULTS_DIR <- "/scratch/alice/n/nhsas1/PTCL/ichorCNA/FGA_results/"

results <- read_csv(paste0(RESULTS_DIR, "FGA_results_all34.csv"),
                    show_col_types = FALSE)

results$tf_group <- factor(results$tf_group,
  levels = c("High (>=10%)", "Low-detectable (3-10%)", "Below floor (<3%)"))

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
