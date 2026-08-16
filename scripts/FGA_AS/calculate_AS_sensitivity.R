library(dplyr)
library(readr)
library(ggplot2)
library(ggrepel)

results_dir <- "/scratch/alice/n/nhsas1/PTCL/ichorCNA/aneuploidy_scores/"

# Already have AS_50 and AS_25 from the corrected calculation
results <- read_csv(file.path(results_dir, "FGA_AS_summary.csv"),
                    show_col_types = FALSE)

comparison <- results %>%
  mutate(
    AS_difference = AS_25 - AS_50,
    direction = case_when(
      AS_difference > 2  ~ "Higher at 25%",
      AS_difference < -2 ~ "Lower at 25%",
      TRUE               ~ "Concordant"
    )
  ) %>%
  arrange(desc(tumor_fraction))

cat("=== AS at 50% vs 25% Threshold (n=20, post-correction) ===\n")
print(as.data.frame(comparison %>%
  select(sample, tumor_fraction, FGA, AS_50, AS_25, AS_difference, direction)))

cat("\n=== Concordance Summary ===\n")
cat(sprintf("Samples concordant (diff <= 2 arms): %d / %d\n",
    sum(comparison$direction == "Concordant"), nrow(comparison)))
cat(sprintf("Samples higher at 25%% threshold:    %d / %d\n",
    sum(comparison$direction == "Higher at 25%"), nrow(comparison)))
cat(sprintf("Samples lower at 25%% threshold:     %d / %d\n",
    sum(comparison$direction == "Lower at 25%"), nrow(comparison)))
cat(sprintf("\nSpearman AS_50 vs AS_25: rho = %.3f\n",
    cor(comparison$AS_50, comparison$AS_25, method = "spearman")))
cat(sprintf("Mean difference (25%% - 50%%):       %.1f arms\n",
    mean(comparison$AS_difference)))
cat(sprintf("Range of difference:                %d to %d arms\n",
    min(comparison$AS_difference), max(comparison$AS_difference)))

write_csv(comparison,
          file.path(results_dir, "AS_threshold_sensitivity_comparison.csv"))

# Plot: AS_50 vs AS_25
p_compare <- ggplot(comparison, aes(x = AS_50, y = AS_25)) +
  geom_abline(slope = 1, intercept = 0,
              linetype = "dashed", colour = "grey60", linewidth = 0.7) +
  geom_point(aes(colour = direction, size = tumor_fraction * 100), alpha = 0.85) +
  geom_text_repel(aes(label = sample), size = 2.8, max.overlaps = 20) +
  scale_colour_manual(values = c(
    "Concordant"    = "#3b6d11",
    "Higher at 25%" = "#854f0b",
    "Lower at 25%"  = "#a32d2d")) +
  scale_size_continuous(name = "TF (%)", range = c(2, 6)) +
  coord_equal() +
  labs(
    title    = "Sensitivity Analysis: AS at 50% vs 25% Arm Coverage Threshold",
    subtitle = "PTCL ctDNA cohort - detectable samples (TF >= 3%, n=20, post-correction)",
    x        = "Aneuploidy Score at 50% threshold",
    y        = "Aneuploidy Score at 25% threshold",
    colour   = "Agreement",
    caption  = "Dashed line = perfect concordance | Points above line = more arms detected at 25%") +
  theme_classic(base_size = 12) +
  theme(
    plot.title    = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(size = 10, colour = "grey40"),
    plot.caption  = element_text(size = 8,  colour = "grey50"))

ggsave(file.path(results_dir, "AS_sensitivity_comparison.pdf"),
       plot = p_compare, width = 9, height = 8)
ggsave(file.path(results_dir, "AS_sensitivity_comparison.png"),
       plot = p_compare, width = 9, height = 8, dpi = 300)

cat("\nPlots and comparison table saved. Done.\n")
