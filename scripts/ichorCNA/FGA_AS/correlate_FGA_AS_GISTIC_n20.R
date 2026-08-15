library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)
library(ggrepel)

fga_as_file   <- "/scratch/alice/n/nhsas1/PTCL/ichorCNA/aneuploidy_scores/FGA_AS_summary.csv"
gistic_counts <- "/scratch/alice/n/nhsas1/PTCL/GISTIC/output/GISTIC_per_sample_counts.csv"
gistic_peaks  <- "/scratch/alice/n/nhsas1/PTCL/GISTIC/output/GISTIC_per_sample_alterations.csv"
results_dir   <- "/scratch/alice/n/nhsas1/PTCL/ichorCNA/aneuploidy_scores/"

fga_as  <- read_csv(fga_as_file,   show_col_types = FALSE)
gcounts <- read_csv(gistic_counts, show_col_types = FALSE)
gpeaks  <- read_csv(gistic_peaks,  show_col_types = FALSE)

# Confirm cohort sizes match before proceeding
cat(sprintf("FGA/AS cohort size: %d\n", nrow(fga_as)))
cat(sprintf("GISTIC samples with >=1 peak: %d\n", nrow(gcounts)))

# GISTIC output has no Amplification column in this run (zero amp peaks found)
# so GISTIC_amp is fixed at 0 for every sample rather than read from file
combined <- fga_as %>%
  select(sample, tumor_fraction, correction_applied, FGA, FGA_gain, FGA_loss,
         AS_50, arms_gained_50, arms_lost_50, AS_25, arms_gained_25, arms_lost_25) %>%
  left_join(gcounts, by = "sample") %>%
  rename(GISTIC_del = Deletion, GISTIC_total = Total) %>%
  mutate(
    GISTIC_del   = ifelse(is.na(GISTIC_del),   0L, GISTIC_del),
    GISTIC_total = ifelse(is.na(GISTIC_total), 0L, GISTIC_total),
    GISTIC_amp   = 0L
  )

cat("\n=== Combined dataset (n=20, GISTIC re-run, deletion peaks only) ===\n")
print(as.data.frame(combined %>%
  select(sample, tumor_fraction, FGA, AS_50, AS_25, GISTIC_del, GISTIC_total) %>%
  arrange(desc(FGA))))

cat("\n=== Spearman correlations ===\n")

cor_fga_gtotal  <- cor.test(combined$FGA,   combined$GISTIC_total, method = "spearman")
cor_as50_gtotal <- cor.test(combined$AS_50, combined$GISTIC_total, method = "spearman")
cor_as25_gtotal <- cor.test(combined$AS_25, combined$GISTIC_total, method = "spearman")

cat(sprintf("FGA  vs GISTIC total (=deletions): rho = %+.3f  p = %.4f\n",
    cor_fga_gtotal$estimate, cor_fga_gtotal$p.value))
cat(sprintf("AS50 vs GISTIC total (=deletions): rho = %+.3f  p = %.4f\n",
    cor_as50_gtotal$estimate, cor_as50_gtotal$p.value))
cat(sprintf("AS25 vs GISTIC total (=deletions): rho = %+.3f  p = %.4f\n",
    cor_as25_gtotal$estimate, cor_as25_gtotal$p.value))

cat("\nNote: no amplification-specific correlations computed, since this\n")
cat("GISTIC run found zero significant amplification peaks in the n=20 cohort.\n")

write_csv(combined, file.path(results_dir, "FGA_AS_GISTIC_integrated_n20.csv"))
cat(sprintf("\nIntegrated table saved: %sFGA_AS_GISTIC_integrated_n20.csv\n", results_dir))

# GISTIC profile simplifies to two categories since no amplification peaks exist
combined <- combined %>%
  mutate(
    gistic_profile = ifelse(GISTIC_del > 0, "Deletion peak(s) present", "No significant peaks"),
    gistic_profile = factor(gistic_profile,
                             levels = c("Deletion peak(s) present", "No significant peaks"))
  )

col_del  <- "#2471A3"
col_none <- "grey70"

p1 <- ggplot(combined, aes(x = FGA, y = GISTIC_total)) +
  geom_smooth(method = "lm", colour = "grey50", linewidth = 0.6,
              linetype = "dashed", se = TRUE, alpha = 0.12) +
  geom_point(aes(colour = gistic_profile, size = tumor_fraction * 100),
             alpha = 0.85) +
  geom_text_repel(aes(label = sample), size = 2.8,
                  max.overlaps = 20, colour = "grey30") +
  scale_colour_manual(values = c(
    "Deletion peak(s) present" = col_del,
    "No significant peaks"     = col_none)) +
  scale_size_continuous(name = "TF (%)", range = c(2, 6)) +
  annotate("text", x = 0.05, y = max(combined$GISTIC_total) - 0.2,
           label = sprintf("Spearman rho = %.3f\np = %.4f",
                           cor_fga_gtotal$estimate, cor_fga_gtotal$p.value),
           hjust = 0, size = 3.2, colour = "grey30") +
  labs(
    title    = "FGA vs GISTIC2 Deletion Burden (n=20, re-run cohort)",
    subtitle = "PTCL ctDNA cohort - post-correction, GISTIC found no significant amp peaks",
    x        = "FGA (Fraction of Genome Altered)",
    y        = "GISTIC significant deletion peaks per patient",
    colour   = "GISTIC profile",
    caption  = "Point size = tumour fraction | Dashed line = linear fit"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title      = element_text(face = "bold", size = 13),
    plot.subtitle   = element_text(size = 10, colour = "grey40"),
    plot.caption    = element_text(size = 8,  colour = "grey50"),
    legend.position = "right"
  )

p2 <- ggplot(combined, aes(x = AS_50, y = GISTIC_total)) +
  geom_smooth(method = "lm", colour = "grey50", linewidth = 0.6,
              linetype = "dashed", se = TRUE, alpha = 0.12) +
  geom_point(aes(colour = gistic_profile, size = tumor_fraction * 100),
             alpha = 0.85) +
  geom_text_repel(aes(label = sample), size = 2.8,
                  max.overlaps = 20, colour = "grey30") +
  scale_colour_manual(values = c(
    "Deletion peak(s) present" = col_del,
    "No significant peaks"     = col_none)) +
  scale_size_continuous(name = "TF (%)", range = c(2, 6)) +
  annotate("text", x = 3, y = max(combined$GISTIC_total) - 0.2,
           label = sprintf("Spearman rho = %.3f\np = %.4f",
                           cor_as50_gtotal$estimate, cor_as50_gtotal$p.value),
           hjust = 0, size = 3.2, colour = "grey30") +
  labs(
    title    = "Taylor AS (50%) vs GISTIC2 Deletion Burden (n=20, re-run cohort)",
    subtitle = "PTCL ctDNA cohort - post-correction, GISTIC found no significant amp peaks",
    x        = "Aneuploidy Score (altered chromosome arms at 50%, max = 39)",
    y        = "GISTIC significant deletion peaks per patient",
    colour   = "GISTIC profile",
    caption  = "Point size = tumour fraction | Dashed line = linear fit"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title      = element_text(face = "bold", size = 13),
    plot.subtitle   = element_text(size = 10, colour = "grey40"),
    plot.caption    = element_text(size = 8,  colour = "grey50"),
    legend.position = "right"
  )

p3 <- ggplot(combined, aes(x = AS_25, y = GISTIC_total)) +
  geom_smooth(method = "lm", colour = "grey50", linewidth = 0.6,
              linetype = "dashed", se = TRUE, alpha = 0.12) +
  geom_point(aes(colour = gistic_profile, size = tumor_fraction * 100),
             alpha = 0.85) +
  geom_text_repel(aes(label = sample), size = 2.8,
                  max.overlaps = 20, colour = "grey30") +
  scale_colour_manual(values = c(
    "Deletion peak(s) present" = col_del,
    "No significant peaks"     = col_none)) +
  scale_size_continuous(name = "TF (%)", range = c(2, 6)) +
  annotate("text", x = 3, y = max(combined$GISTIC_total) - 0.2,
           label = sprintf("Spearman rho = %.3f\np = %.4f",
                           cor_as25_gtotal$estimate, cor_as25_gtotal$p.value),
           hjust = 0, size = 3.2, colour = "grey30") +
  labs(
    title    = "Taylor AS (25%) vs GISTIC2 Deletion Burden (n=20, re-run cohort)",
    subtitle = "PTCL ctDNA cohort - post-correction, GISTIC found no significant amp peaks",
    x        = "Aneuploidy Score (altered chromosome arms at 25%, max = 39)",
    y        = "GISTIC significant deletion peaks per patient",
    colour   = "GISTIC profile",
    caption  = "Point size = tumour fraction | Dashed line = linear fit"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title      = element_text(face = "bold", size = 13),
    plot.subtitle   = element_text(size = 10, colour = "grey40"),
    plot.caption    = element_text(size = 8,  colour = "grey50"),
    legend.position = "right"
  )

ggsave(file.path(results_dir, "FGA_vs_GISTIC_correlation_n20.pdf"), plot = p1, width = 10, height = 7)
ggsave(file.path(results_dir, "FGA_vs_GISTIC_correlation_n20.png"), plot = p1, width = 10, height = 7, dpi = 300)
ggsave(file.path(results_dir, "AS50_vs_GISTIC_correlation_n20.pdf"), plot = p2, width = 10, height = 7)
ggsave(file.path(results_dir, "AS50_vs_GISTIC_correlation_n20.png"), plot = p2, width = 10, height = 7, dpi = 300)
ggsave(file.path(results_dir, "AS25_vs_GISTIC_correlation_n20.pdf"), plot = p3, width = 10, height = 7)
ggsave(file.path(results_dir, "AS25_vs_GISTIC_correlation_n20.png"), plot = p3, width = 10, height = 7, dpi = 300)

cat("\nPlots saved:\n")
cat("  FGA_vs_GISTIC_correlation_n20.png\n")
cat("  AS50_vs_GISTIC_correlation_n20.png\n")
cat("  AS25_vs_GISTIC_correlation_n20.png\n")
cat("\nDone.\n")
