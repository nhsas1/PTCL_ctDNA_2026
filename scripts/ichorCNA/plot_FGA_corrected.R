# Two figures from the curated n=20 cohort: tumour fraction against FGA, and FGA against
# the Taylor Aneuploidy Score. This is the version that produced the committed figures;
# plot_FGA.R is the earlier pre-curation n=34 equivalent.
#
# The TF-versus-FGA plot addresses an obvious objection. If FGA simply tracked tumour
# fraction, it would be measuring assay sensitivity rather than tumour biology, and a
# strong relationship here would undermine using FGA as a biological descriptor. Showing
# the relationship explicitly is what allows that to be argued rather than assumed.
#
# CAVEAT on the second figure, which needs stating if it is used to support a claim: FGA
# and AS_50 are not independent measures. Both are computed from the same segment calls in
# calculate_FGA_AS.R - FGA is the length-weighted altered fraction, AS is a count of arms
# passing a coverage threshold on those same segments. A genome with high FGA necessarily
# has more altered arms. Their correlation is therefore close to tautological and is best
# presented as an internal consistency check, not as an independent finding.

library(dplyr)
library(ggplot2)
library(ggrepel)
library(readr)

RESULTS_DIR <- "/scratch/alice/n/nhsas1/PTCL/ichorCNA/aneuploidy_scores/"

results <- read_csv(paste0(RESULTS_DIR, "FGA_AS_summary.csv"),
                    show_col_types = FALSE)

# Derive batch and TF group from sample names
results <- results %>%
  mutate(
    batch = case_when(
      grepl("^B1_", sample) ~ "Batch1",
      grepl("^B2_", sample) ~ "Batch2",
      grepl("^B3_", sample) ~ "Batch3"
    ),
    # Only two groups here, unlike the n=34 version: every sample in the curated cohort is
    # above the 3% floor by construction, so the below-floor category cannot occur.
    tf_group = case_when(
      tumor_fraction >= 0.10 ~ "High (>=10%)",
      TRUE                   ~ "Low-detectable (3-10%)"
    ),
    tf_group = factor(tf_group,
      levels = c("High (>=10%)", "Low-detectable (3-10%)"))
  )

# Plot 1: TF vs FGA
p1 <- ggplot(results, aes(x = tumor_fraction * 100, y = FGA)) +
  geom_smooth(method = "lm", colour = "black", linewidth = 0.7,
              linetype = "dashed", se = TRUE, alpha = 0.15) +
  geom_point(aes(colour = tf_group, shape = batch),
             size = 3.5, alpha = 0.85) +
  geom_text_repel(aes(label = sample), size = 2.8, max.overlaps = 20) +
  scale_colour_manual(values = c(
    "High (>=10%)"            = "#2166AC",
    "Low-detectable (3-10%)"  = "#F4A582")) +
  scale_shape_manual(values = c(
    "Batch1" = 16, "Batch2" = 17, "Batch3" = 15)) +
  labs(
    title    = "Tumour Fraction vs Fraction of Genome Altered (FGA)",
    subtitle = "PTCL ctDNA cohort - detectable samples (TF >= 3%, n=20, post-correction)",
    x        = "Tumour Fraction (%)",
    y        = "FGA (Fraction of Genome Altered)",
    colour   = "TF Group",
    shape    = "Batch",
    caption  = "Dashed line = linear fit | 3 samples with manual curation applied") +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "right",
    plot.title      = element_text(face = "bold", size = 13),
    plot.subtitle   = element_text(size = 10, colour = "grey40"),
    plot.caption    = element_text(size = 8,  colour = "grey50"))

# Plot 2: FGA vs AS_50 scatter (Figure 4)
# See the caveat in the header. This rho is high by construction because both axes derive
# from the same segment calls, so it confirms the two measures are computed consistently
# rather than demonstrating a relationship between distinct quantities.
cor_fga_as <- cor.test(results$FGA, results$AS_50, method = "spearman")

p2 <- ggplot(results, aes(x = FGA, y = AS_50)) +
  geom_smooth(method = "lm", colour = "grey50", linewidth = 0.6,
              linetype = "dashed", se = TRUE, alpha = 0.12) +
  geom_point(aes(colour = tf_group, shape = batch,
                 size = tumor_fraction * 100), alpha = 0.85) +
  geom_text_repel(aes(label = sample), size = 2.8,
                  max.overlaps = 20, colour = "grey30") +
  scale_colour_manual(values = c(
    "High (>=10%)"            = "#2166AC",
    "Low-detectable (3-10%)"  = "#F4A582")) +
  scale_shape_manual(values = c(
    "Batch1" = 16, "Batch2" = 17, "Batch3" = 15)) +
  scale_size_continuous(name = "TF (%)", range = c(2, 6)) +
  annotate("text", x = 0.10, y = max(results$AS_50) - 1,
           label = sprintf("Spearman rho = %.3f\np = %.1e",
                           cor_fga_as$estimate, cor_fga_as$p.value),
           hjust = 0, size = 3.2, colour = "grey30") +
  labs(
    title    = "FGA vs Aneuploidy Score (Taylor AS, 50% threshold)",
    subtitle = "PTCL ctDNA cohort - detectable samples (TF >= 3%, n=20, post-correction)",
    x        = "FGA (Fraction of Genome Altered)",
    y        = "Aneuploidy Score (max = 39)",
    colour   = "TF Group",
    shape    = "Batch",
    caption  = "Point size = tumour fraction | Dashed line = linear fit") +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "right",
    plot.title      = element_text(face = "bold", size = 13),
    plot.subtitle   = element_text(size = 10, colour = "grey40"),
    plot.caption    = element_text(size = 8,  colour = "grey50"))

ggsave(file.path(RESULTS_DIR, "FGA_TF_correlation.pdf"),
       plot = p1, width = 10, height = 7)
ggsave(file.path(RESULTS_DIR, "FGA_TF_correlation.png"),
       plot = p1, width = 10, height = 7, dpi = 300)
ggsave(file.path(RESULTS_DIR, "FGA_vs_AS_scatter.pdf"),
       plot = p2, width = 10, height = 7)
ggsave(file.path(RESULTS_DIR, "FGA_vs_AS_scatter.png"),
       plot = p2, width = 10, height = 7, dpi = 300)

cat("FGA plots saved (PDF + PNG).\n")
cat(sprintf("FGA vs AS Spearman: rho = %.3f, p = %.1e\n",
    cor_fga_as$estimate, cor_fga_as$p.value))
