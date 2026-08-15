#!/usr/bin/env Rscript
# GenVisR cnFreq and cnSpec — Cohort-Level CNA Visualisation
# 20 detectable PTCL ctDNA samples (TF >= 3%, post-correction)
# Uses Corrected_Copy_Number from ichorCNA seg.txt files

library(GenVisR)
library(ggplot2)

select <- dplyr::select
filter <- dplyr::filter

cat("=== GenVisR cnFreq and cnSpec ===\n\n")

# File paths
batch12_dir   <- "/scratch/alice/n/nhsas1/PTCL/ichorCNA/output"
batch3_dir    <- "/scratch/alice/n/nhsas1/PTCL/ichorCNA/output_batch3"
corrected_dir <- "/scratch/alice/n/nhsas1/PTCL/ichorCNA/corrected_seg_files"
output_dir    <- "/scratch/alice/n/nhsas1/PTCL/ichorCNA/aneuploidy_scores"

samples <- c(
  "B1_S01", "B1_S02", "B1_S05",
  "B2_S01", "B2_S02", "B2_S03", "B2_S04", "B2_S05", "B2_S06", "B2_S07",
  "B2_S08", "B2_S12", "B2_S15", "B2_S16",
  "B3_S01", "B3_S04", "B3_S05", "B3_S06", "B3_S08", "B3_S09"
)

corrected_samples <- c("B2_S08", "B2_S15", "B3_S09")

get_seg_path <- function(s) {
  if (s %in% corrected_samples) {
    file.path(corrected_dir, paste0(s, ".seg.corrected.txt"))
  } else if (grepl("^B3_", s)) {
    file.path(batch3_dir, s, paste0(s, ".seg.txt"))
  } else {
    file.path(batch12_dir, s, paste0(s, ".seg.txt"))
  }
}

# Read and combine seg files
cat("Reading seg files...\n")

all_segs <- do.call(rbind, lapply(samples, function(s) {
  f <- get_seg_path(s)
  df <- read.table(f, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  data.frame(
    sample  = s,
    chrom   = df[[2]],
    start   = df[[3]],
    end     = df[[4]],
    corr_CN = df[[11]],
    stringsAsFactors = FALSE
  )
}))

# Format for GenVisR: chromosome, start, end, segmean, sample
cnData <- data.frame(
  chromosome = paste0("chr", all_segs$chrom),
  start      = all_segs$start,
  end        = all_segs$end,
  segmean    = all_segs$corr_CN,
  sample     = all_segs$sample,
  stringsAsFactors = FALSE
)

# Keep autosomes only
cnData <- cnData[cnData$chromosome %in% paste0("chr", 1:22), ]

cat(sprintf("Loaded %d segments across %d samples\n\n",
            nrow(cnData), length(unique(cnData$sample))))

# Sample ordering by tumour fraction (high to low)
tf_order <- c(
  "B3_S08", "B3_S06", "B2_S05", "B2_S06", "B2_S04",
  "B3_S09", "B2_S16", "B2_S08", "B2_S15", "B2_S07",
  "B2_S01", "B2_S03", "B2_S12", "B3_S05", "B1_S01",
  "B1_S02", "B3_S04", "B2_S02", "B1_S05", "B3_S01"
)

cnData$sample <- factor(cnData$sample, levels = tf_order)

# Mark curated samples with asterisk for cnSpec labels
cnData_spec <- cnData
curated_labels <- setNames(
  ifelse(tf_order %in% corrected_samples,
         paste0(tf_order, " *"),
         tf_order),
  tf_order
)
cnData_spec$sample <- factor(
  curated_labels[as.character(cnData_spec$sample)],
  levels = curated_labels[tf_order]
)

# cnFreq — Frequency plot
cat("Generating cnFreq plot...\n")

freq_layer <- list(
  theme(
    strip.text.x   = element_text(size = 8, face = "bold"),
    panel.spacing  = unit(0.15, "lines"),
    plot.title     = element_text(size = 14, face = "bold"),
    plot.subtitle  = element_text(size = 10, colour = "grey40"),
    axis.title.y   = element_text(size = 11),
    axis.text.y    = element_text(size = 9)
  )
)

pdf(file.path(output_dir, "cnFreq_cohort_n20_v2.pdf"), width = 20, height = 7)
cnFreq(cnData,
       genome         = "hg38",
       CN_low_cutoff  = 1.5,
       CN_high_cutoff = 2.5,
       plotType       = "proportion",
       plotChr        = paste0("chr", 1:22),
       plot_title     = "Copy Number Alteration Frequency Across PTCL Cohort (n = 20)",
       CN_Loss_colour = "#2166AC",
       CN_Gain_colour = "#B2182B",
       x_title_size   = 12,
       y_title_size   = 12,
       plotLayer      = freq_layer)
dev.off()

png(file.path(output_dir, "cnFreq_cohort_n20_v2.png"),
    width = 3000, height = 1050, res = 200)
cnFreq(cnData,
       genome         = "hg38",
       CN_low_cutoff  = 1.5,
       CN_high_cutoff = 2.5,
       plotType       = "proportion",
       plotChr        = paste0("chr", 1:22),
       plot_title     = "Copy Number Alteration Frequency Across PTCL Cohort (n = 20)",
       CN_Loss_colour = "#2166AC",
       CN_Gain_colour = "#B2182B",
       x_title_size   = 12,
       y_title_size   = 12,
       plotLayer      = freq_layer)
dev.off()

cat("cnFreq saved\n\n")

# cnSpec — Copy number heatmap with curated samples marked
cat("Generating cnSpec heatmap...\n")

# hg38 autosome boundaries
hg38_autosomes <- data.frame(
  chromosome = paste0("chr", 1:22),
  start = 1,
  end = c(248956422, 242193529, 198295559, 190214555, 181538259,
          170805979, 159345973, 145138636, 138394717, 133797422,
          135086622, 133275309, 114364328, 107043718, 101991189,
          90338345, 83257441, 80373285, 58617616, 64444167,
          46709983, 50818468),
  stringsAsFactors = FALSE
)

spec_layer <- list(
  theme(
    strip.text.x = element_text(size = 9, face = "bold"),
    axis.text.y  = element_text(size = 8),
    plot.title   = element_text(size = 14, face = "bold"),
    plot.caption = element_text(size = 8, colour = "grey50", hjust = 0),
    legend.text  = element_text(size = 9),
    legend.title = element_text(size = 10, face = "bold")
  ),
  labs(caption = "* = manually curated sample | Ordered by tumour fraction (top = highest)")
)

pdf(file.path(output_dir, "cnSpec_cohort_n20_v2.pdf"), width = 26, height = 12)
cnSpec(cnData_spec,
       y              = hg38_autosomes,
       plot_title     = "Copy Number Spectrum Across PTCL Cohort (n = 20)",
       CN_Loss_colour = "#2166AC",
       CN_Gain_colour = "#B2182B",
       facet_lab_size = 9,
       CNscale        = "absolute",
       plotLayer      = spec_layer)
dev.off()

png(file.path(output_dir, "cnSpec_cohort_n20_v2.png"),
    width = 3900, height = 1800, res = 200)
cnSpec(cnData_spec,
       y              = hg38_autosomes,
       plot_title     = "Copy Number Spectrum Across PTCL Cohort (n = 20)",
       CN_Loss_colour = "#2166AC",
       CN_Gain_colour = "#B2182B",
       facet_lab_size = 9,
       CNscale        = "absolute",
       plotLayer      = spec_layer)
dev.off()

cat("cnSpec saved\n\n")

cat("=== Output files ===\n")
cat(sprintf("  %s/cnFreq_cohort_n20_v2.pdf\n", output_dir))
cat(sprintf("  %s/cnFreq_cohort_n20_v2.png\n", output_dir))
cat(sprintf("  %s/cnSpec_cohort_n20_v2.pdf\n", output_dir))
cat(sprintf("  %s/cnSpec_cohort_n20_v2.png\n", output_dir))
cat("\nDone.\n")
