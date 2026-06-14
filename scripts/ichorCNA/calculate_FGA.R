# =============================================================================
# FGA (Fraction of Genome Altered) Calculator
# PTCL ctDNA sWGS Thesis Analysis
# Input:  ichorCNA .seg.txt files (all 34 samples, 3 batches)
# Output: FGA_results_all34.csv + FGA_TF_correlation.pdf
# Column used: Corrected_Copy_Number (TF-corrected absolute CN)
# =============================================================================

library(dplyr)
library(readr)
library(ggplot2)

# =============================================================================
# 1. DEFINE PATHS
# =============================================================================

OUTPUT_DIR_B1B2   <- "/scratch/alice/n/nhsas1/PTCL/ichorCNA/output/"
OUTPUT_DIR_B3     <- "/scratch/alice/n/nhsas1/PTCL/ichorCNA/output_batch3/"
ICHORCNA_SUMMARY  <- "/scratch/alice/n/nhsas1/ptcl_ctdna_thesis/results/ichorCNA/ichorCNA_summary.csv"
RESULTS_DIR       <- "/scratch/alice/n/nhsas1/PTCL/ichorCNA/FGA_results/"

# Create results directory if it doesn't exist
dir.create(RESULTS_DIR, showWarnings = FALSE, recursive = TRUE)

# Samples to EXCLUDE (re-runs, not canonical)
EXCLUDE_SAMPLES <- c("B2_S05_500kb", "B2_S08_500kb", "B2_S08_ploidy234")

# =============================================================================
# 2. FIND ALL SEG.TXT FILES
# =============================================================================

seg_files_b1b2 <- list.files(
  path       = OUTPUT_DIR_B1B2,
  pattern    = "\\.seg\\.txt$",
  recursive  = TRUE,
  full.names = TRUE
)

seg_files_b3 <- list.files(
  path       = OUTPUT_DIR_B3,
  pattern    = "\\.seg\\.txt$",
  recursive  = TRUE,
  full.names = TRUE
)

# Combine both batches
all_seg_files <- c(seg_files_b1b2, seg_files_b3)

# Extract sample names from file paths
sample_names <- tools::file_path_sans_ext(
  tools::file_path_sans_ext(basename(all_seg_files))
)
# Note: double file_path_sans_ext removes both .seg and .txt

# Remove re-run samples
keep_idx     <- !(sample_names %in% EXCLUDE_SAMPLES)
all_seg_files <- all_seg_files[keep_idx]
sample_names  <- sample_names[keep_idx]

cat("=== Samples found ===\n")
cat(sprintf("Batch 1+2 files found: %d\n", length(seg_files_b1b2)))
cat(sprintf("Batch 3 files found:   %d\n", length(seg_files_b3)))
cat(sprintf("Re-runs excluded:      %d\n", sum(!keep_idx)))
cat(sprintf("Final sample count:    %d\n", length(all_seg_files)))
cat("\nSample list:\n")
print(sort(sample_names))

# Safety check
if (length(all_seg_files) != 34) {
  warning(sprintf(
    "Expected 34 samples but found %d. Please check paths and exclusion list.",
    length(all_seg_files)
  ))
}

# =============================================================================
# 3. READ AND CALCULATE FGA FOR EACH SAMPLE
# =============================================================================

calculate_fga <- function(seg_file, sample_id) {
  
  seg <- tryCatch(
    read_tsv(seg_file, col_types = cols(.default = "c"), show_col_types = FALSE),
    error = function(e) {
      warning(sprintf("Could not read file: %s\nError: %s", seg_file, e$message))
      return(NULL)
    }
  )
  
  if (is.null(seg)) return(NULL)
  
  # Convert numeric columns
  seg <- seg %>%
    mutate(
      start                 = as.numeric(start),
      end                   = as.numeric(end),
      Corrected_Copy_Number = as.numeric(Corrected_Copy_Number)
    ) %>%
    filter(!is.na(start), !is.na(end), !is.na(Corrected_Copy_Number))
  
  # Calculate segment size in bases
  seg <- seg %>%
    mutate(segment_size = end - start)
  
  # FGA using Corrected_Copy_Number
  # Diploid baseline = 2 copies
  # Gain  = Corrected_Copy_Number > 2
  # Loss  = Corrected_Copy_Number < 2
  # Neutral = Corrected_Copy_Number == 2
  
  total_bases   <- sum(seg$segment_size, na.rm = TRUE)
  gain_bases    <- sum(seg$segment_size[seg$Corrected_Copy_Number > 2],  na.rm = TRUE)
  loss_bases    <- sum(seg$segment_size[seg$Corrected_Copy_Number < 2],  na.rm = TRUE)
  altered_bases <- gain_bases + loss_bases
  
  # Count segments by call type
  n_gain <- sum(seg$Corrected_Call == "GAIN",  na.rm = TRUE)
  n_hetd <- sum(seg$Corrected_Call == "HETD",  na.rm = TRUE)
  n_neut <- sum(seg$Corrected_Call == "NEUT",  na.rm = TRUE)
  n_hlamp<- sum(seg$Corrected_Call == "HLAMP", na.rm = TRUE)
  n_total_segments <- nrow(seg)
  
  data.frame(
    sample           = sample_id,
    FGA              = round(altered_bases / total_bases, 4),
    FGA_gain         = round(gain_bases    / total_bases, 4),
    FGA_loss         = round(loss_bases    / total_bases, 4),
    total_bases_Mb   = round(total_bases / 1e6, 1),
    altered_bases_Mb = round(altered_bases / 1e6, 1),
    n_segments       = n_total_segments,
    n_gain_segments  = n_gain,
    n_loss_segments  = n_hetd,
    n_hlamp_segments = n_hlamp,
    n_neut_segments  = n_neut,
    stringsAsFactors = FALSE
  )
}

# Run for all samples
cat("\n=== Calculating FGA ===\n")
fga_list <- mapply(
  calculate_fga,
  seg_file  = all_seg_files,
  sample_id = sample_names,
  SIMPLIFY  = FALSE
)

fga_table <- bind_rows(fga_list) %>%
  arrange(sample)

cat(sprintf("FGA calculated for %d samples\n", nrow(fga_table)))

# =============================================================================
# 4. MERGE WITH ICHORCNA SUMMARY (TF, PLOIDY, MAD)
# =============================================================================

ichorcna <- read_csv(ICHORCNA_SUMMARY, show_col_types = FALSE)

# Standardise column name (handle both 'sample' and 'Sample')
names(ichorcna) <- tolower(names(ichorcna))

results <- left_join(fga_table, ichorcna, by = "sample")

# =============================================================================
# 5. CLASSIFY SAMPLES BY TF CONFIDENCE GROUP
# =============================================================================

results <- results %>%
  mutate(
    batch = case_when(
      grepl("^B1_", sample) ~ "Batch1",
      grepl("^B2_", sample) ~ "Batch2",
      grepl("^B3_", sample) ~ "Batch3"
    ),
    tf_group = case_when(
      tumor_fraction >= 0.10 ~ "High (≥10%)",
      tumor_fraction >= 0.03 ~ "Low-detectable (3–10%)",
      TRUE                   ~ "Below floor (<3%)"
    ),
    tf_group = factor(tf_group,
      levels = c("High (≥10%)", "Low-detectable (3–10%)", "Below floor (<3%)")
    ),
    FGA_reliable = tumor_fraction >= 0.03
  )

# =============================================================================
# 6. PRINT RESULTS TABLE
# =============================================================================

cat("\n=== FGA Results — All 34 Samples ===\n")
results_display <- results %>%
  select(sample, batch, tumor_fraction, ploidy, FGA, FGA_gain, FGA_loss,
         n_segments, tf_group) %>%
  arrange(desc(tumor_fraction))
print(as.data.frame(results_display))

cat("\n=== Summary Statistics ===\n")
cat(sprintf("Samples with FGA calculated: %d\n", nrow(results)))
cat(sprintf("Median FGA (all 34):         %.3f\n", median(results$FGA, na.rm=TRUE)))
cat(sprintf("Median FGA (TF ≥ 3%%):        %.3f\n",
    median(results$FGA[results$FGA_reliable], na.rm=TRUE)))
cat(sprintf("Range FGA:                   %.3f – %.3f\n",
    min(results$FGA, na.rm=TRUE), max(results$FGA, na.rm=TRUE)))

cat("\n=== Samples by TF Group ===\n")
print(results %>% count(tf_group))

# =============================================================================
# 7. SAVE CSV
# =============================================================================

write_csv(results, file.path(RESULTS_DIR, "FGA_results_all34.csv"))
cat(sprintf("\nResults saved to: %sFGA_results_all34.csv\n", RESULTS_DIR))

# =============================================================================
# 8. PLOT — TF vs FGA CORRELATION
# =============================================================================

p <- ggplot(results, aes(x = tumor_fraction * 100, y = FGA)) +
  
  # Shaded region for below-floor sample  annotate("rect",
    xmin = -Inf, xmax = 3, ymin = -Inf, ymax = Inf,
    fill = "grey90", alpha = 0.5
  ) +
  annotate("text",
    x = 1.5, y = max(results$FGA, na.rm=TRUE) * 0.95,
    label = "Below\ndetection\nfloor",
    size = 3, colour = "grey50", hjust = 0.5
  ) +
  
  # Points coloured by TF group
  geom_point(aes(colour = tf_group, shape = batch), size = 3.5, alpha = 0.85) +
  
  # Regression line for detectable samples only
  geom_smooth(
    data   = results %>% filter(FGA_reliable),
    method = "lm",
    colour = "black",
    linewidth = 0.7,
    linetype  = "dashed",
    se     = TRUE,
    alpha  = 0.15
  ) +
  
  # Label samples with TF > 10%
  ggrepel::geom_text_repel(
    data  = results %>% filter(tumor_fraction >= 0.10),
    aes(label = sample),
    size  = 2.8,
    max.overlaps = 15
  ) +
  
  scale_colour_manual(
    values = c(
      "High (≥10%)"            = "#2166AC",
      "Low-detectable (3–10%)" = "#F4A582",
      "Below floor (<3%)"      = "#D9D9D9"
    )
  ) +
  scale_shape_manual(values = c("Batch1" = 16, "Batch2" = 17, "Batch3" = 15)) +
  
  labs(
    title    = "Tumour Fraction vs Fraction of Genome Altered (FGA)",
    subtitle = "PTCL ctDNA cohort (n=34) | ichorCNA Corrected_Copy_Number",
    x        = "Tumour Fraction (%)",
    y        = "FGA (Fraction of Genome Altered)",
    colour   = "TF Group",
    shape    = "Batch",
    caption  = "Grey shading = below 3% detection floor | Dashed line = linear fit (detectable samples only)"
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position  = "right",
    plot.title       = element_text(face = "bold", size = 13),
    plot.subtitle    = element_text(size = 10, colour = "grey40"),
    plot.caption     = element_text(size = 8,  colour = "grey50")
  )

ggsave(
  filename = file.path(RESULTS_DIR, "FGA_TF_correlation.pdf"),
  plot     = p,
  width    = 10, height = 7
)

cat(sprintf("Plot saved to: %sFGA_TF_correlation.pdf\n", RESULTS_DIR))
cat("\n=== Done ===\n")
