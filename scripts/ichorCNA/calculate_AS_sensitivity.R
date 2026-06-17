library(dplyr)
library(readr)
library(ggplot2)
library(ggrepel)

seg_dir_b1b2 <- "/scratch/alice/n/nhsas1/PTCL/ichorCNA/output/"
seg_dir_b3   <- "/scratch/alice/n/nhsas1/PTCL/ichorCNA/output_batch3/"
fga_file     <- "/scratch/alice/n/nhsas1/PTCL/ichorCNA/FGA_results/FGA_results_all34.csv"
as_primary   <- "/scratch/alice/n/nhsas1/PTCL/ichorCNA/FGA_results/FGA_AS_combined_detectable.csv"
results_dir  <- "/scratch/alice/n/nhsas1/PTCL/ichorCNA/FGA_results/"

fga_results <- read_csv(fga_file, show_col_types = FALSE)
primary_as  <- read_csv(as_primary, show_col_types = FALSE)

detectable <- fga_results %>%
  filter(tumor_fraction >= 0.03) %>%
  pull(sample)

hg38_arms <- data.frame(
  chrom = c(
    "1","1","2","2","3","3","4","4","5","5",
    "6","6","7","7","8","8","9","9","10","10",
    "11","11","12","12","16","16","17","17",
    "18","18","19","19","20","20",
    "13","14","15","21","22"
  ),
  arm = c(
    "p","q","p","q","p","q","p","q","p","q",
    "p","q","p","q","p","q","p","q","p","q",
    "p","q","p","q","p","q","p","q",
    "p","q","p","q","p","q",
    "q","q","q","q","q"
  ),
  arm_start = c(
    0,123400000, 0,93900000, 0,90900000, 0,50000000, 0,48800000,
    0,60600000, 0,59900000, 0,45600000, 0,49000000, 0,39800000,
    0,53400000, 0,35500000, 0,36400000, 0,25100000,
    0,18500000, 0,26200000, 0,28100000,
    17700000,17200000,19000000,13200000,15000000
  ),
  arm_end = c(
    123400000,248956422, 93900000,242193529, 90900000,198295559, 50000000,190214555, 48800000,181538259,
    60600000,170805979, 59900000,159345973, 45600000,145138636, 49000000,138394717, 39800000,133797422,
    53400000,135086622, 35500000,133275309, 36400000,90338345, 25100000,83257441,
    18500000,80373285, 26200000,58617616, 28100000,64444167,
    107043718,107043718,102429344,48129895,50818468
  ),
  stringsAsFactors = FALSE
) %>%
  mutate(arm_length = arm_end - arm_start)

seg_files_b1b2 <- list.files(seg_dir_b1b2, pattern = "\\.seg\\.txt$",
                              recursive = TRUE, full.names = TRUE)
seg_files_b3   <- list.files(seg_dir_b3, pattern = "\\.seg\\.txt$",
                              recursive = TRUE, full.names = TRUE)

all_seg_files <- c(seg_files_b1b2, seg_files_b3)
sample_names  <- tools::file_path_sans_ext(
                   tools::file_path_sans_ext(basename(all_seg_files)))

keep          <- sample_names %in% detectable
all_seg_files <- all_seg_files[keep]
sample_names  <- sample_names[keep]

read_seg <- function(filepath, sample_id) {
  read_tsv(filepath, col_types = cols(.default = "c"), show_col_types = FALSE) %>%
    mutate(
      chrom  = as.character(chrom),
      start  = as.numeric(start),
      end    = as.numeric(end),
      cn     = as.numeric(Corrected_Copy_Number),
      sample = sample_id
    ) %>%
    filter(!is.na(start), !is.na(end), !is.na(cn)) %>%
    select(sample, chrom, start, end, cn)
}

seg_combined <- bind_rows(mapply(read_seg,
                                 filepath  = all_seg_files,
                                 sample_id = sample_names,
                                 SIMPLIFY  = FALSE))

ploidy_df <- fga_results %>%
  filter(sample %in% detectable) %>%
  select(sample, ploidy)

score_arm_25pct <- function(seg_sample, arm_row, sample_ploidy) {

  chr        <- arm_row$chrom
  arm_start  <- arm_row$arm_start
  arm_end    <- arm_row$arm_end
  arm_length <- arm_row$arm_length
  rounded_ploidy <- round(sample_ploidy)

  segs_on_arm <- seg_sample %>%
    filter(chrom == chr, end > arm_start, start < arm_end) %>%
    mutate(
      overlap_start  = pmax(start, arm_start),
      overlap_end    = pmin(end,   arm_end),
      overlap_length = pmax(0, overlap_end - overlap_start)
    ) %>%
    filter(overlap_length > 0)

  if (nrow(segs_on_arm) == 0 || arm_length == 0) return(0L)

  gain_frac <- sum(segs_on_arm$overlap_length[segs_on_arm$cn > rounded_ploidy],
                   na.rm = TRUE) / arm_length
  loss_frac <- sum(segs_on_arm$overlap_length[segs_on_arm$cn < rounded_ploidy],
                   na.rm = TRUE) / arm_length

  if (gain_frac >= 0.25 && gain_frac > loss_frac) return(1L)
  if (loss_frac >= 0.25 && loss_frac > gain_frac) return(-1L)
  return(0L)
}

calculate_AS_25 <- function(sample_id, seg_combined, hg38_arms, ploidy_df) {

  seg_sample    <- seg_combined %>% filter(sample == sample_id)
  sample_ploidy <- ploidy_df %>% filter(sample == sample_id) %>% pull(ploidy)
  if (length(sample_ploidy) == 0) sample_ploidy <- 2

  arm_scores <- apply(hg38_arms, 1, function(arm_row) {
    score_arm_25pct(
      seg_sample,
      list(chrom      = arm_row["chrom"],
           arm_start  = as.numeric(arm_row["arm_start"]),
           arm_end    = as.numeric(arm_row["arm_end"]),
           arm_length = as.numeric(arm_row["arm_length"])),
      sample_ploidy
    )
  })

  data.frame(
    sample   = sample_id,
    AS_25pct = sum(abs(arm_scores)),
    AS_25_gain = sum(arm_scores == 1L),
    AS_25_loss = sum(arm_scores == -1L),
    stringsAsFactors = FALSE
  )
}

cat("Calculating AS at 25% threshold for all detectable samples...\n")

as25_list    <- lapply(detectable, calculate_AS_25,
                       seg_combined = seg_combined,
                       hg38_arms    = hg38_arms,
                       ploidy_df    = ploidy_df)
as25_results <- bind_rows(as25_list)

comparison <- primary_as %>%
  select(sample, tumor_fraction, ploidy, FGA, AS, AS_gain, AS_loss, tf_group) %>%
  left_join(as25_results, by = "sample") %>%
  mutate(
    AS_difference  = AS_25pct - AS,
    direction      = case_when(
      AS_difference > 2  ~ "Higher at 25%",
      AS_difference < -2 ~ "Lower at 25%",
      TRUE               ~ "Concordant"
    )
  ) %>%
  arrange(desc(tumor_fraction))

cat("\n=== Primary AS (weighted median) vs Sensitivity AS (25% threshold) ===\n")
print(as.data.frame(comparison %>%
  select(sample, tumor_fraction, FGA, AS, AS_25pct, AS_difference, direction)))

cat("\n=== Concordance Summary ===\n")
cat(sprintf("Samples concordant (diff <= 2 arms): %d / %d\n",
    sum(comparison$direction == "Concordant"), nrow(comparison)))
cat(sprintf("Samples higher at 25%% threshold:    %d / %d\n",
    sum(comparison$direction == "Higher at 25%"), nrow(comparison)))
cat(sprintf("Samples lower at 25%% threshold:     %d / %d\n",
    sum(comparison$direction == "Lower at 25%"), nrow(comparison)))
cat(sprintf("\nSpearman AS vs AS_25pct: rho = %.3f\n",
    cor(comparison$AS, comparison$AS_25pct, method = "spearman")))
cat(sprintf("Mean difference (25%% - primary):    %.1f arms\n",
    mean(comparison$AS_difference)))
cat(sprintf("Range of difference:                 %d to %d arms\n",
    min(comparison$AS_difference), max(comparison$AS_difference)))

write_csv(comparison,
          file.path(results_dir, "AS_threshold_sensitivity_comparison.csv"))
cat(sprintf("\nComparison table saved to: %s\n", results_dir))

p_compare <- ggplot(comparison, aes(x = AS, y = AS_25pct)) +
  geom_abline(slope = 1, intercept = 0,
              linetype = "dashed", colour = "grey60", linewidth = 0.7) +
  geom_point(aes(colour = direction, size = tumor_fraction * 100), alpha = 0.85) +
  geom_text_repel(aes(label = sample), size = 2.8, max.overlaps = 20) +
  scale_colour_manual(
    values = c(
      "Concordant"    = "#3b6d11",
      "Higher at 25%" = "#854f0b",
      "Lower at 25%"  = "#a32d2d"
    )
  ) +
  scale_size_continuous(name = "TF (%)", range = c(2, 6)) +
  labs(
    title    = "Sensitivity Analysis — AS at 50% vs 25% Arm Coverage Threshold",
    subtitle = "PTCL ctDNA cohort — detectable samples (TF >= 3%, n=21)",
    x        = "Primary AS (weighted median, ~50% threshold)",
    y        = "Sensitivity AS (25% coverage threshold)",
    colour   = "Agreement",
    caption  = "Dashed line = perfect concordance | Points above line = higher score at 25% threshold"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title    = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(size = 10, colour = "grey40"),
    plot.caption  = element_text(size = 8,  colour = "grey50")
  )

ggsave(file.path(results_dir, "AS_sensitivity_comparison.pdf"),
       plot = p_compare, width = 9, height = 8)

cat("Sensitivity comparison plot saved. Done.\n")
