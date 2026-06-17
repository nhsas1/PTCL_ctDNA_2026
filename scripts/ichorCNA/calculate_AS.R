library(dplyr)
library(readr)
library(ggplot2)
library(ggrepel)

seg_dir_b1b2 <- "/scratch/alice/n/nhsas1/PTCL/ichorCNA/output/"
seg_dir_b3   <- "/scratch/alice/n/nhsas1/PTCL/ichorCNA/output_batch3/"
fga_file     <- "/scratch/alice/n/nhsas1/PTCL/ichorCNA/FGA_results/FGA_results_all34.csv"
results_dir  <- "/scratch/alice/n/nhsas1/PTCL/ichorCNA/FGA_results/"

dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

fga_results <- read_csv(fga_file, show_col_types = FALSE)

detectable <- fga_results %>%
  filter(tumor_fraction >= 0.03) %>%
  pull(sample)

cat(sprintf("Detectable samples for AS calculation: %d\n", length(detectable)))

hg38_arms <- data.frame(
  chrom = c(
    "1","1",
    "2","2",
    "3","3",
    "4","4",
    "5","5",
    "6","6",
    "7","7",
    "8","8",
    "9","9",
    "10","10",
    "11","11",
    "12","12",
    "16","16",
    "17","17",
    "18","18",
    "19","19",
    "20","20",
    "13",
    "14",
    "15",
    "21",
    "22"
  ),
  arm = c(
    "p","q",
    "p","q",
    "p","q",
    "p","q",
    "p","q",
    "p","q",
    "p","q",
    "p","q",
    "p","q",
    "p","q",
    "p","q",
    "p","q",
    "p","q",
    "p","q",
    "p","q",
    "p","q",
    "p","q",
    "q",
    "q",
    "q",
    "q",
    "q"
  ),
  arm_start = c(
    0, 123400000,
    0, 93900000,
    0, 90900000,
    0, 50000000,
    0, 48800000,
    0, 60600000,
    0, 59900000,
    0, 45600000,
    0, 49000000,
    0, 39800000,
    0, 53400000,
    0, 35500000,
    0, 36400000,
    0, 25100000,
    0, 18500000,
    0, 26200000,
    0, 28100000,
    17700000,
    17200000,
    19000000,
    13200000,
    15000000
  ),
  arm_end = c(
    123400000, 248956422,
    93900000,  242193529,
    90900000,  198295559,
    50000000,  190214555,
    48800000,  181538259,
    60600000,  170805979,
    59900000,  159345973,
    45600000,  145138636,
    49000000,  138394717,
    39800000,  133797422,
    53400000,  135086622,
    35500000,  133275309,
    36400000,  90338345,
    25100000,  83257441,
    18500000,  80373285,
    26200000,  58617616,
    28100000,  64444167,
    107043718,
    107043718,
    102429344,
    48129895,
    50818468
  ),
  stringsAsFactors = FALSE
)
hg38_arms <- hg38_arms %>%
  mutate(arm_length = arm_end - arm_start)

cat(sprintf("Chromosome arms defined: %d (should be 39)\n", nrow(hg38_arms)))

seg_files_b1b2 <- list.files(seg_dir_b1b2, pattern = "\\.seg\\.txt$",
                              recursive = TRUE, full.names = TRUE)
seg_files_b3   <- list.files(seg_dir_b3,    pattern = "\\.seg\\.txt$",
                              recursive = TRUE, full.names = TRUE)

all_seg_files <- c(seg_files_b1b2, seg_files_b3)
sample_names  <- tools::file_path_sans_ext(
                   tools::file_path_sans_ext(basename(all_seg_files)))

keep          <- sample_names %in% detectable
all_seg_files <- all_seg_files[keep]
sample_names  <- sample_names[keep]

cat(sprintf("Seg files loaded: %d\n", length(all_seg_files)))

read_seg <- function(filepath, sample_id) {
  df <- read_tsv(filepath, col_types = cols(.default = "c"), show_col_types = FALSE)
  df %>%
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

seg_list     <- mapply(read_seg, filepath = all_seg_files,
                       sample_id = sample_names, SIMPLIFY = FALSE)
seg_combined <- bind_rows(seg_list)

cat(sprintf("Total segments across detectable samples: %d\n", nrow(seg_combined)))

score_arm <- function(seg_sample, arm_row, sample_ploidy) {

  chr        <- arm_row$chrom
  arm_start  <- arm_row$arm_start
  arm_end    <- arm_row$arm_end
  arm_length <- arm_row$arm_length

  segs_on_arm <- seg_sample %>%
    filter(chrom == chr,
           end   > arm_start,
           start < arm_end) %>%
    mutate(
      overlap_start  = pmax(start, arm_start),
      overlap_end    = pmin(end,   arm_end),
      overlap_length = pmax(0, overlap_end - overlap_start)
    ) %>%
    filter(overlap_length > 0)

  if (nrow(segs_on_arm) == 0) return(0L)

  total_overlap <- sum(segs_on_arm$overlap_length)
  if (total_overlap == 0) return(0L)

  ord            <- order(segs_on_arm$cn)
  cn_ord         <- segs_on_arm$cn[ord]
  wt_ord         <- segs_on_arm$overlap_length[ord]
  cum_wt         <- cumsum(wt_ord)
  half           <- sum(wt_ord) / 2
  weighted_median_cn <- cn_ord[which(cum_wt >= half)[1]]

  rounded_ploidy <- round(sample_ploidy)

  if (weighted_median_cn > rounded_ploidy) return(1L)
  if (weighted_median_cn < rounded_ploidy) return(-1L)
  return(0L)
}
calculate_sample_AS <- function(sample_id, seg_combined, hg38_arms, ploidy_df) {

  seg_sample    <- seg_combined %>% filter(sample == sample_id)
  sample_ploidy <- ploidy_df %>%
                     filter(sample == sample_id) %>%
                     pull(ploidy)

  if (length(sample_ploidy) == 0) {
    warning(sprintf("No ploidy found for %s — using 2", sample_id))
    sample_ploidy <- 2
  }

  arm_scores <- apply(hg38_arms, 1, function(arm_row) {
    score_arm(seg_sample,
              list(chrom      = arm_row["chrom"],
                   arm_start  = as.numeric(arm_row["arm_start"]),
                   arm_end    = as.numeric(arm_row["arm_end"]),
                   arm_length = as.numeric(arm_row["arm_length"])),
              sample_ploidy)
  })

  arm_labels <- paste0(hg38_arms$chrom, "_", hg38_arms$arm)
  names(arm_scores) <- arm_labels

  result <- data.frame(
    sample  = sample_id,
    ploidy_used = sample_ploidy,
    AS      = sum(abs(arm_scores)),
    AS_gain = sum(arm_scores == 1L),
    AS_loss = sum(arm_scores == -1L),
    stringsAsFactors = FALSE
  )

  arm_df        <- as.data.frame(t(arm_scores))
  colnames(arm_df) <- arm_labels
  cbind(result, arm_df)
}

ploidy_df <- fga_results %>%
  filter(sample %in% detectable) %>%
  select(sample, ploidy)

cat("\nCalculating Taylor AS for all detectable samples...\n")

as_list <- lapply(detectable, function(s) {
  cat(sprintf("  Processing %s ...\n", s))
  calculate_sample_AS(s, seg_combined, hg38_arms, ploidy_df)
})

as_results <- bind_rows(as_list)

combined <- left_join(
  fga_results %>% filter(tumor_fraction >= 0.03),
  as_results %>% select(sample, AS, AS_gain, AS_loss),
  by = "sample"
) %>%
  arrange(desc(tumor_fraction))

cat("\n=== Combined FGA + AS Results (detectable samples, n=21) ===\n")
print(as.data.frame(combined %>%
  select(sample, tumor_fraction, ploidy, FGA, AS, AS_gain, AS_loss, tf_group)))

cat("\n=== Summary Statistics ===\n")
cat(sprintf("Median AS:             %.1f\n",  median(combined$AS, na.rm = TRUE)))
cat(sprintf("Range AS:              %d - %d\n",
    min(combined$AS, na.rm = TRUE), max(combined$AS, na.rm = TRUE)))
cat(sprintf("Spearman FGA vs AS:    rho = %.3f, p = %.4f\n",
    cor(combined$FGA, combined$AS, method = "spearman", use = "complete.obs"),
    cor.test(combined$FGA, combined$AS, method = "spearman")$p.value))

write_csv(combined, file.path(results_dir, "FGA_AS_combined_detectable.csv"))

arm_cols <- grep("^[0-9]+_[pq]$", colnames(as_results), value = TRUE)
arm_calls <- as_results %>% select(sample, all_of(arm_cols))
write_csv(arm_calls, file.path(results_dir, "AS_per_arm_calls.csv"))

cat(sprintf("\nCSV results saved to: %s\n", results_dir))

p_scatter <- ggplot(combined, aes(x = AS, y = FGA)) +
  geom_point(aes(colour = tumor_fraction * 100, shape = batch), size = 3.5, alpha = 0.9) +
  geom_smooth(method = "lm", colour = "grey40", linewidth = 0.7,
              linetype = "dashed", se = TRUE, alpha = 0.15) +
  geom_text_repel(aes(label = sample), size = 2.8, max.overlaps = 20) +
  scale_colour_gradient(low = "#f7d488", high = "#1a4f8a", name = "TF (%)") +
  scale_shape_manual(values = c("Batch1" = 16, "Batch2" = 17, "Batch3" = 15)) +
  labs(
    title    = "FGA vs Arm-Level Aneuploidy Score (Taylor AS)",
    subtitle = "PTCL ctDNA cohort — detectable samples (TF >= 3%, n=21)",
    x        = "Aneuploidy Score (altered chromosome arms, max = 39)",
    y        = "FGA (Fraction of Genome Altered)",
    shape    = "Batch",
    caption  = "AS: weighted median CN per arm vs ichorCNA ploidy | hg38 | Corrected_Copy_Number"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title    = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(size = 10, colour = "grey40"),
    plot.caption  = element_text(size = 8,  colour = "grey50")
  )

ggsave(file.path(results_dir, "FGA_vs_AS_scatter.pdf"),
       plot = p_scatter, width = 10, height = 7)

p_bar <- combined %>%
  arrange(desc(AS)) %>%
  mutate(sample = factor(sample, levels = sample)) %>%
  ggplot(aes(x = sample, y = AS)) +
  geom_col(aes(fill = tumor_fraction * 100), width = 0.7) +
  geom_text(aes(label = AS), vjust = -0.4, size = 3.2) +
  scale_fill_gradient(low = "#f7d488", high = "#1a4f8a", name = "TF (%)") +
  labs(
    title    = "Arm-Level Aneuploidy Score per Sample",
    subtitle = "Detectable samples only (TF >= 3%, n=21) | Maximum possible = 39",
    x        = NULL,
    y        = "Aneuploidy Score"
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x   = element_text(angle = 45, hjust = 1, size = 9),
    plot.title    = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(size = 10, colour = "grey40")
  )

ggsave(file.path(results_dir, "AS_barplot.pdf"),
       plot = p_bar, width = 12, height = 6)

cat("All plots saved. Done.\n")
