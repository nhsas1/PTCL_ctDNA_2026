library(dplyr)
library(readr)

seg_dir_b1b2 <- "/scratch/alice/n/nhsas1/PTCL/ichorCNA/output/"
seg_dir_b3   <- "/scratch/alice/n/nhsas1/PTCL/ichorCNA/output_batch3/"
fga_file     <- "/scratch/alice/n/nhsas1/PTCL/ichorCNA/FGA_results/FGA_results_all34.csv"
output_file  <- "/scratch/alice/n/nhsas1/PTCL/GISTIC/input/PTCL_GISTIC_input.seg"

corrected_dir     <- "/scratch/alice/n/nhsas1/PTCL/ichorCNA/corrected_seg_files/"
corrected_samples <- c("B2_S08", "B2_S15", "B3_S09")
excluded_samples  <- c("B3_S12")

fga_results <- read_csv(fga_file, show_col_types = FALSE)

detectable <- fga_results %>%
  filter(tumor_fraction >= 0.03) %>%
  pull(sample)

detectable <- setdiff(detectable, excluded_samples)

cat(sprintf("Detectable samples to include: %d\n", length(detectable)))
cat("Samples:\n")
print(sort(detectable))

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

for (s in corrected_samples) {
  idx <- which(sample_names == s)
  if (length(idx) == 1) {
    corrected_path <- file.path(corrected_dir, paste0(s, ".seg.corrected.txt"))
    if (file.exists(corrected_path)) {
      all_seg_files[idx] <- corrected_path
      cat(sprintf("Using corrected seg file for %s\n", s))
    } else {
      warning(paste("Corrected file not found for", s, "-", corrected_path))
    }
  }
}

cat(sprintf("\nSeg files matched: %d\n", length(all_seg_files)))

read_and_convert <- function(filepath, sample_id) {
  df <- read_tsv(filepath, col_types = cols(.default = "c"), show_col_types = FALSE)
  df <- df %>%
    mutate(
      chrom  = as.character(chrom),
      start  = as.integer(start),
      end    = as.integer(end),
      cn     = as.numeric(Corrected_Copy_Number),
      num_mark = as.integer(num.mark)
    ) %>%
    filter(!is.na(start), !is.na(end), !is.na(cn))
  df <- df %>%
    mutate(
      seg_cn = case_when(
        cn == 0 ~ -5.0,
        TRUE    ~ log2(cn / 2)
      ),
      seg_cn = round(seg_cn, 4)
    )
  df %>%
    mutate(
      Sample    = sample_id,
      Chromosome = chrom,
      Start      = start,
      End        = end,
      Num_Probes = num_mark,
      Segment_Mean = seg_cn
    ) %>%
    select(Sample, Chromosome, Start, End, Num_Probes, Segment_Mean)
}

cat("\nConverting seg files to GISTIC format...\n")
seg_list <- mapply(read_and_convert,
                   filepath  = all_seg_files,
                   sample_id = sample_names,
                   SIMPLIFY  = FALSE)
gistic_input <- bind_rows(seg_list)

cat(sprintf("Total segments in GISTIC input: %d\n", nrow(gistic_input)))
cat(sprintf("Samples represented:            %d\n", n_distinct(gistic_input$Sample)))

cat("\nSegment_Mean distribution (should be centred near 0):\n")
print(summary(gistic_input$Segment_Mean))

cat("\nSegment_Mean value counts (copy number states):\n")
cn_check <- gistic_input %>%
  mutate(approx_cn = round(2 * 2^Segment_Mean)) %>%
  count(approx_cn) %>%
  arrange(approx_cn)
print(cn_check)

cat("\nSample segment counts:\n")
print(gistic_input %>% count(Sample) %>% arrange(Sample) %>% as.data.frame())

write_tsv(gistic_input, output_file, col_names = TRUE)
cat(sprintf("\nGISTIC input file saved to: %s\n", output_file))

cat("\nFirst 5 rows of output file:\n")
print(head(gistic_input, 5))

cat("\nVerifying output file:\n")
verify <- read_tsv(output_file, show_col_types = FALSE)
cat(sprintf("Rows: %d | Columns: %d | Samples: %d\n",
    nrow(verify), ncol(verify), n_distinct(verify$Sample)))
cat("Column names:", paste(colnames(verify), collapse = ", "), "\n")

cat("\nInput preparation complete.\n")
