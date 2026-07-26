# Build the GISTIC2 input file from ichorCNA segment calls.
#
# This script defines the n=20 analysable CNA cohort, so it is where the cohort size in the
# thesis comes from. The derivation is:
#
#   34 samples total
#   - 13 with tumour fraction below the 3% detection floor   -> 21 remain
#   -  1 excluded for ploidy non-identifiability (B3_S12)     -> 20 analysable
#
# Everything downstream that says n=20 traces back to these two filters.
#
# GISTIC2 asks a different question from FGA or the aneuploidy score. Those describe how
# altered an individual genome is; GISTIC identifies regions altered more often across the
# cohort than expected by chance, which is what distinguishes a candidate driver locus from
# a passenger event that happens to be large.

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

# The 3% floor is the tumour fraction below which ichorCNA's copy-number calls are not
# separable from noise. Including such samples would not add signal; it would add segments
# whose calls are arbitrary, which in a recurrence analysis dilutes real peaks and can
# manufacture spurious ones.
detectable <- fga_results %>%
  filter(tumor_fraction >= 0.03) %>%
  pull(sample)

# B3_S12 clears the floor at 0.0903 but is excluded separately: its ploidy is not
# identifiable from coverage alone, so its absolute copy numbers - and therefore every
# gain or loss call derived from them - are not trustworthy.
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

# Substitute the manually curated seg files where they exist. This happens AFTER
# sample_names is derived, which matters: the corrected files are named
# <sample>.seg.corrected.txt, so stripping two extensions from them would yield
# "<sample>.seg" rather than the sample ID. Deriving names from the raw paths first and
# swapping the paths afterwards keeps the identifiers correct.
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

# Fail loudly if any expected sample has no seg file. Without this the cohort silently
# shrinks: a missing seg file simply fails to match, the sample drops out of sample_names,
# and GISTIC runs on fewer samples than intended with nothing but the count above to show
# it. Since GISTIC's q-values depend on the number of samples, that would quietly change
# every significance call rather than producing an obvious error.
missing_segs <- setdiff(detectable, sample_names)
if (length(missing_segs) > 0) {
  stop("No seg file found for ", length(missing_segs), " expected sample(s): ",
       paste(sort(missing_segs), collapse = ", "),
       "\nGISTIC would otherwise run on a silently reduced cohort.")
}

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
      # GISTIC expects a log ratio, so absolute copy number is converted relative to a
      # diploid baseline: CN 2 becomes 0, CN 1 becomes -1, CN 3 becomes +0.585.
      #
      # This uses Corrected_Copy_Number, NOT the raw seg.median.logR. That is the important
      # choice in this script. The corrected value has been adjusted for tumour fraction and
      # ploidy, so a heterozygous loss reads as CN 1 regardless of how diluted the sample
      # is, whereas the raw log ratio for the same event would be close to zero in a low
      # tumour fraction sample and GISTIC would never see it.
      #
      # The CN 0 branch avoids log2(0) = -Inf. It does not currently fire because ichorCNA
      # ran with --includeHOMD False, so no homozygous deletion state exists, but it must
      # stay if that flag is ever changed.
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
