# Turn GISTIC2's all_lesions file into per-sample tables: which samples carry which
# significant peak, how many peaks each sample has, and how often each peak recurs.
#
# GISTIC reports peaks at cohort level; this recovers the per-sample detail behind them,
# which is what allows a peak to be checked for being driven by one or two extreme samples
# rather than by genuine recurrence.

library(dplyr)
library(readr)
library(tidyr)

outdir         <- "/scratch/alice/n/nhsas1/PTCL/GISTIC/output/"
lesions_f      <- paste0(outdir, "all_lesions.conf_99.txt")
seg_input_file <- "/scratch/alice/n/nhsas1/PTCL/GISTIC/input/PTCL_GISTIC_input.seg"

# Take the cohort size from the GISTIC input rather than counting samples in the results.
# This is deliberate and load-bearing: the counts table built below contains only samples
# with at least one peak (15 of 20 here), so deriving n from it would divide by 15 and
# overstate every recurrence frequency by a third.
n_samples <- read_tsv(seg_input_file, show_col_types = FALSE) %>%
  summarise(n = n_distinct(Sample)) %>%
  pull(n)
cat(sprintf("Detected n = %d samples from GISTIC input file\n", n_samples))

# all_lesions.conf_99.txt contains each peak twice: once with binary calls (0 = not
# altered, 1 = past the low amplitude threshold, 2 = past the high one) and once with the
# actual copy number values. Only the binary block is wanted, so the duplicate rows are
# dropped by their label. Reading both would double every peak.
lines <- readLines(lesions_f)
lines <- lines[lines != ""]
binary_lines <- lines[!grepl("CN values|Actual Copy", lines)]
header_line  <- binary_lines[1]
data_lines   <- binary_lines[-1]
col_names    <- strsplit(header_line, "\t")[[1]]

# GISTIC pads rows with a trailing tab inconsistently, so rows can come back one field
# short or long. Forcing the length to match the header pads with NA rather than letting
# rbind recycle values into the wrong columns.
parse_row <- function(line) {
  parts <- strsplit(line, "\t")[[1]]
  length(parts) <- length(col_names)
  parts
}
binary_df <- as.data.frame(
  do.call(rbind, lapply(data_lines, parse_row)),
  stringsAsFactors = FALSE
)
colnames(binary_df) <- col_names

sample_cols <- col_names[!col_names %in% c(
  "Unique Name", "Descriptor", "Wide Peak Limits", "Peak Limits",
  "Region Limits", "q values",
  "Residual q values after removing segments shared with higher peaks",
  "Broad or Focal", "Amplitude Threshold")]

peak_info <- binary_df %>%
  select(`Unique Name`, Descriptor, `q values`) %>%
  mutate(
    peak_type  = ifelse(grepl("Amplification", `Unique Name`), "Amplification", "Deletion"),
    peak_label = trimws(Descriptor),
    q_val      = as.numeric(trimws(`q values`))
  )

cat("\nPeaks found in current results:\n")
print(as.data.frame(peak_info %>% select(peak_label, peak_type, q_val)))

long_df <- binary_df %>%
  select(all_of(c("Unique Name", sample_cols))) %>%
  pivot_longer(cols = all_of(sample_cols),
               names_to = "sample", values_to = "call") %>%
  mutate(call = as.integer(trimws(call))) %>%
  left_join(peak_info %>% select(`Unique Name`, peak_label, peak_type, q_val),
            by = "Unique Name") %>%
  filter(!is.na(call))

# Keep only samples actually carrying each peak. Call 2 means the segment passed GISTIC's
# high amplitude threshold, call 1 the low one; at ~1x coverage and modest tumour fraction
# most real events land at 1, so a peak carried mostly at moderate amplitude is expected
# rather than weak evidence.
per_sample <- long_df %>%
  filter(call > 0) %>%
  mutate(strength = ifelse(call == 2, "High amplitude", "Moderate amplitude")) %>%
  select(sample, peak_label, peak_type, strength, q_val) %>%
  arrange(sample, peak_type, peak_label)

cat("\n=== Per-sample GISTIC alterations ===\n")
print(as.data.frame(per_sample), row.names = FALSE)

# NOTE this table omits samples with no significant peak, because they have no rows to
# count. Five of the twenty are absent here. Anything joining against this file must use a
# left join and fill the gaps with zero, or it will silently drop the low-burden end of the
# cohort - see correlate_FGA_AS_GISTIC_n20.R, which does exactly that.
counts <- per_sample %>%
  count(sample, peak_type) %>%
  pivot_wider(names_from = peak_type, values_from = n, values_fill = 0) %>%
  mutate(Total = rowSums(across(-sample))) %>%
  arrange(desc(Total))

cat("\n=== Count per sample ===\n")
print(as.data.frame(counts), row.names = FALSE)

freq <- per_sample %>%
  count(peak_label, peak_type, q_val) %>%
  mutate(freq_pct = round(n / n_samples * 100, 1)) %>%
  arrange(peak_type, desc(n))

cat(sprintf("\n=== Peak frequency across %d samples ===\n", n_samples))
print(as.data.frame(freq), row.names = FALSE)

for (pl in unique(per_sample$peak_label)) {
  cat(sprintf("\n=== %s ===\n", pl))
  print(as.data.frame(per_sample %>%
    filter(peak_label == pl) %>% select(sample, strength)), row.names = FALSE)
}

write_csv(per_sample, paste0(outdir, "GISTIC_per_sample_alterations.csv"))
write_csv(counts,     paste0(outdir, "GISTIC_per_sample_counts.csv"))
write_csv(freq,       paste0(outdir, "GISTIC_peak_frequency.csv"))

cat("\nAll CSV files saved.\nDone.\n")
