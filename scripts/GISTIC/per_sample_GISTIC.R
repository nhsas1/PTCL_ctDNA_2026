library(dplyr)
library(readr)
library(tidyr)

outdir         <- "/scratch/alice/n/nhsas1/PTCL/GISTIC/output/"
lesions_f      <- paste0(outdir, "all_lesions.conf_99.txt")
seg_input_file <- "/scratch/alice/n/nhsas1/PTCL/GISTIC/input/PTCL_GISTIC_input.seg"

n_samples <- read_tsv(seg_input_file, show_col_types = FALSE) %>%
  summarise(n = n_distinct(Sample)) %>%
  pull(n)
cat(sprintf("Detected n = %d samples from GISTIC input file\n", n_samples))

lines <- readLines(lesions_f)
lines <- lines[lines != ""]
binary_lines <- lines[!grepl("CN values|Actual Copy", lines)]
header_line  <- binary_lines[1]
data_lines   <- binary_lines[-1]
col_names    <- strsplit(header_line, "\t")[[1]]

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

per_sample <- long_df %>%
  filter(call > 0) %>%
  mutate(strength = ifelse(call == 2, "High amplitude", "Moderate amplitude")) %>%
  select(sample, peak_label, peak_type, strength, q_val) %>%
  arrange(sample, peak_type, peak_label)

cat("\n=== Per-sample GISTIC alterations ===\n")
print(as.data.frame(per_sample), row.names = FALSE)

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
