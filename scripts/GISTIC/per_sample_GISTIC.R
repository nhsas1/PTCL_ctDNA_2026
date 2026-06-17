library(dplyr)
library(readr)
library(tidyr)

outdir    <- "/scratch/alice/n/nhsas1/PTCL/GISTIC/output/"
lesions_f <- paste0(outdir, "all_lesions.conf_99.txt")

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
  select(`Unique Name`, `q values`) %>%
  mutate(
    peak_label = case_when(
      grepl("Amplification Peak 1", `Unique Name`) ~ "1q32.1 Amp",
      grepl("Amplification Peak 2", `Unique Name`) ~ "8q24.3 Amp (MYC)",
      grepl("Amplification Peak 3", `Unique Name`) ~ "11q24.3 Amp (ETS1/FLI1)",
      grepl("Deletion Peak 1",      `Unique Name`) ~ "1q21.3 Del (broad)",
      grepl("Deletion Peak 2",      `Unique Name`) ~ "2q37.3 Del (PDCD1)",
      grepl("Deletion Peak 3",      `Unique Name`) ~ "4q35.1 Del (FAT1/CASP3/IRF2)",
      grepl("Deletion Peak 4",      `Unique Name`) ~ "7q11.22 Del",
      TRUE ~ `Unique Name`
    ),
    q_val = as.numeric(trimws(`q values`))
  )

long_df <- binary_df %>%
  select(all_of(c("Unique Name", sample_cols))) %>%
  pivot_longer(cols = all_of(sample_cols),
               names_to = "sample", values_to = "call") %>%
  mutate(call = as.integer(trimws(call))) %>%
  left_join(peak_info %>% select(`Unique Name`, peak_label, q_val),
            by = "Unique Name") %>%
  filter(!is.na(call))

per_sample <- long_df %>%
  filter(call > 0) %>%
  mutate(
    strength  = ifelse(call == 2, "High amplitude", "Moderate amplitude"),
    peak_type = ifelse(grepl("Amp", peak_label), "Amplification", "Deletion")
  ) %>%
  select(sample, peak_label, peak_type, strength, q_val) %>%
  arrange(sample, peak_type, peak_label)

cat("=== Per-sample GISTIC alterations ===\n")
print(as.data.frame(per_sample), row.names = FALSE)

counts <- per_sample %>%
  count(sample, peak_type) %>%
  pivot_wider(names_from = peak_type, values_from = n, values_fill = 0) %>%
  mutate(Total = Amplification + Deletion) %>%
  arrange(desc(Total))
cat("\n=== Count per sample ===\n")
print(as.data.frame(counts), row.names = FALSE)

freq <- per_sample %>%
  count(peak_label, peak_type, q_val) %>%
  mutate(freq_pct = round(n / 21 * 100, 1)) %>%
  arrange(peak_type, desc(n))
cat("\n=== Peak frequency across 21 samples ===\n")
print(as.data.frame(freq), row.names = FALSE)

cat("\n=== 4q35.1 deletion (FAT1/CASP3/IRF2) ===\n")
print(as.data.frame(per_sample %>%
  filter(grepl("4q35", peak_label)) %>% select(sample, strength)))

cat("\n=== 8q24.3 amplification (MYC) ===\n")
print(as.data.frame(per_sample %>%
  filter(grepl("8q24", peak_label)) %>% select(sample, strength)))

cat("\n=== 2q37.3 deletion (PDCD1) ===\n")
print(as.data.frame(per_sample %>%
  filter(grepl("2q37", peak_label)) %>% select(sample, strength)))

write_csv(per_sample, paste0(outdir, "GISTIC_per_sample_alterations.csv"))
write_csv(counts,     paste0(outdir, "GISTIC_per_sample_counts.csv"))
write_csv(freq,       paste0(outdir, "GISTIC_peak_frequency.csv"))
cat("\nAll CSV files saved.\nDone.\n")
