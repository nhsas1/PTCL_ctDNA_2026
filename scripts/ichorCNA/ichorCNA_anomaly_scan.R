
# Anomaly Types Detected:
#   A — GAIN call with negative seg.median.logR
#   B — HETD (loss) call with positive seg.median.logR
#   C — logR_Copy_Number and Corrected_Copy_Number strongly disagree (>1.0)
#   D — HLAMP call with logR below TF-adjusted expected minimum
#
# Key concept for Type D:
#   The expected logR for any copy number state is not fixed — it scales with
#   tumour fraction (TF) using the ichorCNA mixture model formula:
#     expected_logR = log2( (TF * cn / ploidy) + (1 - TF) )
#   A fixed threshold of e.g. 0.3 is WRONG for low-TF samples.
#   This script computes a per-sample threshold from actual TF and ploidy.
#
# Outputs:
#   1. ichorCNA_anomaly_report.csv    — full anomaly table (one row per flagged segment)
#   2. ichorCNA_anomaly_summary.csv   — per-sample count of each anomaly type
#   3. ichorCNA_expected_logR.csv     — expected logR values per CN state per sample
#   4. Console report                 — readable summary for supervisor meeting
#


library(dplyr)
library(readr)


# 1. CONFIGURATION — Edit these paths to match your ALICE setup

# --- Paths for seg.txt files (two locations: Batch 1+2, Batch 3)
SEG_DIR_B1B2  <- "/scratch/alice/n/nhsas1/PTCL/ichorCNA/output/"
SEG_DIR_B3    <- "/scratch/alice/n/nhsas1/PTCL/ichorCNA/output_batch3/"

# --- ichorCNA summary file (contains TF, ploidy, MAD per sample)
SUMMARY_FILE  <- "/scratch/alice/n/nhsas1/ptcl_ctdna_thesis/results/ichorCNA/ichorCNA_summary.csv"

# --- Output directory for anomaly report files
OUTPUT_DIR    <- "/scratch/alice/n/nhsas1/PTCL/ichorCNA/anomaly_report/"

# --- Samples to exclude (re-runs or non-canonical)
EXCLUDE_SAMPLES <- c("B2_S05_500kb", "B2_S08_500kb", "B2_S08_ploidy234")

# --- Anomaly detection parameters
# Type D: HLAMP threshold factor
#   A segment is flagged if: seg.median.logR < expected_logR_HLAMP * HLAMP_FACTOR
#   0.5 means "logR is less than half of what physics predicts for this TF"
#   Use 0.5 as a conservative flag — you can discuss tightening this with supervisor
HLAMP_FACTOR <- 0.5

# Type C: Copy number discordance threshold
#   Flagged if abs(logR_Copy_Number - Corrected_Copy_Number) > CN_DISCORD_THRESHOLD
CN_DISCORD_THRESHOLD <- 1.0


# 2. HELPER FUNCTION — Expected logR from ichorCNA mixture model

# The ichorCNA mixture model formula:
#   expected_logR = log2( (TF * cn / ploidy) + (1 - TF) )
#
# Where:
#   TF     = tumour fraction (0 to 1)
#   cn     = copy number being modelled (1, 2, 3, 4, 5)
#   ploidy = estimated tumour ploidy from ichorCNA params
#
# This formula gives the logR you would EXPECT to see in the seg.txt if the
# copy number call were correct. Comparing this to the actual seg.median.logR
# tells you whether the call is physically consistent with the data.

expected_logR <- function(TF, cn, ploidy) {
  log2((TF * cn / ploidy) + (1 - TF))
}


# 3. FIND ALL SEG.TXT FILES


cat("=================================================================\n")
cat("ichorCNA Anomaly Scanner\n")
cat("=================================================================\n\n")

seg_files_b1b2 <- list.files(SEG_DIR_B1B2,
                               pattern    = "\\.seg\\.txt$",
                               recursive  = TRUE,
                               full.names = TRUE)

seg_files_b3   <- list.files(SEG_DIR_B3,
                               pattern    = "\\.seg\\.txt$",
                               recursive  = TRUE,
                               full.names = TRUE)

all_seg_files  <- c(seg_files_b1b2, seg_files_b3)

# Extract sample names from filenames
sample_names <- tools::file_path_sans_ext(
  tools::file_path_sans_ext(basename(all_seg_files))  # removes .seg then .txt
)

# Remove re-run samples
keep_idx      <- !(sample_names %in% EXCLUDE_SAMPLES)
all_seg_files <- all_seg_files[keep_idx]
sample_names  <- sample_names[keep_idx]

cat(sprintf("Seg files found:   %d\n", length(all_seg_files)))
cat(sprintf("Samples to scan:   %d\n\n", length(all_seg_files)))


# 4. LOAD SUMMARY FILE (TF + PLOIDY + MAD PER SAMPLE)


summary_df <- read_csv(SUMMARY_FILE, show_col_types = FALSE)

cat("ichorCNA summary loaded:", nrow(summary_df), "samples\n\n")


# 5. COMPUTE EXPECTED LOGRFOR ALL CN STATES — TABLE FOR REFERENCE

cat("Computing expected logR per CN state per sample...\n\n")

expected_table <- summary_df %>%
  rowwise() %>%
  mutate(
    expected_logR_cn1  = expected_logR(tumor_fraction, 1, ploidy),
    expected_logR_cn2  = expected_logR(tumor_fraction, 2, ploidy),
    expected_logR_cn3  = expected_logR(tumor_fraction, 3, ploidy),
    expected_logR_cn4  = expected_logR(tumor_fraction, 4, ploidy),
    expected_logR_cn5  = expected_logR(tumor_fraction, 5, ploidy),
    hlamp_flag_threshold = expected_logR(tumor_fraction, 5, ploidy) * HLAMP_FACTOR
  ) %>%
  ungroup() %>%
  select(sample, tumor_fraction, ploidy, MAD,
         expected_logR_cn1, expected_logR_cn2, expected_logR_cn3,
         expected_logR_cn4, expected_logR_cn5, hlamp_flag_threshold)

cat("Expected logR table (showing all samples):\n")
cat(sprintf("%-10s | %6s | %6s | %6s | %8s | %8s | %8s | %8s | %8s | %10s\n",
            "Sample", "TF", "Ploidy", "MAD",
            "E[cn=1]", "E[cn=2]", "E[cn=3]", "E[cn=4]", "E[cn=5]", "HLAMP_flag"))
cat(strrep("-", 105), "\n")

for (i in 1:nrow(expected_table)) {
  r <- expected_table[i, ]
  cat(sprintf("%-10s | %6.4f | %6.3f | %6.5f | %8.4f | %8.4f | %8.4f | %8.4f | %8.4f | %10.4f\n",
              r$sample, r$tumor_fraction, r$ploidy, r$MAD,
              r$expected_logR_cn1, r$expected_logR_cn2, r$expected_logR_cn3,
              r$expected_logR_cn4, r$expected_logR_cn5, r$hlamp_flag_threshold))
}
cat("\n")


# 6. MAIN SCAN LOOP — LOAD EACH SEG.TXT AND FLAG ANOMALIES


cat("=================================================================\n")
cat("Scanning all samples for anomalies...\n")
cat("=================================================================\n\n")

all_anomalies <- list()
summary_rows  <- list()

for (i in seq_along(all_seg_files)) {

  seg_file   <- all_seg_files[i]
  sample_id  <- sample_names[i]

  cat(sprintf("[%02d/%02d] %s\n", i, length(all_seg_files), sample_id))

  # Get TF and ploidy for this sample
  tf_row <- filter(summary_df, sample == sample_id)

  if (nrow(tf_row) == 0) {
    cat("  WARNING: No TF/ploidy data found in summary — skipping\n")
    next
  }

  tf     <- tf_row$tumor_fraction
  ploidy <- tf_row$ploidy
  mad    <- tf_row$MAD

  # Compute per-sample expected logR values
  e_cn1  <- expected_logR(tf, 1, ploidy)
  e_cn3  <- expected_logR(tf, 3, ploidy)
  e_cn5  <- expected_logR(tf, 5, ploidy)
  hlamp_threshold <- e_cn5 * HLAMP_FACTOR

  # Load seg.txt
  seg <- tryCatch(
    read.table(seg_file, header = TRUE, sep = "\t",
               stringsAsFactors = FALSE, quote = ""),
    error = function(e) {
      cat(sprintf("  ERROR reading file: %s\n", conditionMessage(e)))
      return(NULL)
    }
  )

  if (is.null(seg) || nrow(seg) == 0) {
    cat("  WARNING: Empty or unreadable seg.txt\n")
    next
  }

  # Rename columns for safety — some ichorCNA versions use slightly different names
  # Expected: ID, chrom, start, end, num.mark, seg.median.logR, copy.number,
  #           call, subclone.status, logR_Copy_Number, Corrected_Copy_Number, Corrected_Call
  seg$sample_id <- sample_id
  seg$TF        <- tf
  seg$ploidy    <- ploidy
  seg$MAD       <- mad

  n_total <- nrow(seg)


  # ANOMALY A: GAIN call with NEGATIVE logR

  anomaly_A <- seg %>%
    filter(call %in% c("GAIN", "AMP", "HLAMP") & seg.median.logR < 0) %>%
    mutate(
      anomaly_type        = "A",
      anomaly_description = "GAIN/AMP/HLAMP called but logR is NEGATIVE",
      expected_logR_for_call = case_when(
        call == "GAIN"  ~ expected_logR(tf, 3, ploidy),
        call == "AMP"   ~ expected_logR(tf, 4, ploidy),
        call == "HLAMP" ~ expected_logR(tf, 5, ploidy),
        TRUE ~ NA_real_
      )
    )


  # ANOMALY B: HETD (loss) call with POSITIVE logR

  anomaly_B <- seg %>%
    filter(call == "HETD" & seg.median.logR > 0) %>%
    mutate(
      anomaly_type        = "B",
      anomaly_description = "HETD (loss) called but logR is POSITIVE",
      expected_logR_for_call = expected_logR(tf, 1, ploidy)
    )


  # ANOMALY C: logR_Copy_Number and Corrected_Copy_Number strongly disagree
  # logR_Copy_Number = raw conversion from logR (no TF correction)
  # Corrected_Copy_Number = TF-corrected final call
  # A large gap means the model has applied a large correction — suspicious

  anomaly_C <- seg %>%
    filter(abs(logR_Copy_Number - Corrected_Copy_Number) > CN_DISCORD_THRESHOLD) %>%
    mutate(
      anomaly_type        = "C",
      anomaly_description = sprintf(
        "logR_CN=%.2f vs Corrected_CN=%d — discordance=%.2f (threshold %.1f)",
        logR_Copy_Number, Corrected_Copy_Number,
        abs(logR_Copy_Number - Corrected_Copy_Number),
        CN_DISCORD_THRESHOLD
      ),
      expected_logR_for_call = expected_logR(tf, Corrected_Copy_Number, ploidy)
    )


  # ANOMALY D: HLAMP call with logR below TF-adjusted minimum
  # Physical minimum for HLAMP (cn=5) at this TF:
  #   expected = log2((TF*5/ploidy) + (1-TF))
  # Flag if observed logR < expected * HLAMP_FACTOR (default 0.5)
  # This is the KEY anomaly type requiring TF adjustment

  anomaly_D <- seg %>%
    filter(call == "HLAMP" & seg.median.logR < hlamp_threshold) %>%
    mutate(
      anomaly_type        = "D",
      anomaly_description = sprintf(
        "HLAMP called but logR=%.4f < TF-adjusted threshold=%.4f (expected for cn=5 at TF=%.3f: %.4f)",
        seg.median.logR, hlamp_threshold, tf, e_cn5
      ),
      expected_logR_for_call = e_cn5
    )


  # COMBINED ANOMALY TABLE FOR THIS SAMPLE

  sample_anomalies <- bind_rows(anomaly_A, anomaly_B, anomaly_C, anomaly_D)

  # Add severity rating
  if (nrow(sample_anomalies) > 0) {
    sample_anomalies <- sample_anomalies %>%
      mutate(
        severity = case_when(
          # Impossible calls — logR direction contradicts call direction
          anomaly_type == "A" & seg.median.logR < -0.1  ~ "HIGH",
          anomaly_type == "B" & seg.median.logR >  0.1  ~ "HIGH",
          anomaly_type == "D" & seg.median.logR < 0     ~ "HIGH",
          # Subclonal anomalies — model flagged uncertainty itself
          subclone.status == TRUE                         ~ "MEDIUM",
          # Other flagged anomalies
          TRUE                                            ~ "LOW"
        )
      )
  }

  # Report per sample
  cat(sprintf("  Total segments: %d\n", n_total))
  cat(sprintf("  Type A (GAIN + neg logR):     %d\n", nrow(anomaly_A)))
  cat(sprintf("  Type B (HETD + pos logR):     %d\n", nrow(anomaly_B)))
  cat(sprintf("  Type C (CN discordance >%.1f): %d\n", CN_DISCORD_THRESHOLD, nrow(anomaly_C)))
  cat(sprintf("  Type D (HLAMP below TF-adj):  %d  [threshold: %.4f, expected cn5: %.4f]\n",
              nrow(anomaly_D), hlamp_threshold, e_cn5))
  cat(sprintf("  Total anomalies flagged:      %d\n\n",
              nrow(sample_anomalies)))

  if (nrow(sample_anomalies) > 0) {
    all_anomalies[[sample_id]] <- sample_anomalies
  }

  summary_rows[[sample_id]] <- data.frame(
    sample            = sample_id,
    TF                = tf,
    ploidy            = ploidy,
    MAD               = mad,
    n_segments        = n_total,
    n_anomaly_A       = nrow(anomaly_A),
    n_anomaly_B       = nrow(anomaly_B),
    n_anomaly_C       = nrow(anomaly_C),
    n_anomaly_D       = nrow(anomaly_D),
    n_total_anomalies = nrow(sample_anomalies),
    expected_logR_cn5 = e_cn5,
    hlamp_threshold   = hlamp_threshold,
    stringsAsFactors  = FALSE
  )
}


# 7. COMBINE AND SAVE OUTPUTS


dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# --- Full anomaly report (all flagged segments)
if (length(all_anomalies) > 0) {
  full_anomaly_df <- bind_rows(all_anomalies)

  # Select key columns for readability
  output_cols <- c("sample_id", "TF", "ploidy", "MAD",
                   "chrom", "start", "end", "num.mark",
                   "seg.median.logR", "copy.number", "call",
                   "subclone.status", "logR_Copy_Number",
                   "Corrected_Copy_Number", "Corrected_Call",
                   "anomaly_type", "anomaly_description",
                   "expected_logR_for_call", "severity")

  output_cols_present <- intersect(output_cols, colnames(full_anomaly_df))
  full_anomaly_df     <- full_anomaly_df[, output_cols_present]

  anomaly_file <- file.path(OUTPUT_DIR, "ichorCNA_anomaly_report.csv")
  write.csv(full_anomaly_df, anomaly_file, row.names = FALSE)
  cat(sprintf("Anomaly report saved: %s\n", anomaly_file))
  cat(sprintf("Total anomalous segments across cohort: %d\n\n", nrow(full_anomaly_df)))
} else {
  cat("No anomalies detected across cohort.\n\n")
}

# --- Per-sample summary
summary_df_out <- bind_rows(summary_rows)
summary_file   <- file.path(OUTPUT_DIR, "ichorCNA_anomaly_summary.csv")
write.csv(summary_df_out, summary_file, row.names = FALSE)
cat(sprintf("Summary saved: %s\n", summary_file))

# --- Expected logR reference table
expected_file <- file.path(OUTPUT_DIR, "ichorCNA_expected_logR.csv")
write.csv(expected_table, expected_file, row.names = FALSE)
cat(sprintf("Expected logR table saved: %s\n\n", expected_file))


# 8. CONSOLE SUMMARY REPORT


cat("=================================================================\n")
cat("COHORT ANOMALY SUMMARY\n")
cat("=================================================================\n\n")

cat(sprintf("%-10s | %6s | %6s | %5s | %8s | %8s | %8s | %8s | %8s\n",
            "Sample", "TF", "Ploidy", "Segs",
            "TypeA", "TypeB", "TypeC", "TypeD", "TOTAL"))
cat(strrep("-", 90), "\n")

for (i in 1:nrow(summary_df_out)) {
  r <- summary_df_out[i, ]
  flag <- if (r$n_total_anomalies > 0) " ***" else ""
  cat(sprintf("%-10s | %6.4f | %6.3f | %5d | %8d | %8d | %8d | %8d | %8d%s\n",
              r$sample, r$TF, r$ploidy, r$n_segments,
              r$n_anomaly_A, r$n_anomaly_B, r$n_anomaly_C,
              r$n_anomaly_D, r$n_total_anomalies, flag))
}

cat("\n*** = samples with anomalies requiring review\n\n")

cat("=================================================================\n")
cat("HIGH SEVERITY ANOMALIES (immediate review required)\n")
cat("=================================================================\n\n")

if (length(all_anomalies) > 0 && "severity" %in% colnames(full_anomaly_df)) {
  high_severity <- full_anomaly_df %>% filter(severity == "HIGH")
  if (nrow(high_severity) > 0) {
    for (i in 1:nrow(high_severity)) {
      r <- high_severity[i, ]
      cat(sprintf("Sample: %s | Chr: %s | logR: %.4f | Call: %s | Expected: %.4f\n",
                  r$sample_id, r$chrom, r$seg.median.logR, r$call,
                  r$expected_logR_for_call))
      cat(sprintf("  → %s\n\n", r$anomaly_description))
    }
  } else {
    cat("No HIGH severity anomalies found.\n\n")
  }
}

cat("=================================================================\n")
cat("TYPE D — HLAMP TF-ADJUSTED THRESHOLD EXPLANATION (per sample)\n")
cat("=================================================================\n\n")

cat("Formula: expected_logR = log2( (TF * cn / ploidy) + (1 - TF) )\n")
cat(sprintf("HLAMP flag threshold = expected_logR_cn5 * %.1f\n\n", HLAMP_FACTOR))

cat(sprintf("%-10s | %6s | %6s | %10s | %12s\n",
            "Sample", "TF", "Ploidy", "E[logR cn=5]", "HLAMP flag at"))
cat(strrep("-", 55), "\n")

for (i in 1:nrow(summary_df_out)) {
  r <- summary_df_out[i, ]
  cat(sprintf("%-10s | %6.4f | %6.3f | %10.4f | %12.4f\n",
              r$sample, r$TF, r$ploidy,
              r$expected_logR_cn5, r$hlamp_threshold))
}

cat("\n")
cat("Note: Samples with TF < 0.05 have expected logR for cn=5 < 0.1\n")
cat("      Any HLAMP call in these samples deserves close scrutiny.\n\n")

cat("=================================================================\n")
cat("Script complete.\n")
cat(sprintf("Output directory: %s\n", OUTPUT_DIR))
cat("=================================================================\n")
