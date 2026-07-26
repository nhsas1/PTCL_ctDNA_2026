#!/usr/bin/env Rscript

# Fraction of Genome Altered (FGA) and Aneuploidy Score (AS) from ichorCNA calls.
#
# FGA is a continuous measure of how much of the genome carries a copy-number change,
# weighted by segment length. AS is the count of whole chromosome arms called as gained or
# lost, following Taylor et al. 2018. The two capture different things: FGA is dominated by
# focal and broad events alike, whereas AS deliberately ignores focal events and measures
# only arm-scale aneuploidy. Reporting both distinguishes a genome with many small changes
# from one with whole-arm imbalance.
#
# Both are computed at two AS thresholds. The 50% threshold is the Taylor definition (an
# arm is called if more than half of it is altered). The 25% threshold is a more sensitive
# variant included because at ~1x coverage and low tumour fraction, real arm-level events
# are frequently segmented into pieces that individually fall short of 50%.
#
# Inputs are ichorCNA seg files, preferring the manually curated versions where they
# exist. Outputs are FGA_AS_summary.csv (one row per sample) and
# arm_level_scores_detail.csv (one row per sample per arm).

library(dplyr)

cat("=== FGA and Aneuploidy Score Calculation ===\n\n")

# 1. Sample list and file paths
#
# This is the n=20 analysable CNA cohort, not the full 34. Samples below the ~3% tumour
# fraction detection floor are excluded because their copy-number calls are not
# distinguishable from noise, and FGA computed on noise is meaningless. Tumour fractions
# are carried here so the output table can be sorted and correlated without re-reading the
# ichorCNA params files.
batch12_dir  <- "/scratch/alice/n/nhsas1/PTCL/ichorCNA/output"
batch3_dir   <- "/scratch/alice/n/nhsas1/PTCL/ichorCNA/output_batch3"
corrected_dir <- "/scratch/alice/n/nhsas1/PTCL/ichorCNA/corrected_seg_files"

samples <- data.frame(
  sample = c(
    "B1_S01", "B1_S02", "B1_S05",
    "B2_S01", "B2_S02", "B2_S03", "B2_S04", "B2_S05", "B2_S06", "B2_S07",
    "B2_S08", "B2_S12", "B2_S15", "B2_S16",
    "B3_S01", "B3_S04", "B3_S05", "B3_S06", "B3_S08", "B3_S09"
  ),
  tumor_fraction = c(
    0.06348, 0.05354, 0.03841,
    0.1718, 0.04475, 0.08393, 0.4552, 0.5769, 0.5066, 0.2236,
    0.2838, 0.07153, 0.2478, 0.3281,
    0.03636, 0.04876, 0.06386, 0.6033, 0.6228, 0.3612
  ),
  correction_applied = c(
    "none", "none", "none",
    "none", "none", "none", "none", "none", "none", "none",
    "manual_curation", "none", "manual_curation", "none",
    "none", "none", "none", "none", "none", "manual_curation"
  ),
  stringsAsFactors = FALSE
)

# Three samples had ichorCNA solutions that required manual review and correction; those
# read from corrected_seg_files rather than the raw ichorCNA output. Curation altered only
# the corrected copy number and call columns (11 and 12), which are exactly the columns
# read below, so downstream results reflect the curated calls.
samples$seg_file <- sapply(samples$sample, function(s) {
  if (s == "B2_S08") {
    file.path(corrected_dir, "B2_S08.seg.corrected.txt")
  } else if (s == "B2_S15") {
    file.path(corrected_dir, "B2_S15.seg.corrected.txt")
  } else if (s == "B3_S09") {
    file.path(corrected_dir, "B3_S09.seg.corrected.txt")
  } else if (grepl("^B3_", s)) {
    file.path(batch3_dir, s, paste0(s, ".seg.txt"))
  } else {
    file.path(batch12_dir, s, paste0(s, ".seg.txt"))
  }
})

cat("Checking file paths...\n")
missing <- samples$seg_file[!file.exists(samples$seg_file)]
if (length(missing) > 0) {
  cat("ERROR - Missing files:\n")
  cat(paste(" ", missing, collapse = "\n"), "\n")
  stop("Cannot proceed with missing files.")
} else {
  cat("All 20 seg.txt files found.\n\n")
}

# 2. Chromosome arm coordinates (hg38)
#
# 39 arms, which is the Taylor et al. denominator: 22 autosomes give 44 arms, minus the
# five acrocentric p-arms (13p, 14p, 15p, 21p, 22p) which carry only repetitive satellite
# DNA and are not assayable. Sex chromosomes are absent throughout because ichorCNA was run
# with --chrs "c(1:22)", so no X or Y segments exist to score.
#
# IMPORTANT: these coordinates are aligned to ichorCNA's 1Mb bin grid, not to true hg38 arm
# boundaries. p-arms start at 1e6 rather than 1, and centromere positions are rounded to
# the bin edge. arm_length below is therefore the assayed span of the arm, not its
# reference length. This is internally consistent - the same coordinates define both the
# numerator (overlapping altered bases) and the denominator - so the fractions are valid,
# but the values should not be quoted as hg38 arm lengths.
arms <- data.frame(
  chr = c(
    1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10, 10,
    11, 11, 12, 12, 13, 14, 15, 16, 16, 17, 17, 18, 18, 19, 19, 20, 20, 21, 22
  ),
  arm = c(
    "1p", "1q", "2p", "2q", "3p", "3q", "4p", "4q", "5p", "5q",
    "6p", "6q", "7p", "7q", "8p", "8q", "9p", "9q", "10p", "10q",
    "11p", "11q", "12p", "12q", "13q", "14q", "15q",
    "16p", "16q", "17p", "17q", "18p", "18q", "19p", "19q", "20p", "20q",
    "21q", "22q"
  ),
  start = c(
    1e6, 123400000, 1e6, 93900000, 1e6, 90900000, 1e6, 50000000,
    1e6, 48800000, 1e6, 58500000, 1e6, 60100000, 1e6, 45200000,
    1e6, 43000000, 1e6, 39800000, 1e6, 53400000, 1e6, 35500000,
    20000000, 20000000, 25000000, 1e6, 36800000, 1e6, 25100000,
    1e6, 18500000, 1e6, 26200000, 1e6, 28100000, 14000000, 19000000
  ),
  end = c(
    123400000, 248956000, 93900000, 242193000, 90900000, 198295000,
    50000000, 190214000, 48800000, 181538000, 58500000, 170805000,
    60100000, 159345000, 45200000, 145138000, 43000000, 138394000,
    39800000, 133797000, 53400000, 135086000, 35500000, 133275000,
    113000000, 105000000, 101000000, 36800000, 90338000, 25100000,
    83257000, 18500000, 80373000, 26200000, 58617000, 28100000,
    64444000, 46709000, 50818000
  ),
  stringsAsFactors = FALSE
)

arms$arm_length <- arms$end - arms$start
cat(sprintf("Chromosome arms defined: %d arms (max AS = %d)\n\n", nrow(arms), nrow(arms)))

# 3. Read all seg files
cat("Reading seg.txt files...\n")

all_segs <- do.call(rbind, lapply(1:nrow(samples), function(i) {
  s <- samples$sample[i]
  f <- samples$seg_file[i]
  seg <- read.table(f, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  # Columns are taken by position from the ichorCNA seg format. Columns 7-8 are the raw
  # copy number and call; columns 11-12 are the purity- and ploidy-corrected versions.
  # Everything downstream uses the corrected copy number, because the raw call is relative
  # and does not account for the fact that a low tumour fraction sample dilutes every
  # change towards diploid. orig_CN is retained only for the curation check at the end.
  data.frame(
    sample     = s,
    chr        = seg[[2]],
    start      = seg[[3]],
    end        = seg[[4]],
    num_mark   = seg[[5]],
    logR       = seg[[6]],
    orig_CN    = seg[[7]],
    orig_call  = seg[[8]],
    corr_CN    = seg[[11]],
    corr_call  = seg[[12]],
    stringsAsFactors = FALSE
  )
}))

cat(sprintf("Loaded %d segments across %d samples.\n\n",
            nrow(all_segs), length(unique(all_segs$sample))))

# 4. Calculate FGA
cat("Calculating FGA...\n")

# FGA is length-weighted: each segment contributes its span in bases, so one 50Mb deletion
# counts far more than ten 100kb ones. The denominator is the total segmented span for that
# sample rather than a fixed genome size, which means FGA is the altered fraction of the
# assayed genome. Altered is defined strictly as corrected copy number not equal to 2, so
# copy-neutral loss of heterozygosity is not counted - sWGS cannot detect it without allele
# fractions.
fga_results <- all_segs %>%
  mutate(
    segment_size = end - start,
    is_altered   = corr_CN != 2,
    is_gain      = corr_CN > 2,
    is_loss      = corr_CN < 2
  ) %>%
  group_by(sample) %>%
  summarise(
    total_bases   = sum(segment_size),
    altered_bases = sum(segment_size[is_altered]),
    gain_bases    = sum(segment_size[is_gain]),
    loss_bases    = sum(segment_size[is_loss]),
    FGA           = altered_bases / total_bases,
    FGA_gain      = gain_bases / total_bases,
    FGA_loss      = loss_bases / total_bases,
    n_segments    = n(),
    .groups = "drop"
  )

cat("FGA calculated.\n\n")

# 5. Calculate arm-level aneuploidy score (both 50% and 25% thresholds)
cat("Calculating Aneuploidy Score...\n")

# Bases shared between a segment and an arm. Segments can straddle the centromere, so a
# single segment may contribute to both the p and q arm of a chromosome; taking the
# intersection rather than assigning the whole segment to one arm avoids double-counting.
# Returns 0 when the intervals do not meet.
calc_overlap <- function(seg_start, seg_end, arm_start, arm_end) {
  overlap_start <- max(seg_start, arm_start)
  overlap_end   <- min(seg_end, arm_end)
  max(0, overlap_end - overlap_start)
}

arm_scores_list <- list()

for (s in unique(all_segs$sample)) {
  sample_segs <- all_segs[all_segs$sample == s, ]

  for (j in 1:nrow(arms)) {
    arm_chr    <- arms$chr[j]
    arm_name   <- arms$arm[j]
    arm_start  <- arms$start[j]
    arm_end    <- arms$end[j]
    arm_length <- arms$arm_length[j]

    chr_segs <- sample_segs[sample_segs$chr == arm_chr, ]

    gain_bases <- 0
    loss_bases <- 0

    for (k in 1:nrow(chr_segs)) {
      overlap <- calc_overlap(chr_segs$start[k], chr_segs$end[k], arm_start, arm_end)
      if (overlap > 0) {
        cn <- chr_segs$corr_CN[k]
        if (cn > 2) gain_bases <- gain_bases + overlap
        if (cn < 2) loss_bases <- loss_bases + overlap
      }
    }

    # Fractions of the arm gained and lost. These are independent quantities and can both
    # be substantial on the same arm, when part of it is gained and another part lost.
    frac_gain <- gain_bases / arm_length
    frac_loss <- loss_bases / arm_length

    # 50% threshold, the Taylor et al. definition. Only one of the two fractions can exceed
    # 0.5, so the order of these branches cannot affect the outcome.
    if (frac_gain > 0.5) {
      status_50 <- "GAIN"
    } else if (frac_loss > 0.5) {
      status_50 <- "LOSS"
    } else {
      status_50 <- "NEUTRAL"
    }

    # 25% threshold, the sensitive variant. Unlike the 50% case, both fractions can exceed
    # 0.25 on the same arm, so the direction must be decided by which is larger rather than
    # by branch order. Previously this was an if/else-if chain testing gain first, which
    # returned GAIN whenever both fractions cleared the threshold - including arms that
    # were overwhelmingly deleted. Ties are broken towards loss, which is the conservative
    # choice for lymphoma where arm-level deletions predominate.
    if (frac_gain > 0.25 && frac_gain > frac_loss) {
      status_25 <- "GAIN"
    } else if (frac_loss > 0.25) {
      status_25 <- "LOSS"
    } else {
      status_25 <- "NEUTRAL"
    }

    arm_scores_list[[length(arm_scores_list) + 1]] <- data.frame(
      sample     = s,
      arm        = arm_name,
      chr        = arm_chr,
      frac_gain  = round(frac_gain, 4),
      frac_loss  = round(frac_loss, 4),
      status_50  = status_50,
      scored_50  = as.integer(status_50 != "NEUTRAL"),
      status_25  = status_25,
      scored_25  = as.integer(status_25 != "NEUTRAL"),
      stringsAsFactors = FALSE
    )
  }
}

arm_scores <- do.call(rbind, arm_scores_list)

# Summarise AS per sample at both thresholds
as_results <- arm_scores %>%
  group_by(sample) %>%
  summarise(
    AS_50           = sum(scored_50),
    arms_gained_50  = sum(status_50 == "GAIN"),
    arms_lost_50    = sum(status_50 == "LOSS"),
    AS_25           = sum(scored_25),
    arms_gained_25  = sum(status_25 == "GAIN"),
    arms_lost_25    = sum(status_25 == "LOSS"),
    .groups = "drop"
  )

cat("Aneuploidy Score calculated.\n\n")

# 6. Merge and output
final <- samples %>%
  select(sample, tumor_fraction, correction_applied) %>%
  left_join(fga_results, by = "sample") %>%
  left_join(as_results, by = "sample") %>%
  arrange(desc(tumor_fraction))

cat("=== RESULTS ===\n\n")
cat(sprintf("%-8s  TF%%    FGA    FGA_G  FGA_L  AS50 G50 L50 AS25 G25 L25 Corr\n", "Sample"))
cat(paste(rep("-", 88), collapse = ""), "\n")

for (i in 1:nrow(final)) {
  r <- final[i, ]
  cat(sprintf("%-8s  %5.1f  %5.3f  %5.3f  %5.3f  %3d  %2d  %2d  %3d  %2d  %2d  %s\n",
              r$sample,
              r$tumor_fraction * 100,
              r$FGA,
              r$FGA_gain,
              r$FGA_loss,
              r$AS_50,
              r$arms_gained_50,
              r$arms_lost_50,
              r$AS_25,
              r$arms_gained_25,
              r$arms_lost_25,
              r$correction_applied))
}

# 7. Save output files
output_dir <- "/scratch/alice/n/nhsas1/PTCL/ichorCNA/aneuploidy_scores"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

write.csv(final,
          file.path(output_dir, "FGA_AS_summary.csv"),
          row.names = FALSE)

write.csv(arm_scores,
          file.path(output_dir, "arm_level_scores_detail.csv"),
          row.names = FALSE)

cat(sprintf("\nOutput saved to: %s\n", output_dir))
cat("  FGA_AS_summary.csv          - one row per sample\n")
cat("  arm_level_scores_detail.csv - per-arm breakdown (both thresholds)\n")

# 8. Sanity checks
#
# These guard against the failure modes that would otherwise be invisible in a summary
# table: an FGA outside 0-1 or an AS above 39 would mean the arm coordinates or the
# overlap logic are wrong; AS_25 below AS_50 would be arithmetically impossible, since any
# arm passing 50% necessarily passes 25%; and FGA not equalling gain plus loss would mean
# a segment was counted in neither or both directions.
cat("\n=== SANITY CHECKS ===\n")

if (any(final$FGA < 0 | final$FGA > 1)) {
  cat("WARNING: FGA values outside 0-1 range!\n")
} else {
  cat("CHECK 1 PASS: All FGA values in 0-1 range\n")
}

if (any(final$AS_50 < 0 | final$AS_50 > 39)) {
  cat("WARNING: AS_50 values outside 0-39 range!\n")
} else {
  cat("CHECK 2a PASS: All AS_50 values in 0-39 range\n")
}

if (any(final$AS_25 < 0 | final$AS_25 > 39)) {
  cat("WARNING: AS_25 values outside 0-39 range!\n")
} else {
  cat("CHECK 2b PASS: All AS_25 values in 0-39 range\n")
}

if (any(final$AS_25 < final$AS_50)) {
  cat("WARNING: AS_25 < AS_50 for some samples (impossible)!\n")
} else {
  cat("CHECK 2c PASS: AS_25 >= AS_50 for all samples (expected)\n")
}

tf_fga_cor <- cor(final$tumor_fraction, final$FGA, method = "spearman")
cat(sprintf("CHECK 3 INFO: TF vs FGA Spearman correlation = %.3f\n", tf_fga_cor))

fga_check <- abs(final$FGA - (final$FGA_gain + final$FGA_loss))
if (any(fga_check > 0.001)) {
  cat("WARNING: FGA != FGA_gain + FGA_loss for some samples!\n")
} else {
  cat("CHECK 4 PASS: FGA = FGA_gain + FGA_loss for all samples\n")
}

# Confirms the curated samples actually differ from their raw calls, so a silent failure to
# read the corrected file would show up here as zero mismatches.
cat("\nCHECK 5: Corrected sample details:\n")
for (s in c("B2_S08", "B2_S15", "B3_S09")) {
  s_segs <- all_segs[all_segs$sample == s, ]
  mismatches <- sum(s_segs$orig_CN != s_segs$corr_CN)
  cat(sprintf("  %s: %d segments with orig_CN != corr_CN\n", s, mismatches))
}

cat("\n=== COMPLETE ===\n")
