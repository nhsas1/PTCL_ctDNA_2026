# RASCAL Batch Analysis — Batch 3 (13 Samples)
# PTCL ctDNA Thesis | RASCAL v0.7.0 | R 4.3.1 | ALICE HPC
# Based on: run_rascal_batch.R (Batches 1 and 2)
# Changes from original:
#   - QDNASEQ_DIR points to QDNAseq_output_batch3
#   - OUTPUT_DIR points to RASCAL/output_batch3
#   - ichorcna_tf updated with Batch 3 TF values from ichorCNA
#   - all_samples updated to 13 Batch 3 sample IDs
# All other parameters identical to original pipeline.
#
# WHAT THIS IS FOR: RASCAL is an independent estimator of absolute copy number, fitting
# ploidy and cellularity by searching for the combination that places relative copy numbers
# closest to integers. It is run here as a cross-check on ichorCNA, not as a source of
# results - the thesis CNA calls come from ichorCNA. Nothing downstream reads the
# *_abs_cn.seg files this script writes.
#
# INTERPRETING THE OUTPUT - read this before quoting anything from the summary CSV:
#
# RASCAL needs enough tumour signal to find integer structure. Eleven of these thirteen
# samples sit below the 3% ichorCNA detection floor, where no such structure exists, so
# RASCAL fits whatever minimises distance on noise. That is why the summary reports
# cellularity of 0.80-0.94 for samples whose ichorCNA tumour fraction is under 2%, and why
# all thirteen come out DISCORDANT.
#
# Those rows are not evidence that the two methods disagree. They are evidence that RASCAL
# cannot fit a sample with no detectable tumour, which is expected and is the reason the CNA
# cohort is restricted to detectable samples in the first place. Only B3_S06, B3_S08 and
# B3_S09 have enough tumour content for the comparison to mean anything, and those three are
# re-run in run_rascal_constrained_batch3.R.

.libPaths(c("/home/n/nhsas1/R/library", .libPaths()))

library(rascal)
library(dplyr)
library(ggplot2)

# Global paths
QDNASEQ_DIR <- "/scratch/alice/n/nhsas1/PTCL/QDNAseq_output_batch3"
OUTPUT_DIR  <- "/scratch/alice/n/nhsas1/PTCL/RASCAL/output_batch3"
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ichorCNA tumour fraction reference — Batch 3
ichorcna_tf <- data.frame(
  sample = c("B3_S01","B3_S02","B3_S03","B3_S04","B3_S05",
             "B3_S06","B3_S07","B3_S08","B3_S09","B3_S10",
             "B3_S11","B3_S12","B3_S13"),
  ichorcna_tf = c(0.03636, 0.01195, 0.01944, 0.04876, 0.06386,
                  0.6033,  0.008588, 0.6228, 0.3612,  0.01431,
                  0.009693, 0.0903, 0.01993),
  stringsAsFactors = FALSE
)

# Convert a RASCAL (cellularity, ploidy) fit into the tumour fraction ichorCNA reports.
#
# These are different quantities and must not be compared directly. Cellularity is the
# proportion of CELLS that are tumour; ichorCNA's tumour fraction is the proportion of DNA
# derived from tumour. A tumour cell at ploidy p contributes p genome-equivalents of DNA
# while a normal cell contributes 2, so
#
#   tumour fraction = (cellularity * ploidy) / (cellularity * ploidy + (1 - cellularity) * 2)
#
# The two coincide only when ploidy is exactly 2. This is the same conversion used in
# run_rascal_constrained_batch3.R.
calc_tf <- function(cellularity, ploidy) {
  (cellularity * ploidy) / (cellularity * ploidy + (1 - cellularity) * 2)
}

# Sample list — all 13 Batch 3
all_samples <- c(
  "B3_S01","B3_S02","B3_S03","B3_S04","B3_S05",
  "B3_S06","B3_S07","B3_S08","B3_S09","B3_S10",
  "B3_S11","B3_S12","B3_S13"
)

cat("═══════════════════════════════════════════════════════\n")
cat("RASCAL BATCH — BATCH 3 (13 SAMPLES)\n")
cat("═══════════════════════════════════════════════════════\n")
cat("Samples to process:", length(all_samples), "\n")
cat("Output directory:", OUTPUT_DIR, "\n\n")

# Summary collector
summary_rows <- list()

# MAIN LOOP
for (SAMPLE_ID in all_samples) {

  cat("───────────────────────────────────────────────────────\n")
  cat("Processing:", SAMPLE_ID, "\n")
  cat("───────────────────────────────────────────────────────\n")

  sample_out <- file.path(OUTPUT_DIR, SAMPLE_ID)
  dir.create(sample_out, recursive = TRUE, showWarnings = FALSE)

  seg_file <- file.path(QDNASEQ_DIR, SAMPLE_ID,
                        paste0(SAMPLE_ID, ".seg"))
  igv_file <- file.path(QDNASEQ_DIR, SAMPLE_ID,
                        paste0(SAMPLE_ID, ".igv"))

  tf_row      <- ichorcna_tf[ichorcna_tf$sample == SAMPLE_ID, ]
  ICHORCNA_TF <- if (nrow(tf_row) > 0) tf_row$ichorcna_tf else NA

  cat("  ichorCNA TF:", ICHORCNA_TF, "\n")

  result <- tryCatch({

    cat("  Step 1: Reading seg file...\n")
    seg_data <- read.table(seg_file, header = TRUE,
                           sep = "\t", stringsAsFactors = FALSE)
    cat("    Segments:", nrow(seg_data), "\n")

    cat("  Step 2: Converting log2 to linear...\n")
    seg_data <- seg_data %>%
      mutate(
        chromosome  = as.character(CHROMOSOME),
        start       = as.numeric(START),
        end         = as.numeric(STOP),
        log2_ratio  = LOG2_RATIO_MEAN,
        relative_cn = 2 ^ log2_ratio
      ) %>%
      select(chromosome, start, end,
             log2_ratio, relative_cn, DATAPOINTS)

    cat("  Step 3: Filtering outliers...\n")
    n_before <- nrow(seg_data)
    seg_data <- seg_data %>%
      filter(log2_ratio > -3, log2_ratio < 3)
    n_after  <- nrow(seg_data)
    cat("    Removed:", n_before - n_after,
        "| Remaining:", n_after, "\n")

    cn_vector <- seg_data$relative_cn
    wt_vector <- seg_data$DATAPOINTS

    # Search ploidy-cellularity space for the combination placing relative copy numbers
    # closest to integers. The distance being minimised is a weighted median absolute
    # deviation from integer, with segments weighted by their bin count so long segments
    # count more than short noisy ones.
    #
    # CAVEAT on min_ploidy = 1.5: this floor is binding. Five of the thirteen samples
    # (B3_S06, B3_S07, B3_S08, B3_S10, B3_S12) return a fitted ploidy of exactly 1.50,
    # meaning the optimum lies at or below the edge of the search space and the true minimum
    # was never reached. A fit sitting on a grid boundary is not a converged fit and should
    # not be quoted as a ploidy estimate. RASCAL's own default lower bound is 1.25.
    #
    # Widening the floor would change every fit, so it is left as-is; see RERUN_REQUIRED.md.
    cat("  Step 4: Grid search...\n")
    grid <- absolute_copy_number_distance_grid(
      relative_copy_numbers = cn_vector,
      weights               = wt_vector,
      min_ploidy            = 1.5,
      max_ploidy            = 4.5,
      ploidy_step           = 0.05,
      min_cellularity       = 0.05,
      max_cellularity       = 1.00,
      cellularity_step      = 0.01
    )
    cat("    Best MAD in grid:", round(min(grid$distance), 4), "\n")

    # find_best_fit_solutions applies RASCAL's quality filters and can return several
    # acceptable solutions, or none. More than one means the data genuinely cannot
    # distinguish between e.g. a near-diploid and a doubled interpretation, which is the
    # same non-identifiability that excludes B3_S12 from the ichorCNA cohort.
    cat("  Step 5: Finding best solutions...\n")
    best <- find_best_fit_solutions(
      relative_copy_numbers = cn_vector,
      weights               = wt_vector,
      min_ploidy            = 1.5,
      max_ploidy            = 4.5,
      ploidy_step           = 0.05,
      min_cellularity       = 0.05,
      max_cellularity       = 1.00,
      cellularity_step      = 0.01
    )

    if (nrow(best) == 0) {
      cat("    WARNING: No solutions passed filters",
          "— using grid minimum\n")
      best <- grid %>%
        arrange(distance) %>%
        slice(1) %>%
        select(ploidy, cellularity, distance)
    }

    best_ploidy      <- best$ploidy[1]
    best_cellularity <- best$cellularity[1]
    best_distance    <- best$distance[1]
    n_solutions      <- nrow(best)

    cat(sprintf("    BEST: ploidy=%.3f | cellularity=%.3f | MAD=%.4f\n",
                best_ploidy, best_cellularity, best_distance))
    cat("    Solutions:", n_solutions,
        "| Ambiguous:", n_solutions > 1, "\n")

    rascal_expected_tf <- calc_tf(best_cellularity, best_ploidy)

    if (!is.na(ICHORCNA_TF)) {
      cat(sprintf("    Implied TF: %.3f | TF diff vs ichorCNA: %.3f\n",
                  rascal_expected_tf,
                  abs(rascal_expected_tf - ICHORCNA_TF)))
    }

    cat("  Step 6: Saving heatmap...\n")
    heatmap_file <- file.path(sample_out,
                              paste0(SAMPLE_ID, "_heatmap.pdf"))
    pdf(heatmap_file, width = 8, height = 6)
    p_heat <- distance_heatmap(
      distances            = grid,
      low_distance_colour  = "darkgreen",
      high_distance_colour = "white"
    ) + ggtitle(paste0(
        SAMPLE_ID,
        " | ploidy=",      round(best_ploidy, 2),
        " | cellularity=", round(best_cellularity, 2),
        " | MAD=",         round(best_distance, 4),
        " | n_solutions=", n_solutions,
        " | ichorCNA_TF=", ICHORCNA_TF))
    print(p_heat)
    dev.off()
    cat("    Saved:", heatmap_file, "\n")

    cat("  Step 7: Converting to absolute CN...\n")
    seg_data$absolute_cn <- relative_to_absolute_copy_number(
      relative_copy_numbers = cn_vector,
      ploidy                = best_ploidy,
      cellularity           = best_cellularity
    )

    cat("  Step 8: Building output segments...\n")
    output_segs <- seg_data %>%
      mutate(
        sample      = SAMPLE_ID,
        copy_number = absolute_cn
      ) %>%
      select(sample, chromosome, start, end,
             copy_number, log2_ratio, DATAPOINTS) %>%
      arrange(as.integer(chromosome), start)

    cat("  Step 9: Saving segment file...\n")
    seg_out <- file.path(sample_out,
                         paste0(SAMPLE_ID, "_abs_cn.seg"))
    write.table(output_segs, seg_out,
                sep = "\t", quote = FALSE, row.names = FALSE)
    cat("    Saved:", seg_out,
        "| Segments:", nrow(output_segs), "\n")

    cat("  Step 10: Saving genome plot...\n")
    igv_data <- read.table(igv_file, header = TRUE,
                           sep = "\t", stringsAsFactors = FALSE,
                           skip = 2, comment.char = "")

    igv_col  <- colnames(igv_data)[5]
    igv_data <- igv_data %>%
      rename(log2_ratio = all_of(igv_col)) %>%
      select(chromosome, start, end, log2_ratio) %>%
      mutate(chromosome = as.integer(chromosome)) %>%
      filter(!is.na(chromosome),
             chromosome %in% 1:22,
             is.finite(log2_ratio),
             log2_ratio > -3,
             log2_ratio < 3)

    seg_lines <- seg_data %>%
      mutate(chromosome = as.integer(chromosome)) %>%
      filter(!is.na(chromosome), chromosome %in% 1:22) %>%
      select(chromosome, start, end, log2_ratio)

    chrom_lengths <- igv_data %>%
      group_by(chromosome) %>%
      summarise(length = max(end), .groups = "drop") %>%
      arrange(chromosome) %>%
      mutate(offset = lag(cumsum(as.numeric(length)),
                          default = 0))

    igv_plot <- igv_data %>%
      left_join(chrom_lengths %>% select(chromosome, offset),
                by = "chromosome") %>%
      mutate(position = (start + end) / 2 + offset)

    seg_plot <- seg_lines %>%
      left_join(chrom_lengths %>% select(chromosome, offset),
                by = "chromosome") %>%
      mutate(start_pos = start + offset,
             end_pos   = end   + offset)

    chrom_mids     <- chrom_lengths %>%
      mutate(mid = offset + length / 2)
    chrom_dividers <- chrom_lengths$offset[-1]

    plot_file <- file.path(sample_out,
                           paste0(SAMPLE_ID, "_genome_plot.pdf"))
    pdf(plot_file, width = 14, height = 4)
    gp <- ggplot() +
      geom_vline(xintercept = chrom_dividers,
                 colour = "grey85", linewidth = 0.3) +
      geom_point(data = igv_plot,
                 aes(x = position, y = log2_ratio),
                 colour = "grey60", alpha = 0.25,
                 size = 0.2) +
      geom_hline(yintercept = 0, colour = "black",
                 linewidth = 0.3, linetype = "dashed") +
      geom_segment(data = seg_plot,
                   aes(x = start_pos, xend = end_pos,
                       y = log2_ratio, yend = log2_ratio),
                   colour = "red", linewidth = 1.0,
                   alpha = 0.9) +
      scale_x_continuous(
        breaks = chrom_mids$mid,
        labels = chrom_mids$chromosome,
        expand = expansion(mult = 0.01)
      ) +
      scale_y_continuous(
        limits = c(-2, 2),
        breaks = c(-2, -1, 0, 1, 2),
        expand = expansion(mult = 0)
      ) +
      labs(
        x     = "chromosome",
        y     = "log2 ratio",
        title = paste0(SAMPLE_ID,
                       " | ploidy=",      round(best_ploidy, 2),
                       " | cellularity=", round(best_cellularity, 2),
                       " | MAD=",         round(best_distance, 4),
                       " | ichorCNA_TF=", ICHORCNA_TF)
      ) +
      theme_bw() +
      theme(
        axis.title         = element_text(size = 11),
        axis.text.x        = element_text(size = 8),
        axis.text.y        = element_text(size = 10),
        axis.ticks.x       = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor   = element_blank(),
        plot.title         = element_text(size = 10)
      )
    print(gp)
    dev.off()
    cat("    Saved:", plot_file, "\n")

    list(
      sample             = SAMPLE_ID,
      rascal_ploidy      = round(best_ploidy, 3),
      rascal_cellularity = round(best_cellularity, 3),
      rascal_MAD         = round(best_distance, 4),
      n_solutions        = n_solutions,
      ambiguous          = n_solutions > 1,
      segments_used      = n_after,
      output_segments    = nrow(output_segs),
      ichorcna_tf        = ICHORCNA_TF,
      # Compare like with like: the tumour fraction implied by the RASCAL fit, not its
      # cellularity. Previously this subtracted cellularity from ichorCNA tumour fraction
      # directly, which is only valid at ploidy 2.
      rascal_expected_tf = round(rascal_expected_tf, 3),
      tf_difference      = round(abs(rascal_expected_tf - ICHORCNA_TF), 3),
      status             = "SUCCESS"
    )

  }, error = function(e) {
    cat("  ERROR:", conditionMessage(e), "\n")
    list(
      sample             = SAMPLE_ID,
      rascal_ploidy      = NA,
      rascal_cellularity = NA,
      rascal_MAD         = NA,
      n_solutions        = NA,
      ambiguous          = NA,
      segments_used      = NA,
      output_segments    = NA,
      ichorcna_tf        = ICHORCNA_TF,
      rascal_expected_tf = NA,
      tf_difference      = NA,
      status             = paste("FAILED:", conditionMessage(e))
    )
  })

  summary_rows[[SAMPLE_ID]] <- result
  cat("  Completed:", SAMPLE_ID, "— Status:", result$status, "\n\n")
}

# MASTER SUMMARY
cat("═══════════════════════════════════════════════════════\n")
cat("Building master summary...\n")

summary_df <- bind_rows(summary_rows)

summary_df <- summary_df %>%
  mutate(tf_agreement = case_when(
    is.na(rascal_cellularity)  ~ "UNRESOLVABLE",
    tf_difference <= 0.05      ~ "GOOD (<5%)",
    tf_difference <= 0.10      ~ "ACCEPTABLE (5-10%)",
    TRUE                       ~ "DISCORDANT (>10%)"
  ))

summary_file <- file.path(OUTPUT_DIR, "RASCAL_batch3_summary.csv")
write.csv(summary_df, summary_file, row.names = FALSE)
cat("Master summary saved:", summary_file, "\n\n")

cat("═══════════════════════════════════════════════════════\n")
cat("RASCAL BATCH 3 COMPLETE\n")
cat("═══════════════════════════════════════════════════════\n")
print(summary_df[, c("sample", "rascal_ploidy",
                     "rascal_cellularity", "rascal_MAD",
                     "ambiguous", "ichorcna_tf",
                     "tf_difference", "tf_agreement",
                     "status")])
cat("\n")
cat("Successful:               ",
    sum(summary_df$status == "SUCCESS"), "/ 13\n")
cat("Ambiguous (need curation):",
    sum(summary_df$ambiguous == TRUE, na.rm = TRUE), "\n")
cat("UNRESOLVABLE (low TF):    ",
    sum(summary_df$tf_agreement == "UNRESOLVABLE",
        na.rm = TRUE), "\n")
cat("Output directory:", OUTPUT_DIR, "\n")
cat("Summary CSV:", summary_file, "\n")
cat("═══════════════════════════════════════════════════════\n")
