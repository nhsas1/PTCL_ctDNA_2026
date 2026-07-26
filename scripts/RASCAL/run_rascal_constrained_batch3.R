# RASCAL Constrained Run — Batch 3 (3 Samples)
# PTCL ctDNA Thesis | RASCAL v0.7.0 | R 4.3.1 | ALICE HPC
# Samples: B3_S06, B3_S08, B3_S09
# Cellularity search restricted using ichorCNA TF ± 0.15
# Based on: run_rascal_constrained.R (Batches 1 and 2)
# Changes from original:
#   - QDNASEQ_DIR points to QDNAseq_output_batch3
#   - OUTPUT_DIR points to RASCAL/output_batch3_constrained
#   - sample_table updated with 3 Batch 3 samples only
#
# WHY THESE THREE SAMPLES: the operative criterion is tumour fraction above 0.1, which
# selects exactly B3_S06 (0.603), B3_S08 (0.623) and B3_S09 (0.361). These are the only
# Batch 3 samples with enough tumour content for a constrained fit to be meaningful.
#
# The header previously gave the criterion as "ploidy=1.5 boundary hit or high TF
# discordance", which does not describe this selection. Five samples hit the ploidy floor
# (B3_S06, B3_S07, B3_S08, B3_S10, B3_S12) and all thirteen were labelled discordant, so
# that criterion would have selected five or thirteen samples rather than three. B3_S09 is
# included here despite fitting at ploidy 2.05, well off the floor.
#
# The practical consequence: B3_S07, B3_S10 and B3_S12 hit the ploidy floor and were never
# re-run. That is defensible, since all three are below the detection floor and a
# constrained fit would be no more meaningful than an unconstrained one, but it should be
# stated as the reason rather than left implied.
#
# NOTE the ploidy search is NOT constrained here - only cellularity is. The unconstrained
# run's ploidy boundary artefact is therefore not directly addressed by this script; the
# fitted ploidies move off the floor as a side effect of restricting cellularity.

.libPaths(c("/home/n/nhsas1/R/library", .libPaths()))

library(rascal)
library(dplyr)
library(ggplot2)

# Paths
QDNASEQ_DIR <- "/scratch/alice/n/nhsas1/PTCL/QDNAseq_output_batch3"
OUTPUT_DIR  <- "/scratch/alice/n/nhsas1/PTCL/RASCAL/output_batch3_constrained"
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Sample table — 3 Batch 3 samples
# ichorCNA values from verified params files
# Unconstrained RASCAL values from Run 1 summary CSV
sample_table <- data.frame(
  sample       = c("B3_S06", "B3_S08", "B3_S09"),
  ichorcna_tf  = c(0.6033,   0.6228,   0.3612),
  ichorcna_ply = c(2.002,    2.005,    2.002),
  rascal_unc_ploidy = c(1.50, 1.50, 2.05),
  rascal_unc_cell   = c(0.94, 0.94, 0.74),
  rascal_unc_MAD    = c(0.0704, 0.0623, 0.0741),
  stringsAsFactors = FALSE
)

# Constraint width
CONSTRAINT <- 0.15   # search TF ± 0.15

cat("═══════════════════════════════════════════════════════\n")
cat("RASCAL CONSTRAINED RUN — BATCH 3 (3 SAMPLES)\n")
cat("Cellularity window: ichorCNA TF ±", CONSTRAINT, "\n")
cat("Output:", OUTPUT_DIR, "\n")
cat("═══════════════════════════════════════════════════════\n\n")

# Summary collector
summary_rows <- list()

# MAIN LOOP
for (i in 1:nrow(sample_table)) {

  SAMPLE_ID    <- sample_table$sample[i]
  ICHORCNA_TF  <- sample_table$ichorcna_tf[i]
  ICHORCNA_PLY <- sample_table$ichorcna_ply[i]

  # KNOWN LIMITATION - the window is built on the wrong variable.
  #
  # The +/-0.15 band is derived from ichorCNA's TUMOUR FRACTION but applied as bounds on
  # RASCAL's CELLULARITY. Those quantities differ whenever ploidy is not 2, by exactly the
  # conversion calc_tf defines further down this script. Constraining cellularity to
  # [TF-0.15, TF+0.15] therefore does not constrain the implied tumour fraction to that
  # band.
  #
  # The effect is visible in the committed output. B3_S09 searched cellularity in
  # [0.21, 0.51], settled at cellularity 0.451 with ploidy 3.10, and produced an implied
  # tumour fraction of 0.560 against ichorCNA's 0.361 - a difference of 0.199, from a
  # constraint labelled +/-0.15. B3_S06 and B3_S08 are unaffected in practice because they
  # fit near ploidy 2, where the two quantities nearly coincide.
  #
  # Doing this properly means back-solving the cellularity bounds through calc_tf at each
  # candidate ploidy, which changes the searched space and so changes the fits. Left
  # unchanged here because that is a new analysis rather than a correction; see
  # RERUN_REQUIRED.md.
  cell_min <- max(0.05, ICHORCNA_TF - CONSTRAINT)
  cell_max <- min(1.00, ICHORCNA_TF + CONSTRAINT)

  cat("───────────────────────────────────────────────────────\n")
  cat("Sample:", SAMPLE_ID, "\n")
  cat("  ichorCNA TF:", ICHORCNA_TF,
      "| ploidy:", ICHORCNA_PLY, "\n")
  cat("  Cellularity search range: [",
      round(cell_min, 3), ",", round(cell_max, 3), "]\n")
  cat("  Unconstrained was: ploidy=",
      sample_table$rascal_unc_ploidy[i],
      "cellularity=", sample_table$rascal_unc_cell[i], "\n\n")

  sample_out <- file.path(OUTPUT_DIR, SAMPLE_ID)
  dir.create(sample_out, recursive = TRUE, showWarnings = FALSE)

  seg_file <- file.path(QDNASEQ_DIR, SAMPLE_ID,
                        paste0(SAMPLE_ID, ".seg"))
  igv_file <- file.path(QDNASEQ_DIR, SAMPLE_ID,
                        paste0(SAMPLE_ID, ".igv"))

  result <- tryCatch({

    cat("  Step 1: Reading seg file...\n")
    seg_data <- read.table(seg_file, header = TRUE,
                           sep = "\t", stringsAsFactors = FALSE)

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

    n_before <- nrow(seg_data)
    seg_data <- seg_data %>%
      filter(log2_ratio > -3, log2_ratio < 3)
    n_after  <- nrow(seg_data)
    cat("    Segments: total=", n_before,
        "kept=", n_after,
        "removed=", n_before - n_after, "\n")

    cn_vector <- seg_data$relative_cn
    wt_vector <- seg_data$DATAPOINTS

    cat("  Step 2: Constrained grid search...\n")
    cat("    Ploidy: 1.5–4.5 (unchanged)\n")
    cat("    Cellularity:", round(cell_min, 3),
        "–", round(cell_max, 3), "\n")

    grid_con <- absolute_copy_number_distance_grid(
      relative_copy_numbers = cn_vector,
      weights               = wt_vector,
      min_ploidy            = 1.5,
      max_ploidy            = 4.5,
      ploidy_step           = 0.05,
      min_cellularity       = cell_min,
      max_cellularity       = cell_max,
      cellularity_step      = 0.01
    )
    cat("    Combinations tested:", nrow(grid_con), "\n")
    cat("    Best MAD in constrained grid:",
        round(min(grid_con$distance), 4), "\n")
    cat("    Top 5 constrained solutions:\n")
    print(head(arrange(grid_con, distance), 5))

    cat("  Step 3: Finding best constrained solution...\n")
    best_con <- find_best_fit_solutions(
      relative_copy_numbers = cn_vector,
      weights               = wt_vector,
      min_ploidy            = 1.5,
      max_ploidy            = 4.5,
      ploidy_step           = 0.05,
      min_cellularity       = cell_min,
      max_cellularity       = cell_max,
      cellularity_step      = 0.01
    )

    if (nrow(best_con) == 0) {
      cat("    WARNING: No solutions passed filters",
          "— using constrained grid minimum\n")
      best_con <- grid_con %>%
        arrange(distance) %>%
        slice(1) %>%
        select(ploidy, cellularity, distance)
    }

    best_ploidy      <- best_con$ploidy[1]
    best_cellularity <- best_con$cellularity[1]
    best_MAD         <- best_con$distance[1]
    n_solutions      <- nrow(best_con)

    cat(sprintf("    Constrained best: ploidy=%.3f cellularity=%.3f MAD=%.4f\n",
                best_ploidy, best_cellularity, best_MAD))

    calc_tf <- function(c, p) {
      (c * p) / (c * p + (1 - c) * 2)
    }

    tf_from_constrained   <- calc_tf(best_cellularity, best_ploidy)
    tf_from_unconstrained <- calc_tf(
      sample_table$rascal_unc_cell[i],
      sample_table$rascal_unc_ploidy[i]
    )

    diff_constrained   <- abs(tf_from_constrained - ICHORCNA_TF)
    diff_unconstrained <- abs(tf_from_unconstrained - ICHORCNA_TF)
    improvement        <- diff_unconstrained - diff_constrained

    cat(sprintf("    Expected TF (constrained):   %.3f | diff: %.3f\n",
                tf_from_constrained, diff_constrained))
    cat(sprintf("    Expected TF (unconstrained): %.3f | diff: %.3f\n",
                tf_from_unconstrained, diff_unconstrained))
    cat(sprintf("    Improvement: %.3f\n", improvement))

    cat("  Step 4: Saving constrained heatmap...\n")
    heatmap_file <- file.path(sample_out,
                              paste0(SAMPLE_ID,
                                     "_constrained_heatmap.pdf"))
    pdf(heatmap_file, width = 8, height = 6)
    p_heat <- distance_heatmap(
      distances            = grid_con,
      low_distance_colour  = "darkgreen",
      high_distance_colour = "white"
    ) + ggtitle(paste0(
        SAMPLE_ID, " CONSTRAINED",
        " | ploidy=",      round(best_ploidy, 2),
        " | cellularity=", round(best_cellularity, 2),
        " | MAD=",         round(best_MAD, 4),
        " | ichorCNA_TF=", ICHORCNA_TF,
        " | window=±",     CONSTRAINT))
    print(p_heat)
    dev.off()
    cat("    Saved:", heatmap_file, "\n")

    cat("  Step 5: Converting to absolute CN...\n")
    seg_data$absolute_cn <- relative_to_absolute_copy_number(
      relative_copy_numbers = cn_vector,
      ploidy                = best_ploidy,
      cellularity           = best_cellularity
    )

    output_segs <- seg_data %>%
      mutate(
        sample      = SAMPLE_ID,
        copy_number = absolute_cn
      ) %>%
      select(sample, chromosome, start, end,
             copy_number, log2_ratio, DATAPOINTS) %>%
      arrange(as.integer(chromosome), start)

    seg_out <- file.path(sample_out,
                         paste0(SAMPLE_ID,
                                "_constrained_abs_cn.seg"))
    write.table(output_segs, seg_out,
                sep = "\t", quote = FALSE, row.names = FALSE)
    cat("    Seg file saved:", seg_out, "\n")

    cat("  Step 6: Saving genome plot...\n")
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
                           paste0(SAMPLE_ID,
                                  "_constrained_genome_plot.pdf"))
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
        title = paste0(SAMPLE_ID, " CONSTRAINED",
                       " | ploidy=",      round(best_ploidy, 2),
                       " | cellularity=", round(best_cellularity, 2),
                       " | MAD=",         round(best_MAD, 4),
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
    cat("    Plot saved:", plot_file, "\n")

    list(
      sample                  = SAMPLE_ID,
      ichorcna_tf             = ICHORCNA_TF,
      ichorcna_ploidy         = ICHORCNA_PLY,
      unc_ploidy              = sample_table$rascal_unc_ploidy[i],
      unc_cellularity         = sample_table$rascal_unc_cell[i],
      unc_MAD                 = sample_table$rascal_unc_MAD[i],
      unc_expected_tf         = round(tf_from_unconstrained, 3),
      unc_tf_diff             = round(diff_unconstrained, 3),
      con_ploidy              = round(best_ploidy, 3),
      con_cellularity         = round(best_cellularity, 3),
      con_MAD                 = round(best_MAD, 4),
      con_expected_tf         = round(tf_from_constrained, 3),
      con_tf_diff             = round(diff_constrained, 3),
      improvement             = round(improvement, 3),
      con_n_solutions         = n_solutions,
      cell_range_used         = paste0("[", round(cell_min, 2),
                                       ",", round(cell_max, 2), "]"),
      status                  = "SUCCESS"
    )

  }, error = function(e) {
    cat("  ERROR:", conditionMessage(e), "\n")
    list(
      sample          = SAMPLE_ID,
      ichorcna_tf     = ICHORCNA_TF,
      ichorcna_ploidy = ICHORCNA_PLY,
      unc_ploidy      = sample_table$rascal_unc_ploidy[i],
      unc_cellularity = sample_table$rascal_unc_cell[i],
      unc_MAD         = sample_table$rascal_unc_MAD[i],
      unc_expected_tf = NA,
      unc_tf_diff     = NA,
      con_ploidy      = NA,
      con_cellularity = NA,
      con_MAD         = NA,
      con_expected_tf = NA,
      con_tf_diff     = NA,
      improvement     = NA,
      con_n_solutions = NA,
      cell_range_used = paste0("[", round(cell_min, 2),
                               ",", round(cell_max, 2), "]"),
      status          = paste("FAILED:", conditionMessage(e))
    )
  })

  summary_rows[[SAMPLE_ID]] <- result
  cat("  Completed:", SAMPLE_ID, "\n\n")
}

# MASTER COMPARISON TABLE
cat("═══════════════════════════════════════════════════════\n")
cat("Building comparison table...\n")

summary_df <- bind_rows(summary_rows)

# The IMPROVED branch is only reached after both agreement thresholds have failed, so a
# sample landing there still disagrees with ichorCNA by more than 0.10 - it has merely moved
# closer than the unconstrained fit was. The label alone read as a pass, which mattered for
# B3_S09: it improved by 0.184 but ends 0.199 away from ichorCNA, outside the acceptance
# band. The wording now states that explicitly.
summary_df <- summary_df %>%
  mutate(interpretation = case_when(
    is.na(con_MAD)         ~ "FAILED",
    con_tf_diff <= 0.05    ~ "GOOD — constrained agrees with ichorCNA",
    con_tf_diff <= 0.10    ~ "ACCEPTABLE — constrained close to ichorCNA",
    improvement > 0.10     ~ "IMPROVED but still discordant — closer than unconstrained, still >0.10 from ichorCNA",
    TRUE                   ~ "POOR — neither solution agrees with ichorCNA"
  ))

summary_file <- file.path(OUTPUT_DIR,
                           "RASCAL_batch3_constrained_summary.csv")
write.csv(summary_df, summary_file, row.names = FALSE)
cat("Summary saved:", summary_file, "\n\n")

cat("═══════════════════════════════════════════════════════\n")
cat("SUPERVISOR COMPARISON TABLE — BATCH 3 CONSTRAINED\n")
cat("═══════════════════════════════════════════════════════\n\n")

cat(sprintf("%-8s | %6s | %10s | %10s | %8s | %8s | %s\n",
    "Sample", "iCNA_TF",
    "UNC_exp_TF", "CON_exp_TF",
    "UNC_diff", "CON_diff",
    "Interpretation"))
cat(strrep("-", 90), "\n")

for (i in 1:nrow(summary_df)) {
  r <- summary_df[i, ]
  cat(sprintf("%-8s | %6.3f | %10.3f | %10.3f | %8.3f | %8.3f | %s\n",
      r$sample, r$ichorcna_tf,
      ifelse(is.na(r$unc_expected_tf), NA, r$unc_expected_tf),
      ifelse(is.na(r$con_expected_tf), NA, r$con_expected_tf),
      ifelse(is.na(r$unc_tf_diff), NA, r$unc_tf_diff),
      ifelse(is.na(r$con_tf_diff), NA, r$con_tf_diff),
      r$interpretation))
}

cat("\n")
cat("Successful:", sum(summary_df$status == "SUCCESS"), "/ 3\n")
cat("GOOD/ACCEPTABLE constrained solutions:",
    sum(summary_df$interpretation %in%
          c("GOOD — constrained agrees with ichorCNA",
            "ACCEPTABLE — constrained close to ichorCNA"),
        na.rm = TRUE), "\n")
cat("═══════════════════════════════════════════════════════\n")
