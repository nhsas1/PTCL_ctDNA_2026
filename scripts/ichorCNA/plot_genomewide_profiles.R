#!/usr/bin/env Rscript
# Genome-wide copy number plots for manually curated samples, matching
# ichorCNA's native bin-level scatter style (log2 ratio per 1Mb bin),
# but coloured according to the CORRECTED segment calls rather than
# the original erroneous ichorCNA calls

library(ggplot2)
library(dplyr)
library(readr)

depth_paths <- list(
  B2_S08 = "/scratch/alice/n/nhsas1/PTCL/ichorCNA/output/B2_S08/B2_S08.correctedDepth.txt",
  B2_S15 = "/scratch/alice/n/nhsas1/PTCL/ichorCNA/output/B2_S15/B2_S15.correctedDepth.txt",
  B3_S09 = "/scratch/alice/n/nhsas1/PTCL/ichorCNA/output_batch3/B3_S09/B3_S09.correctedDepth.txt"
)

corrected_dir <- "/scratch/alice/n/nhsas1/PTCL/ichorCNA/corrected_seg_files"
output_dir <- "/scratch/alice/n/nhsas1/PTCL/ichorCNA/corrected_seg_files/genomewide_plots"
dir.create(output_dir, showWarnings = FALSE)

# hg38 chromosome lengths for cumulative x-axis positioning
chr_lengths <- c(
  248956422, 242193529, 198295559, 190214555, 181538259,
  170805979, 159345973, 145138636, 138394717, 133797422,
  135086622, 133275309, 114364328, 107043718, 101991189,
  90338345, 83257441, 80373285, 58617616, 64444167,
  46709983, 50818468
)
names(chr_lengths) <- 1:22
chr_offset <- c(0, cumsum(chr_lengths)[1:21])
names(chr_offset) <- 1:22

# Assign each 1Mb bin to the corrected segment that contains it
assign_bin_calls <- function(bins, seg) {
  bins$corrected_call <- NA_character_
  for (i in 1:nrow(seg)) {
    match_idx <- which(
      bins$chr == as.character(seg$chrom[i]) &
      bins$start >= seg$start[i] &
      bins$end <= seg$end[i]
    )
    bins$corrected_call[match_idx] <- seg$Corrected_Call[i]
  }
  bins
}

plot_sample <- function(sample_id) {
  # Bin-level log2 ratios (ichorCNA's own corrected depth output)
  bins <- read_tsv(depth_paths[[sample_id]], show_col_types = FALSE) %>%
    filter(chr %in% as.character(1:22), !is.na(log2_TNratio_corrected)) %>%
    mutate(chr = as.character(chr))

  # Corrected segment calls
  seg <- read.table(file.path(corrected_dir, paste0(sample_id, ".seg.corrected.txt")),
                     header = TRUE, sep = "\t", stringsAsFactors = FALSE) %>%
    filter(chrom %in% 1:22) %>%
    mutate(chrom = as.character(chrom))

  bins <- assign_bin_calls(bins, seg)
  bins <- bins %>% filter(!is.na(corrected_call))

  bins <- bins %>%
    mutate(
      x_pos = chr_offset[chr] + (start + end) / 2,
      point_colour = case_when(
        corrected_call == "HETD" ~ "#2166AC",
        corrected_call == "NEUT" ~ "grey55",
        corrected_call == "GAIN" ~ "#F4A582",
        corrected_call %in% c("AMP", "HLAMP") ~ "#B2182B",
        TRUE ~ "grey70"
      )
    )

  chr_mid <- chr_offset + chr_lengths / 2
  chr_boundaries <- cumsum(chr_lengths)

  y_range <- range(bins$log2_TNratio_corrected, na.rm = TRUE)
  y_limit <- max(abs(y_range)) * 1.1

  p <- ggplot(bins, aes(x = x_pos, y = log2_TNratio_corrected)) +
    geom_vline(xintercept = chr_boundaries, colour = "grey85", linewidth = 0.3) +
    geom_hline(yintercept = 0, colour = "grey70", linewidth = 0.3, linetype = "dashed") +
    geom_point(aes(colour = point_colour), size = 0.6, alpha = 0.75) +
    scale_colour_identity() +
    scale_x_continuous(
      breaks = chr_mid, labels = 1:22,
      expand = expansion(mult = 0.01)
    ) +
    scale_y_continuous(limits = c(-y_limit, y_limit)) +
    labs(
      title = paste0(sample_id, " \u2014 Corrected Genome-Wide Copy Number Profile"),
      subtitle = "1Mb bin-level log2 ratio, recoloured by manually curated segment calls",
      x = "Chromosome", y = "log2 Ratio (Tumour/Normal, Corrected)",
      caption = "Blue = loss (HETD) | Grey = neutral | Orange = gain | Red = amplification"
    ) +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(size = 9, colour = "grey40"),
      plot.caption = element_text(size = 8, colour = "grey50"),
      panel.grid.major.x = element_blank()
    )

  ggsave(file.path(output_dir, paste0(sample_id, "_corrected_genomewide_v2.png")),
         plot = p, width = 14, height = 4, dpi = 300)
  ggsave(file.path(output_dir, paste0(sample_id, "_corrected_genomewide_v2.pdf")),
         plot = p, width = 14, height = 4)

  cat(sprintf("%s: %d bins plotted, plot saved\n", sample_id, nrow(bins)))
}

for (s in names(depth_paths)) plot_sample(s)

cat(sprintf("\nAll plots saved to: %s\n", output_dir))
