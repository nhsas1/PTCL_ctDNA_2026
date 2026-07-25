library(dplyr)
library(readr)
library(ggplot2)
library(ggrepel)
library(patchwork)

scores_file   <- "/scratch/alice/n/nhsas1/PTCL/GISTIC/output/scores.gistic"
lesions_file  <- "/scratch/alice/n/nhsas1/PTCL/GISTIC/output/all_lesions.conf_99.txt"
seg_input_file <- "/scratch/alice/n/nhsas1/PTCL/GISTIC/input/PTCL_GISTIC_input.seg"
results_dir   <- "/scratch/alice/n/nhsas1/PTCL/GISTIC/output/"

n_samples <- read_tsv(seg_input_file, show_col_types = FALSE) %>%
  summarise(n = n_distinct(Sample)) %>%
  pull(n)
cat(sprintf("Detected n = %d samples from GISTIC input file\n", n_samples))

scores <- read_tsv(scores_file, show_col_types = FALSE) %>%
  rename(
    type      = Type,
    chr       = Chromosome,
    start     = Start,
    end       = End,
    log10q    = `-log10(q-value)`,
    gscore    = `G-score`,
    amplitude = `average amplitude`,
    frequency = frequency
  ) %>%
  mutate(
    chr    = as.integer(trimws(chr)),
    log10q = as.numeric(log10q),
    gscore = as.numeric(gscore)
  ) %>%
  filter(!is.na(chr), chr %in% 1:22)

chr_sizes <- scores %>%
  group_by(chr) %>%
  summarise(chr_len = max(end), .groups = "drop") %>%
  arrange(chr) %>%
  mutate(cum_start = lag(cumsum(chr_len), default = 0))

scores <- scores %>%
  left_join(chr_sizes, by = "chr") %>%
  mutate(genome_pos = (start + end) / 2 + cum_start)

chr_labels <- chr_sizes %>%
  mutate(mid = cum_start + chr_len / 2)

chr_rects <- chr_sizes %>%
  mutate(
    xmin = cum_start,
    xmax = cum_start + chr_len,
    fill = ifelse(chr %% 2 == 0, "grey97", "white")
  )

amp_data <- scores %>% filter(type == "Amp")
del_data <- scores %>% filter(type == "Del")

# Peak labels read dynamically from all_lesions.conf_99.txt
# (previously hardcoded to the old n=21 peak set)
lesions_raw <- read_tsv(lesions_file, show_col_types = FALSE)

peaks <- lesions_raw %>%
  filter(!grepl("CN values", `Unique Name`)) %>%
  mutate(
    peak_type = ifelse(grepl("Amplification", `Unique Name`), "Amp", "Del"),
    cytoband  = trimws(Descriptor),
    qval      = as.numeric(`q values`),
    chr       = as.integer(sub("chr(\\d+):.*", "\\1", `Wide Peak Limits`)),
    coords    = sub("^chr\\d+:", "", `Wide Peak Limits`),
    coords    = trimws(sub("\\(.*", "", coords)),
    peak_start = as.numeric(sub("-.*", "", coords)),
    peak_end   = as.numeric(sub(".*-", "", coords)),
    mid        = (peak_start + peak_end) / 2
  ) %>%
  filter(!is.na(chr), chr %in% 1:22) %>%
  left_join(chr_sizes, by = "chr") %>%
  mutate(genome_pos = mid + cum_start,
         log10q     = -log10(qval)) %>%
  select(peak_type, chr, cytoband, genome_pos, log10q, qval)

cat("\nPeaks that will be labelled on the plots:\n")
print(as.data.frame(peaks))

amp_peaks_label <- peaks %>% filter(peak_type == "Amp")
del_peaks_label <- peaks %>% filter(peak_type == "Del")

if (nrow(amp_peaks_label) == 0) {
  cat("\nNo significant amplification peaks - amplification plot will show the curve with no labelled peaks.\n")
}

palette_amp <- "#C0392B"
palette_del <- "#2471A3"
sig_line_col <- "grey40"

qthresh_25 <- -log10(0.25)
qthresh_10 <- -log10(0.10)
qthresh_01 <- -log10(0.01)

y_max <- 8.5

p_amp <- ggplot() +
  geom_rect(data = chr_rects,
            aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = y_max, fill = fill),
            inherit.aes = FALSE) +
  scale_fill_identity() +
  geom_hline(yintercept = qthresh_25, linetype = "dashed",
             colour = sig_line_col, linewidth = 0.4) +
  geom_hline(yintercept = qthresh_01, linetype = "dotted",
             colour = "black", linewidth = 0.5) +
  geom_line(data = amp_data,
            aes(x = genome_pos, y = pmin(log10q, y_max - 0.2)),
            colour = palette_amp, linewidth = 0.55, alpha = 0.85) +
  geom_area(data = amp_data,
            aes(x = genome_pos, y = pmin(log10q, y_max - 0.2)),
            fill = palette_amp, alpha = 0.18) +
  geom_point(data = amp_peaks_label,
             aes(x = genome_pos, y = pmin(log10q, y_max - 0.3)),
             colour = palette_amp, size = 2, shape = 18) +
  geom_text_repel(
    data          = amp_peaks_label,
    aes(x = genome_pos, y = pmin(log10q, y_max - 0.3), label = cytoband),
    size          = 2.8,
    colour        = palette_amp,
    fontface      = "bold",
    nudge_y       = 1.5,
    nudge_x       = 8e7,
    segment.colour = palette_amp,
    segment.size  = 0.3,
    segment.alpha = 0.6,
    max.overlaps  = 30
  ) +
  annotate("text",
           x     = max(chr_sizes$cum_start + chr_sizes$chr_len) * 0.01,
           y     = qthresh_25 + 0.18,
           label = "q = 0.25",
           size  = 2.5, colour = sig_line_col, hjust = 0) +
  annotate("text",
           x     = max(chr_sizes$cum_start + chr_sizes$chr_len) * 0.01,
           y     = qthresh_01 + 0.18,
           label = "q = 0.01",
           size  = 2.5, colour = "grey30", hjust = 0) +
  scale_x_continuous(breaks = chr_labels$mid,
                     labels = chr_labels$chr,
                     expand = c(0.005, 0)) +
  scale_y_continuous(limits = c(0, y_max),
                     breaks = c(0, 2, 4, 6, 8),
                     labels = c("0", "2", "4", "6", "8"),
                     expand = c(0, 0)) +
  labs(x = NULL,
       y = expression(-log[10](q)),
       title = "Amplifications") +
  theme_classic(base_size = 11) +
  theme(
    axis.text.x   = element_text(size = 8, colour = "grey30"),
    axis.text.y   = element_text(size = 9),
    axis.title.y  = element_text(size = 10),
    axis.line.x   = element_blank(),
    axis.ticks.x  = element_blank(),
    plot.title    = element_text(face = "bold", size = 11,
                                 colour = palette_amp, hjust = 0),
    panel.grid    = element_blank(),
    plot.margin   = margin(t = 6, r = 10, b = 0, l = 10)
  )

p_del <- ggplot() +
  geom_rect(data = chr_rects,
            aes(xmin = xmin, xmax = xmax, ymin = -y_max, ymax = 0, fill = fill),
            inherit.aes = FALSE) +
  scale_fill_identity() +
  geom_hline(yintercept = -qthresh_25, linetype = "dashed",
             colour = sig_line_col, linewidth = 0.4) +
  geom_hline(yintercept = -qthresh_01, linetype = "dotted",
             colour = "black", linewidth = 0.5) +
  geom_line(data = del_data,
            aes(x = genome_pos, y = -pmin(log10q, y_max - 0.2)),
            colour = palette_del, linewidth = 0.55, alpha = 0.85) +
  geom_area(data = del_data,
            aes(x = genome_pos, y = -pmin(log10q, y_max - 0.2)),
            fill = palette_del, alpha = 0.18) +
  geom_point(data = del_peaks_label,
             aes(x = genome_pos, y = -pmin(log10q, y_max - 0.3)),
             colour = palette_del, size = 2, shape = 18) +
  geom_text_repel(
    data          = del_peaks_label,
    aes(x = genome_pos, y = -pmin(log10q, y_max - 0.3), label = cytoband),
    size          = 2.8,
    colour        = palette_del,
    fontface      = "bold",
    nudge_y       = -1.2,
    nudge_x       = 8e7,
    segment.colour = palette_del,
    segment.size  = 0.3,
    segment.alpha = 0.6,
    max.overlaps  = 30
  ) +
  annotate("text",
           x     = max(chr_sizes$cum_start + chr_sizes$chr_len) * 0.01,
           y     = -qthresh_25 - 0.2,
           label = "q = 0.25",
           size  = 2.5, colour = sig_line_col, hjust = 0) +
  annotate("text",
           x     = max(chr_sizes$cum_start + chr_sizes$chr_len) * 0.01,
           y     = -qthresh_01 - 0.2,
           label = "q = 0.01",
           size  = 2.5, colour = "grey30", hjust = 0) +
  scale_x_continuous(breaks = chr_labels$mid,
                     labels = chr_labels$chr,
                     expand = c(0.005, 0)) +
  scale_y_continuous(limits = c(-y_max, 0),
                     breaks = c(-8, -6, -4, -2, 0),
                     labels = c("8", "6", "4", "2", "0"),
                     expand = c(0, 0)) +
  labs(x       = "Chromosome",
       y       = expression(-log[10](q)),
       title   = "Deletions",
       caption = sprintf("GISTIC2 | hg38 | q < 0.25 threshold | n = %d detectable PTCL ctDNA samples", n_samples)) +
  theme_classic(base_size = 11) +
  theme(
    axis.text.x   = element_text(size = 8, colour = "grey30"),
    axis.text.y   = element_text(size = 9),
    axis.title.y  = element_text(size = 10),
    axis.title.x  = element_text(size = 10, colour = "grey30"),
    plot.title    = element_text(face = "bold", size = 11,
                                 colour = palette_del, hjust = 0),
    panel.grid    = element_blank(),
    plot.caption  = element_text(size = 7.5, colour = "grey50", hjust = 1),
    plot.margin   = margin(t = 0, r = 10, b = 6, l = 10)
  )

p_combined <- p_amp / p_del +
  plot_annotation(
    title    = "GISTIC2 Copy Number Significance Landscape",
    subtitle = sprintf("PTCL ctDNA cohort | sWGS plasma cfDNA | detectable samples (TF >= 3%%, n = %d)", n_samples),
    theme    = theme(
      plot.title    = element_text(face = "bold", size = 14, hjust = 0),
      plot.subtitle = element_text(size = 10, colour = "grey40", hjust = 0),
      plot.margin   = margin(12, 12, 8, 12)
    )
  )

ggsave(filename = file.path(results_dir, "GISTIC_landscape.pdf"),
       plot = p_combined, width = 14, height = 8, device = "pdf")

ggsave(filename = file.path(results_dir, "GISTIC_landscape.png"),
       plot = p_combined, width = 14, height = 8, dpi = 300, device = "png")

p_freq <- ggplot() +
  geom_rect(data = chr_rects,
            aes(xmin = xmin, xmax = xmax,
                ymin = -1, ymax = 1, fill = fill),
            inherit.aes = FALSE) +
  scale_fill_identity() +
  geom_hline(yintercept = 0, colour = "grey50", linewidth = 0.4) +
  geom_line(data = amp_data,
            aes(x = genome_pos, y = frequency),
            colour = palette_amp, linewidth = 0.5, alpha = 0.8) +
  geom_area(data = amp_data,
            aes(x = genome_pos, y = frequency),
            fill = palette_amp, alpha = 0.2) +
  geom_line(data = del_data,
            aes(x = genome_pos, y = -frequency),
            colour = palette_del, linewidth = 0.5, alpha = 0.8) +
  geom_area(data = del_data,
            aes(x = genome_pos, y = -frequency),
            fill = palette_del, alpha = 0.2) +
  scale_x_continuous(breaks = chr_labels$mid,
                     labels = chr_labels$chr,
                     expand = c(0.005, 0)) +
  scale_y_continuous(limits = c(-1, 1),
                     breaks = c(-1, -0.5, 0, 0.5, 1),
                     labels = c("100%", "50%", "0", "50%", "100%"),
                     expand = c(0, 0)) +
  labs(x       = "Chromosome",
       y       = "Sample frequency",
       title   = "Copy Number Alteration Frequency",
       subtitle = sprintf("Red = gain  |  Blue = loss  |  PTCL ctDNA cohort (n = %d)", n_samples),
       caption = "GISTIC2 | hg38 | ichorCNA Corrected_Copy_Number") +
  theme_classic(base_size = 11) +
  theme(
    axis.text.x   = element_text(size = 8, colour = "grey30"),
    axis.text.y   = element_text(size = 9),
    axis.title    = element_text(size = 10),
    panel.grid    = element_blank(),
    plot.title    = element_text(face = "bold", size = 13, hjust = 0),
    plot.subtitle = element_text(size = 9.5, colour = "grey40", hjust = 0),
    plot.caption  = element_text(size = 7.5, colour = "grey50", hjust = 1),
    plot.margin   = margin(10, 12, 8, 12)
  )

ggsave(filename = file.path(results_dir, "GISTIC_frequency.pdf"),
       plot = p_freq, width = 14, height = 5, device = "pdf")

ggsave(filename = file.path(results_dir, "GISTIC_frequency.png"),
       plot = p_freq, width = 14, height = 5, dpi = 300, device = "png")

cat("Plots saved:\n")
cat(sprintf("  %sGISTIC_landscape.pdf\n", results_dir))
cat(sprintf("  %sGISTIC_frequency.pdf\n", results_dir))
cat("Done.\n")
