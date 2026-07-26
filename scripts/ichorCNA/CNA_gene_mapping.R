#!/usr/bin/env Rscript
# ============================================================================
# CNA-Gene Mapping — Patient-Driven Segment Annotation
# ============================================================================
# For each patient's altered segments, identifies cancer-relevant genes
# residing within the segment boundaries and annotates with published evidence.
#
# Input:  20 ichorCNA seg.txt files (corrected where applicable)
#         PTCL Gene List v2 (xlsx), cBioPortal CNA + mutation data
#         FGA/AS summary, arm-level scores
# Output: Per-segment gene annotation table + per-gene frequency summary
# ============================================================================

library(dplyr)
library(readr)
library(readxl)
library(biomaRt)
library(tidyr)
# --- FIX: reclaim select()/filter() from AnnotationDbi/biomaRt masking ---
select <- dplyr::select
filter <- dplyr::filter
rename <- dplyr::rename 

cat("=== CNA-Gene Mapping Pipeline ===\n\n")

# ============================================================================
# 1. FILE PATHS
# ============================================================================

batch12_dir   <- "/scratch/alice/n/nhsas1/PTCL/ichorCNA/output"
batch3_dir    <- "/scratch/alice/n/nhsas1/PTCL/ichorCNA/output_batch3"
corrected_dir <- "/scratch/alice/n/nhsas1/PTCL/ichorCNA/corrected_seg_files"
scores_dir    <- "/scratch/alice/n/nhsas1/PTCL/ichorCNA/aneuploidy_scores"
output_dir    <- "/scratch/alice/n/nhsas1/PTCL/ichorCNA/aneuploidy_scores"

gene_list_file <- "/scratch/alice/n/nhsas1/PTCL/ichorCNA/PTCL_Gene_List_v2.xlsx"
cna_genes_file <- "/scratch/alice/n/nhsas1/PTCL/ichorCNA/CNA_Genes.txt"
mut_genes_file <- "/scratch/alice/n/nhsas1/PTCL/ichorCNA/Mutated_Genes.txt"
fga_as_file    <- file.path(scores_dir, "FGA_AS_summary.csv")
arm_detail_file <- file.path(scores_dir, "arm_level_scores_detail.csv")

# ============================================================================
# 2. DEFINE SAMPLES AND PATHS
# ============================================================================

samples <- c(
  "B1_S01", "B1_S02", "B1_S05",
  "B2_S01", "B2_S02", "B2_S03", "B2_S04", "B2_S05", "B2_S06", "B2_S07",
  "B2_S08", "B2_S12", "B2_S15", "B2_S16",
  "B3_S01", "B3_S04", "B3_S05", "B3_S06", "B3_S08", "B3_S09"
)

corrected_samples <- c("B2_S08", "B2_S15", "B3_S09")

get_seg_path <- function(s) {
  if (s %in% corrected_samples) {
    file.path(corrected_dir, paste0(s, ".seg.corrected.txt"))
  } else if (grepl("^B3_", s)) {
    file.path(batch3_dir, s, paste0(s, ".seg.txt"))
  } else {
    file.path(batch12_dir, s, paste0(s, ".seg.txt"))
  }
}

# ============================================================================
# 3. READ ALL SEG FILES AND FILTER TO ALTERED SEGMENTS
# ============================================================================

cat("Reading seg files...\n")

all_segs <- do.call(rbind, lapply(samples, function(s) {
  f <- get_seg_path(s)
  seg <- read.table(f, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  data.frame(
    sample    = s,
    chr       = seg[[2]],
    start     = seg[[3]],
    end       = seg[[4]],
    num_mark  = seg[[5]],
    logR      = seg[[6]],
    orig_CN   = seg[[7]],
    orig_call = seg[[8]],
    corr_CN   = seg[[11]],
    corr_call = seg[[12]],
    stringsAsFactors = FALSE
  )
}))

altered <- all_segs %>%
  filter(corr_call != "NEUT") %>%
  mutate(
    segment_size_mb = round((end - start) / 1e6, 1),
    segment_label   = paste0("chr", chr, ":", round(start / 1e6), "-", round(end / 1e6), "Mb")
  )

cat(sprintf("Total segments: %d | Altered segments: %d\n\n",
            nrow(all_segs), nrow(altered)))

# ============================================================================
# 4. ADD PATIENT-LEVEL CONTEXT (TF, FGA, AS)
# ============================================================================

cat("Adding patient-level context...\n")

fga_as <- read_csv(fga_as_file, show_col_types = FALSE) %>%
  select(sample, tumor_fraction, FGA, AS_50)

altered <- altered %>%
  left_join(fga_as, by = "sample")

# ============================================================================
# 5. ADD ARM-LEVEL VALIDATION
# ============================================================================

cat("Adding arm-level validation...\n")

# hg38 centromere positions for arm assignment
centromeres <- data.frame(
  chr = 1:22,
  centromere = c(
    123400000, 93900000, 90900000, 50000000, 48800000, 58500000,
    60100000, 45200000, 43000000, 39800000, 53400000, 35500000,
    17000000, 17000000, 19000000, 36800000, 25100000, 18500000,
    26200000, 28100000, 12000000, 15000000
  )
)

# Determine which arm each segment's midpoint falls in
assign_arm <- function(chr_val, start_val, end_val) {
  midpoint <- (start_val + end_val) / 2
  acrocentric <- c(13, 14, 15, 21, 22)
  if (chr_val %in% acrocentric) return(paste0(chr_val, "q"))
  cent <- centromeres$centromere[centromeres$chr == chr_val]
  if (length(cent) == 0) return(NA_character_)
  if (midpoint < cent) paste0(chr_val, "p") else paste0(chr_val, "q")
}

altered$arm <- mapply(assign_arm, altered$chr, altered$start, altered$end)

arm_detail <- read_csv(arm_detail_file, show_col_types = FALSE) %>%
  select(sample, arm, status_50) %>%
  rename(arm_status_50 = status_50)

altered <- altered %>%
  left_join(arm_detail, by = c("sample", "arm")) %>%
  mutate(
    segment_arm_concordance = case_when(
      is.na(arm_status_50) ~ "N/A",
      corr_call %in% c("GAIN", "AMP", "HLAMP") & arm_status_50 == "GAIN" ~ "YES",
      corr_call == "HETD" & arm_status_50 == "LOSS" ~ "YES",
      arm_status_50 == "NEUTRAL" ~ "Focal (arm neutral)",
      TRUE ~ "Discordant"
    )
  )

# ============================================================================
# 6. BUILD COMBINED GENE REFERENCE
# ============================================================================

cat("Building gene reference from three sources...\n")

# Source 1: PTCL Gene List v2 (xlsx)
ptcl_raw <- read_excel(gene_list_file, sheet = 1, skip = 2)
ptcl_genes <- ptcl_raw %>%
  filter(!is.na(`Gene Symbol`)) %>%
  transmute(
    gene = `Gene Symbol`,
    source_ptcl_list = TRUE,
    ptcl_functional_category = `Functional Category`,
    ptcl_cna_type = `CNA Type\n& Cytoband`,
    ptcl_alteration_type = `Alteration Type`
  )

# Source 2: cBioPortal CNA genes
cna_raw <- read_tsv(cna_genes_file, show_col_types = FALSE)
cna_genes <- cna_raw %>%
  transmute(
    gene = Gene,
    source_cna_cbioportal = TRUE,
    cbioportal_cna_type = CNA,
    cbioportal_cna_cytoband = Cytoband,
    cbioportal_cna_freq = Freq
  )

# Source 3: cBioPortal mutated genes
mut_raw <- read_tsv(mut_genes_file, show_col_types = FALSE)
mut_genes <- mut_raw %>%
  transmute(
    gene = Gene,
    source_mut_cbioportal = TRUE,
    cbioportal_mut_freq = Freq
  )

# Combine unique gene symbols
all_gene_symbols <- unique(c(ptcl_genes$gene, cna_genes$gene, mut_genes$gene))
cat(sprintf("Unique genes across all sources: %d\n", length(all_gene_symbols)))

# ============================================================================
# 7. FETCH hg38 COORDINATES VIA biomaRt
# ============================================================================

cat("Fetching hg38 coordinates from Ensembl via biomaRt...\n")

# The default Ensembl mart is GRCh38, which is the correct build to match hg38 segment
# coordinates. GRCh37 would require host = "grch37.ensembl.org", so the build is right.
#
# REPRODUCIBILITY GAP: no Ensembl release is pinned, so this queries whatever version is
# live at run time. Gene coordinates and HGNC symbol assignments change between releases,
# which means re-running this script months apart can produce different gene lists from
# identical segment input. For a thesis result that should be reproducible, either pin the
# release (useEnsembl(..., version = NNN)) or record the release used alongside the output.
# The committed gene mapping does not state which release produced it.
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

coords <- getBM(
  attributes = c("hgnc_symbol", "chromosome_name", "start_position", "end_position"),
  filters    = "hgnc_symbol",
  values     = all_gene_symbols,
  mart       = ensembl
)

# Filter to standard chromosomes only (remove patches, scaffolds)
standard_chr <- c(as.character(1:22), "X", "Y")
coords <- coords %>%
  filter(chromosome_name %in% standard_chr) %>%
  rename(gene = hgnc_symbol, gene_chr = chromosome_name,
         gene_start = start_position, gene_end = end_position) %>%
  # Remove duplicates (keep longest transcript region per gene per chromosome)
  group_by(gene, gene_chr) %>%
  summarise(gene_start = min(gene_start), gene_end = max(gene_end), .groups = "drop") %>%
  # X is remapped to 23 so the column can be integer for joining against segment chr.
  #
  # Note Y is admitted by the filter above but has no such remapping, so as.integer("Y")
  # yields NA with a coercion warning and any Y gene is silently unmatchable. This has no
  # effect on results: ichorCNA ran with --chrs "c(1:22)", so no X or Y segments exist and
  # neither X nor Y genes could ever match regardless. Including them in standard_chr is
  # therefore redundant rather than harmful, but the Y path is the reason for any
  # "NAs introduced by coercion" warning in the log.
  mutate(gene_chr = as.integer(ifelse(gene_chr == "X", "23", gene_chr)))

cat(sprintf("Coordinates retrieved for %d genes on standard chromosomes\n",
            nrow(coords)))

# Identify genes that failed to map
missing <- setdiff(all_gene_symbols, coords$gene)
if (length(missing) > 0) {
  cat(sprintf("WARNING: %d genes not found in Ensembl: %s\n",
              length(missing), paste(head(missing, 10), collapse = ", ")))
}

# ============================================================================
# 8. FIND OVERLAPPING GENES FOR EACH ALTERED SEGMENT
# ============================================================================

cat("Finding gene overlaps for altered segments...\n")

# For each altered segment, find genes whose coordinates overlap
find_overlapping_genes <- function(seg_chr, seg_start, seg_end, gene_coords) {
  matches <- gene_coords %>%
    filter(gene_chr == seg_chr,
           gene_start < seg_end,
           gene_end > seg_start)
  if (nrow(matches) == 0) return(NULL)
  return(matches)
}

# Build annotation for each segment
annotation_list <- list()

for (i in 1:nrow(altered)) {
  seg <- altered[i, ]
  matches <- find_overlapping_genes(seg$chr, seg$start, seg$end, coords)

  if (is.null(matches) || nrow(matches) == 0) {
    annotation_list[[i]] <- data.frame(
      row_id = i,
      overlapping_genes = "none",
      n_genes = 0,
      gene_details = "-",
      gene_cytobands = "-",
      published_cna_freq = "-",
      published_mut_freq = "-",
      evidence_tier = "No known cancer gene",
      in_ptcl_list = FALSE,
      in_cbioportal_cna = FALSE,
      in_cbioportal_mut = FALSE,
      stringsAsFactors = FALSE
    )
  } else {
    # Annotate each matched gene
    gene_annotations <- lapply(matches$gene, function(g) {
      in_ptcl <- g %in% ptcl_genes$gene
      in_cna  <- g %in% cna_genes$gene
      in_mut  <- g %in% mut_genes$gene

      # Get published frequencies
      cna_freq <- if (in_cna) {
        cf <- cna_genes %>% filter(gene == g)
        paste(cf$cbioportal_cna_type, cf$cbioportal_cna_freq, collapse = "; ")
      } else "-"

      mut_freq <- if (in_mut) {
        mf <- mut_genes %>% filter(gene == g)
        paste(unique(mf$cbioportal_mut_freq), collapse = "; ")
      } else "-"

      # Get functional category
      func_cat <- if (in_ptcl) {
        pc <- ptcl_genes %>% filter(gene == g)
        pc$ptcl_functional_category[1]
      } else if (in_cna || in_mut) {
        "OncoKB cancer gene"
      } else "-"

      # Get cytoband from cBioPortal CNA data
      cytoband <- if (in_cna) {
        cb <- cna_genes %>% filter(gene == g)
        cb$cbioportal_cna_cytoband[1]
      } else "-"

      # Assign evidence tier
      tier <- if (in_ptcl) {
        "Tier 1 (PTCL gene list)"
      } else if (in_cna) {
        "Tier 2 (cBioPortal CNA)"
      } else if (in_mut) {
        "Tier 3 (cBioPortal mutated)"
      } else {
        "Tier 4 (Ensembl cancer gene)"
      }

      list(gene = g, func_cat = func_cat, cytoband = cytoband,
           cna_freq = cna_freq, mut_freq = mut_freq, tier = tier,
           in_ptcl = in_ptcl, in_cna = in_cna, in_mut = in_mut)
    })

    # Collapse into single row
    genes_str     <- paste(sapply(gene_annotations, `[[`, "gene"), collapse = ", ")
    details_str   <- paste(sapply(gene_annotations, function(x)
      paste0(x$gene, " (", x$func_cat, ")")), collapse = "; ")
    cytobands_str <- paste(unique(sapply(gene_annotations, `[[`, "cytoband")), collapse = "; ")
    cna_freq_str  <- paste(sapply(gene_annotations, function(x)
      if (x$cna_freq != "-") paste0(x$gene, ": ", x$cna_freq) else NULL), collapse = "; ")
    mut_freq_str  <- paste(sapply(gene_annotations, function(x)
      if (x$mut_freq != "-") paste0(x$gene, ": ", x$mut_freq) else NULL), collapse = "; ")

    # Clean empty strings
    if (cna_freq_str == "") cna_freq_str <- "-"
    if (mut_freq_str == "") mut_freq_str <- "-"

    # Determine highest evidence tier
    tiers <- sapply(gene_annotations, `[[`, "tier")
    best_tier <- if (any(grepl("Tier 1", tiers))) "Tier 1 (PTCL gene)"
                 else if (any(grepl("Tier 2", tiers))) "Tier 2 (cBioPortal CNA)"
                 else if (any(grepl("Tier 3", tiers))) "Tier 3 (cBioPortal mutated)"
                 else "Tier 4 (Ensembl only)"

    annotation_list[[i]] <- data.frame(
      row_id = i,
      overlapping_genes = genes_str,
      n_genes = nrow(matches),
      gene_details = details_str,
      gene_cytobands = cytobands_str,
      published_cna_freq = cna_freq_str,
      published_mut_freq = mut_freq_str,
      evidence_tier = best_tier,
      in_ptcl_list = any(sapply(gene_annotations, `[[`, "in_ptcl")),
      in_cbioportal_cna = any(sapply(gene_annotations, `[[`, "in_cna")),
      in_cbioportal_mut = any(sapply(gene_annotations, `[[`, "in_mut")),
      stringsAsFactors = FALSE
    )
  }
}

annotations <- do.call(rbind, annotation_list)

# ============================================================================
# 9. COMBINE INTO FINAL TABLE
# ============================================================================

cat("Building final output table...\n")

final <- altered %>%
  mutate(row_id = row_number()) %>%
  left_join(annotations, by = "row_id") %>%
  select(
    # Patient block
    sample, tumor_fraction, FGA, AS_50,
    # Segment block
    chr, start, end, segment_size_mb, segment_label, logR,
    corr_CN, corr_call, orig_CN, orig_call,
    # Arm validation
    arm, arm_status_50, segment_arm_concordance,
    # Gene annotation
    overlapping_genes, n_genes, gene_details, gene_cytobands,
    published_cna_freq, published_mut_freq,
    evidence_tier, in_ptcl_list, in_cbioportal_cna, in_cbioportal_mut
  ) %>%
  arrange(sample, chr, start)

# ============================================================================
# 10. GENERATE PER-GENE FREQUENCY SUMMARY
# ============================================================================

cat("Generating per-gene frequency summary...\n")

# Extract individual genes from comma-separated overlapping_genes column
gene_patient_pairs <- final %>%
  filter(overlapping_genes != "none") %>%
  select(sample, overlapping_genes, corr_call, tumor_fraction) %>%
  separate_rows(overlapping_genes, sep = ", ") %>%
  rename(gene = overlapping_genes)

# Calculate per-gene CNA frequency
gene_freq <- gene_patient_pairs %>%
  mutate(
    direction = case_when(
      corr_call %in% c("GAIN", "AMP", "HLAMP") ~ "Gain",
      corr_call == "HETD" ~ "Loss",
      TRUE ~ "Other"
    )
  ) %>%
  group_by(gene) %>%
  summarise(
    n_patients_altered = n_distinct(sample),
    freq_pct = round(n_distinct(sample) / 20 * 100, 1),
    n_gain = n_distinct(sample[direction == "Gain"]),
    n_loss = n_distinct(sample[direction == "Loss"]),
    patients = paste(sort(unique(sample)), collapse = ", "),
    .groups = "drop"
  ) %>%
  arrange(desc(n_patients_altered))

# Add source information
gene_freq <- gene_freq %>%
  mutate(
    in_ptcl_list = gene %in% ptcl_genes$gene,
    in_cbioportal_cna = gene %in% cna_genes$gene,
    in_cbioportal_mut = gene %in% mut_genes$gene,
    evidence_tier = case_when(
      in_ptcl_list ~ "Tier 1 (PTCL gene list)",
      in_cbioportal_cna ~ "Tier 2 (cBioPortal CNA)",
      in_cbioportal_mut ~ "Tier 3 (cBioPortal mutated)",
      TRUE ~ "Tier 4 (Ensembl only)"
    )
  )

# ============================================================================
# 11. SAVE OUTPUTS
# ============================================================================

write_csv(final,
          file.path(output_dir, "CNA_gene_mapping_all_segments.csv"))

write_csv(gene_freq,
          file.path(output_dir, "CNA_gene_frequency_summary.csv"))

cat(sprintf("\n=== OUTPUT FILES ===\n"))
cat(sprintf("  CNA_gene_mapping_all_segments.csv  — %d altered segments annotated\n",
            nrow(final)))
cat(sprintf("  CNA_gene_frequency_summary.csv     — %d genes with CNA in >= 1 patient\n",
            nrow(gene_freq)))
cat(sprintf("  Saved to: %s\n\n", output_dir))

# ============================================================================
# 12. SUMMARY STATISTICS
# ============================================================================

cat("=== SUMMARY ===\n\n")

cat(sprintf("Altered segments: %d across %d patients\n", nrow(final), length(samples)))
cat(sprintf("Segments with cancer gene overlap: %d (%.1f%%)\n",
            sum(final$n_genes > 0),
            sum(final$n_genes > 0) / nrow(final) * 100))
cat(sprintf("Segments with no known cancer gene: %d (%.1f%%)\n",
            sum(final$n_genes == 0),
            sum(final$n_genes == 0) / nrow(final) * 100))

cat(sprintf("\nUnique cancer genes in altered regions: %d\n", nrow(gene_freq)))
cat(sprintf("  Tier 1 (PTCL gene list):     %d\n", sum(gene_freq$in_ptcl_list)))
cat(sprintf("  Tier 2 (cBioPortal CNA):     %d\n",
            sum(gene_freq$in_cbioportal_cna & !gene_freq$in_ptcl_list)))
cat(sprintf("  Tier 3 (cBioPortal mutated):  %d\n",
            sum(gene_freq$in_cbioportal_mut & !gene_freq$in_cbioportal_cna & !gene_freq$in_ptcl_list)))

cat("\n=== TOP 20 GENES BY CNA FREQUENCY IN COHORT ===\n\n")
cat(sprintf("%-12s  Patients  Freq   Gain  Loss  Tier\n", "Gene"))
cat(paste(rep("-", 65), collapse = ""), "\n")

for (i in 1:min(20, nrow(gene_freq))) {
  r <- gene_freq[i, ]
  cat(sprintf("%-12s  %2d/20     %4.1f%%  %2d    %2d    %s\n",
              r$gene, r$n_patients_altered, r$freq_pct,
              r$n_gain, r$n_loss, r$evidence_tier))
}

# Arm concordance check
cat("\n=== ARM CONCORDANCE CHECK ===\n")
conc_table <- table(final$segment_arm_concordance)
for (status in names(conc_table)) {
  cat(sprintf("  %s: %d segments (%.1f%%)\n",
              status, conc_table[status],
              conc_table[status] / nrow(final) * 100))
}

cat("\n=== COMPLETE ===\n")
