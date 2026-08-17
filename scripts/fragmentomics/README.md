# Fragmentomics

Fragment-length metric extraction and tumour fraction correlation analysis for all 34 plasma samples.

Fragment lengths reflect nucleosomal architecture of the cell of origin; tumour-derived cfDNA fragments are systematically shorter than background cfDNA. Four metrics were calculated per sample: median fragment length, peak fragment length, short:long ratio (100–150 bp / 151–220 bp), and long-fragment fraction (>250 bp). Batch 3 BAMs were supplied without duplicate marking and required post-hoc Picard MarkDuplicates before fragmentomics analysis (thesis Section 2.16.3). B2_S11 was identified as an outlier across three independent lines of evidence and sensitivity-tested with and without inclusion.

## Execution order

1. `run_batch3_markdup.slurm` — Picard MarkDuplicates for Batch 3 BAMs only (Batches 1+2 pre-supplied)
2. `run_fragment_extraction.slurm` — fragment-length histogram extraction, all 34 samples
3. `calculate_fragment_metrics.R` / `run_fragment_metrics.slurm` — per-sample metrics and Batch 3 deduplication sensitivity check
4. `generate_publication_figures.R` / `run_publication_figures.slurm` — publication-quality figures
5. `tf_correlation.R` / `run_tf_correlation.slurm` — Spearman TF correlation, BH-corrected (n=34 primary; n=33 sensitivity)

## Inputs
BAM files (not included in this repository as it is patient-level data;); `fragment_sample_sheet.tsv` — sample to BAM path lookup

## Notes
SLURM scripts are specific to the ALICE HPC environment and will require path updates for other systems. Batch 3 duplicate marking was applied for fragmentomics only; the CNA arm uses original Batch 3 BAMs as supplied, as CNA pipline does not affected by this.

## Outputs
`results/fragmentomics/` : per-sample metric tables, deduplication sensitivity figures, TF correlation plots and tables
