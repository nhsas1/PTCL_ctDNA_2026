# PTCL ctDNA — Shallow WGS Analysis

**Thesis title:** Shallow whole genome sequencing of ctDNA in
T-cell lymphoma focusing on copy number variations and fragmentomics

**Author:** Noor Shaban
**Institution:** University of Leicester
**Year:** 2026

## Project Overview

This repository contains the analysis scripts and QC documentation
for shallow whole genome sequencing (sWGS) of circulating tumour DNA
(ctDNA) from peripheral T-cell lymphoma (PTCL) patients.

Two independent analyses run on the same libraries. **Copy number alteration (CNA)
profiling** asks which parts of the genome are gained or lost, and is restricted to the
subset of samples with detectable tumour. **Fragmentomics** asks whether the physical
length distribution of cfDNA fragments shifts with tumour burden, and uses all samples,
since it does not depend on calling a copy number state.

## Cohort

34 plasma cfDNA samples across three sequencing batches. Paired-end 150bp, aligned with
BWA-MEM to hg38. Coverage 0.87–1.77x.

| Batch | n | Deduplicated at source |
|---|---|---|
| 1 | 5  | yes |
| 2 | 16 | yes |
| 3 | 13 | **no** |

**The Batch 3 deduplication asymmetry explains differences throughout the pipeline** and is
worth reading once here rather than rediscovering per script. Because Batch 3 BAMs carry no
duplicate flags, SAMtools filtering differs: `-F 3328` for Batches 1–2 (secondary +
supplementary + duplicate) against `-F 2304` for the original Batch 3 BAMs (secondary +
supplementary only).

The fragmentomics arm additionally produces a deduplicated Batch 3 set with Picard
(`run_batch3_markdup.slurm`) and runs a sensitivity analysis against it. The CNA arm reads
the original non-deduplicated BAMs. **The two arms therefore treat Batch 3 differently**;
this was a choice rather than a constraint, and should be described as such.

## Cohort derivation for the CNA arm

The n=20 figure used throughout the CNA analysis comes from two filters applied in
`scripts/GISTIC/prepare_GISTIC_input.R`:

```
 34  total samples
-13  tumour fraction below the 3% ichorCNA detection floor
- 1  B3_S12, excluded for ploidy non-identifiability
=20  analysable for CNA
```

Below the 3% floor, copy number calls are not separable from noise; in a recurrence
analysis, arbitrary calls dilute real peaks and can manufacture spurious ones. B3_S12
clears the floor at TF 0.0903 but is excluded separately, because its ploidy cannot be
identified from coverage alone and so its absolute copy numbers are untrustworthy.

Note the phrasing "11 below floor + 2 borderline + 1 ploidy" used elsewhere describes the
**same 14 exclusions**; the two borderline samples (B1_S04 at 0.0292 and B2_S14 at 0.0251)
sit just under the cutoff and are counted separately in that narrative. Both framings give
20 analysable samples.

The fragmentomics arm uses all 34 samples and does not apply these exclusions.

## Pipeline execution order

Dependency order is not obvious from the filenames. Two ordering hazards are marked below.

### CNA arm

1. `QDNAseq/run_QDNAseq.slurm` and `run_QDNAseq_batch3.slurm` — bin read counts, GC
   correction, segmentation
2. `ichorCNA/run_readCounter_all.slurm` and `run_readCounter_batch3.slurm` — 1Mb WIG files
3. `ichorCNA/run_ichorCNA_all.slurm` and `run_ichorCNA_batch3.slurm` — joint tumour
   fraction and ploidy estimation; this is the primary CNA caller
4. `ichorCNA/run_anomaly_scan.slurm` — screens for segments physically inconsistent with
   the fitted solution, producing the shortlist for manual curation
5. Manual curation of B2_S08, B2_S15, B3_S09 → `results/ichorCNA/corrected_seg_files/`
6. `ichorCNA/calculate_FGA_AS.R` — fraction genome altered and Taylor aneuploidy score
7. `GISTIC/prepare_GISTIC_input.R` → GISTIC2 (Apptainer) → `GISTIC/per_sample_GISTIC.R`
8. `ichorCNA/CNA_gene_mapping.R` and `plot_cnFreq_cohort_v2.R` — gene mapping and cohort
   figures

Validation runs alongside: `RASCAL/` cross-checks ploidy and cellularity,
`ichorCNA/plot_bin_size_comparison.R` checks 1Mb against 500kb bins, and
`ichorCNA/plot_qdnaseq_ichorcna_concordance_v2.R` compares the two callers.

### Fragmentomics arm

1. `run_fragment_extraction.slurm` — fragment length histograms from TLEN
2. `run_fragment_metrics.slurm` — median, peak, short:long ratio, long-fragment fraction
3. `run_group_comparison.slurm` — Kruskal-Wallis plus pairwise Wilcoxon, BH-corrected
4. `run_tf_correlation.slurm` — Spearman against ichorCNA tumour fraction
5. `compute_gc_content_bins.sh` → `run_genomewide_binning.slurm` →
   `run_genomewide_profile.slurm` — DELFI-style 5Mb profiling with LOESS GC correction
6. `run_final_profile_plot.slurm` (5d) — centromere exclusion and z-score recomputation
7. `run_perbin_significance.slurm` (5e) — per-bin testing, BH across all retained bins
8. `run_telomere_check.slurm` (5f)

**Ordering hazard 1:** step 5c (`run_fix_profile_plot.slurm`) corrects chromosome ordering
in the step 5 figure but writes the same filenames, so running step 5 afterwards silently
restores the incorrect version.

**Ordering hazard 2:** step 5e consumes step 5d's z-scores and must be re-run after it.

## Environment

Modules on ALICE: `samtools/1.18`, `picard/3.0.0`, `R/4.3.1`, GISTIC2 via Apptainer,
`gcc/12.3.0`.

R packages are bootstrapped by `scripts/QDNAseq/install_packages.R`, which installs the
Bioconductor and CRAN dependencies into a user library. **`rascal` is not on CRAN or
Bioconductor and must be installed separately** before anything in `scripts/RASCAL/` runs.

All BAMs and raw data live only on ALICE (`/scratch/alice/n/nhsas1/PTCL/BAMs/` and
`BAMs_batch3/`) and are excluded from git. Nothing in this repository can be regenerated
locally.

## Stated limitations

- **No mappability correction.** QDNAseq runs with `residual=FALSE`, `blacklist=FALSE` and
  no mappability input, so low-mappability bins carry artificially low counts.
- **Batch 3 duplicates retained in the CNA arm**, though a deduplicated set exists and the
  fragmentomics arm uses it.
- **Two mapping-quality thresholds.** QDNAseq uses its default MAPQ 37; fragmentomics
  filters at `-q 20`.
- **Ploidy estimated from coverage alone.** Without allele fractions, near-diploid and
  near-triploid/WGD solutions can be indistinguishable — the basis for excluding B3_S12.
- **chrX and chrY excluded throughout.** ichorCNA ran with `--chrs "c(1:22)"`, and QDNAseq
  excludes them implicitly via an unset `applyFilters(chromosomes=)` default.
- **Copy-neutral loss of heterozygosity is not detectable** by sWGS and is not counted as
  altered in FGA.

### Reproducibility gaps

- `hg38_bins_15kb_annotated.rds` is required by six scripts and **has no generating script
  in this repository**.
- `CNA_gene_mapping.R` queries Ensembl without pinning a release, so gene coordinates can
  differ between runs.
- Outputs of the B2_S08, B2_S15 and B3_S12 diagnostic ichorCNA runs are not committed, so
  three curation decisions and one exclusion cannot be verified from this repository.

## Repository Structure

```
scripts/
  QDNAseq/        relative copy number, bin-level and segmented
  ichorCNA/       primary CNA caller; FGA, aneuploidy score, gene mapping, figures
  RASCAL/         independent ploidy and cellularity cross-check
  GISTIC/         cohort-level recurrence analysis
  fragmentomics/  fragment size analysis, five sequential steps
results/
  ichorCNA/       segments, params, curated files, FGA/AS, GenVisR, concordance
  GISTIC/         input seg, lesion tables, peak frequencies, figures
  RASCAL/         ploidy and cellularity summaries
  fragmentomics/  metrics, statistical tables, figures
QC_results/       BAM QC, coverage, duplicate rates, sample metadata
```

## Further reading in this repository

- [`RERUN_REQUIRED.md`](RERUN_REQUIRED.md) — outputs that predate a bug fix and need
  regenerating, separated from those verified unaffected, plus open data-integrity
  questions.
- [`CLEANUP_PROPOSAL.md`](CLEANUP_PROPOSAL.md) — proposed structural changes, with the
  applied items marked.
- [`CLAUDE.md`](CLAUDE.md) — project conventions and key technical facts.
