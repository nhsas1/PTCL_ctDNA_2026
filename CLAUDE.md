# Project Context: PTCL ctDNA Thesis Repo

MSc Bioinformatics thesis project
Shallow whole-genome sequencing (sWGS) of plasma cfDNA from Peripheral T-cell Lymphoma (PTCL) patients.
Thesis title: "Shallow Whole Genome Sequencing of ctDNA in T-Cell Lymphoma."

## Cohort
34 plasma cfDNA samples, three batches:
- Batch 1 (n=5) and Batch 2 (n=16): paired-end 150bp, BWA-MEM/hg38, Picard MarkDuplicates already applied at source.
- Batch 3 (n=13): same alignment, NOT deduplicated at source.
- This batch asymmetry means SAMtools filtering flags differ: -F 3328 for Batches 1-2, -F 2304 for original Batch 3 BAMs.
- Coverage range 0.87-1.77x.

## Two analysis arms
1. **CNA (copy number alteration) profiling** — complete pipeline: ichorCNA -> manual curation -> FGA/AS (Taylor Aneuploidy Score) recalculation -> GISTIC2 -> CNA-gene mapping -> GenVisR visualisation. Final analysable cohort: 20/34 samples (14 excluded: 11 below 3% tumour fraction detection floor, 2 borderline with no anomalies, 1 unresolved ploidy ambiguity, B3_S12).
2. **Fragmentomics** — self-directed, no supervisor guidance. Five sequential steps: fragment size extraction, per-sample metrics (median length, peak position, short:long ratio, long-fragment fraction), three-group statistical comparisons (Kruskal-Wallis + pairwise Wilcoxon BH-corrected), continuous Spearman correlations vs ichorCNA tumour fraction, DELFI-style genome-wide fragmentation profiling (5Mb bins, LOESS GC correction, centromere exclusion).

## Key technical facts
- GISTIC2 uses Corrected_Copy_Number -> log2(CN/2) as input, not raw seg.median.logR — GISTIC validity is independent of ichorCNA curation.
- Mappability correction is absent from the QDNAseq pipeline (confirmed, documented as a stated limitation).
- B3_S12 excluded due to genuine ploidy non-identifiability (near-diploid vs near-triploid/WGD ambiguity unresolvable from coverage-only sWGS without BAF/orthogonal purity data) — not a pipeline failure.
- Sample B2_S11 retained in primary fragmentomics dataset despite suspected pre-analytical contamination (elevated duplicate rate, abnormal 300-400bp fragment tail, long-fragment fraction 59x cohort's next-highest); results reported both with and without it.

## Conventions
- Copy-first protocol: never modify original files in place; back up as `.ORIGINAL` before editing.
- Heredoc syntax (`cat > file << 'EOF'`) preferred for deploying scripts; Python `open()` when heredoc unsuitable.
- No decorative comment separators (`#=====` style) — plain `#` comments only.
- PNG preferred over PDF for figures.
- Script style: minimal changes to working originals (paths, sample lists, batch labels only).
- All BAMs/raw data live only on ALICE HPC (`/scratch/alice/n/nhsas1/PTCL/BAMs/` and `BAMs_batch3/`), excluded from git via .gitignore — never expected locally.

## Tools
samtools/1.18, Picard 3.0.0, R/4.3.1, QDNAseq, ichorCNA, RASCAL, GISTIC2 (Apptainer), GenVisR 1.31.1, biomaRt. SLURM scheduler on ALICE.
