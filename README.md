# Shallow Whole Genome Sequencing of ctDNA in T-Cell Lymphoma

**Thesis:** Shallow whole genome sequencing of ctDNA in T-cell lymphoma: copy number alterations and fragmentomics  
**Student:** Noor Shaban
**Institution:** University of Leicester - MSc Bioinformatics - 



## Overview

Analysis pipeline for plasma sWGS ctDNA profiling in peripheral T-cell lymphoma (PTCL), across two analytical arms:

- **CNA profiling** : tumour fraction estimation (ichorCNA), manual curation, absolute copy number scaling (RASCAL), recurrent CNA detection (GISTIC2), arm-level aneuploidy scoring (FGA/AS)
- **Fragmentomics** :  fragment-length metrics, tumour fraction correlation

A sub-analysis compares plasma ctDNA against matched FFPE tissue sWGS for 15 evaluable patient pairs.

Raw BAM files are not included (patient-level data).



## Cohort

| | |
|---|---|
| Total plasma samples | 34 (Batch 1: n=5, Batch 2: n=16, Batch 3: n=13) |
| Sequencing | Illumina NovaSeq, GRCh38 |
| CNA cohort | n=20 (TF ≥3%; 13 below detection floor, 1 excluded for unresolvable ploidy ambiguity) |
| Fragmentomics cohort | n=33 (B2_S11 excluded: outlier across fragment-length metrics) |
| Matched FFPE pairs | 16 total; 14 evaluable (1 excluded: degraded library |



## Repository Structure

```
PTCL_ctDNA_2026/
├── scripts/
│   ├── QDNAseq/         # Relative copy number profiling
│   ├── ichorCNA/        # Tumour fraction estimation and CNA calling
│   ├── RASCAL/          # Absolute copy number scaling
│   ├── GISTIC/          # Cohort-level recurrent CNA detection
│   ├── FGA_AS/          # Fraction genome altered and aneuploidy scoring
│   ├── GenVisR/         # Cohort CNA frequency visualisation
│   └── fragmentomics/   # Fragment-length metrics
├── results/             # Summary tables and figures per analysis arm
├── QC_results/          # SAMtools flagstat output for all 34 plasma and 16 FFPE BAMs
├── data/
│   └── PTCL_Gene_List.xlsx  # Curated PTCL driver gene list (literature + cBioPortal)
└── README.md
```



## Dependencies

| Tool | Version |
|---|---|
| R | 4.3.1 |
| QDNAseq | 1.36.0 (Bioconductor 3.17) |
| ichorCNA | broadinstitute/ichorCNA (GitHub) |
| RASCAL | 0.7.0 (crukci-bioinformatics/rascal) |
| GISTIC2 | Apptainer container |
| SAMtools | 1.18 |
| Picard MarkDuplicates | 2.18.23 |
| HMMcopy readCounter | Built from source (cmake 3.27.7) |

All scripts were run on the University of Leicester ALICE HPC cluster (Rocky Linux 9, SLURM). SLURM scripts contain ALICE-specific paths and will require updating for other environments.



## Execution Order

### CNA Arm
1. `scripts/QDNAseq/run_QDNAseq.R` — relative copy number profiling (15 kb bins)
2. `scripts/ichorCNA/run_readCounter_*.slurm` → `run_ichorCNA_*.slurm` — TF estimation, all 34 samples
3. `scripts/ichorCNA/ichorCNA_anomaly_scan.R` — flag sign-discordant segments for manual review
4. `scripts/RASCAL/run_rascal_constrained.R` — absolute copy number scaling (n=20)
5. `scripts/GISTIC/prepare_GISTIC_input.R` → `run_GISTIC.slurm` → `plot_GISTIC.R`
6. `scripts/FGA_AS/calculate_FGA_AS.R` → `correlate_FGA_AS_GISTIC_n20.R`
7. `scripts/GenVisR/plot_cnFreq_cohort.R` — cohort CNA frequency plot

### Fragmentomics Arm
1. `scripts/fragmentomics/run_batch3_markdup.slurm` : Picard MarkDuplicates (Batch 3 only)
2. `scripts/fragmentomics/run_fragment_extraction.slurm` : fragment-length histograms
3. `scripts/fragmentomics/calculate_fragment_metrics.R` : per-sample metrics
4. `scripts/fragmentomics/tf_correlation.R` : Spearman TF correlation, BH-corrected



## Key Results

**CNA:** ichorCNA detected tumour signal in 20/34 samples (TF range 3.6–62.3%). GISTIC2 identified three significant deletion peaks: 4q35.1 (*FAT1*, *IRF2*; q=3.9×10⁻⁸), 1p35.3 (q=0.007), and 1q43 (q=0.045). No amplification peaks reached significance. These deletions are consistent with established tumour suppressor loss in PTCL.

**Fragmentomics:** Three metrics correlated significantly with tumour fraction (Spearman, BH-corrected): median fragment length (ρ=−0.515, p=0.0022), short:long ratio (ρ=+0.479, p=0.0053), and peak length (ρ=−0.459, p=0.0073), supporting the utility of fragment-length profiling as a non-invasive tumour burden proxy.

**FFPE:** 9/14 evaluable pairs concordant (64%); tissue sWGS recovered ctDNA signal in 4/5 patients below the plasma detection floor.



## Important Notes

Three ichorCNA samples required manual curation (B2_S08, B2_S15, B3_S09) and one was excluded for unresolvable ploidy ambiguity (B3_S12). Batch 3 BAMs were supplied without duplicate marking; Picard MarkDuplicates was applied post-hoc for the fragmentomics arm only.



## Data Availability

Raw BAM files are not publicly available (patient-level genomic data). The curated PTCL gene list (`data/PTCL_Gene_List.xlsx`) is compiled from published literature and cBioPortal.



## Licence

Code released under the MIT Licence. Raw sequencing data and patient metadata are not included.
