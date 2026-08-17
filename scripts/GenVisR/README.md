# GenVisR

Cohort-level copy number alteration visualisation using GenVisR (v1.31.1, Bioconductor 3.17).

Generates genome-wide copy number frequency plots stratified by PTCL subtype (cnFreq) and per-sample copy number spectrum plots (cnSpec) for the final n=20 CNA cohort (thesis Section 2.15).

## Execution order

1. `plot_cnFreq_cohort.R` — cohort CNA frequency plot and per-sample spectrum

## Inputs
Corrected ichorCNA seg files (`results/ichorCNA/seg_files/` and `seg_files_corrected/`)

## Outputs
`results/GenVisR/`: genome-wide frequency plots by subtype, per-sample copy number spectrum
