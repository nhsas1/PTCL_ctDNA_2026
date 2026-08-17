# FGA and Aneuploidy Score

Quantification of genome-wide copy number burden per sample using fraction of genome altered (FGA) and the Taylor Aneuploidy Score (AS; Taylor et al., 2018).

Applied to the final n=20 CNA cohort using ichorCNA corrected seg files. FGA was defined as the proportion of mappable autosomal genome bases within altered segments. AS was calculated at both 50% (AS_50) and 25% (AS_25) arm-coverage thresholds. Correlation between FGA/AS and GISTIC2 deletion burden was assessed as an internal consistency check (thesis Section 2.14).

## Execution order

1. `calculate_FGA_AS.R` — FGA and AS calculation for n=20 cohort
2. `calculate_AS_sensitivity.R` — AS_25 vs AS_50 threshold sensitivity comparison
3. `plot_FGA.R` — FGA distribution and TF correlation plots
4. `correlate_FGA_AS_GISTIC_n20.R` — Spearman correlation between FGA/AS and GISTIC2 deletion burden

## Inputs
Corrected ichorCNA seg files (`results/ichorCNA/seg_files/` and `seg_files_corrected/`)

## Outputs
`results/FGA_AS/`: summary tables, figures, GISTIC correlation plots and CSV
