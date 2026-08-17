# GISTIC2

Cohort-level recurrent copy number alteration detection using GISTIC2 (n=20 final CNA cohort).

GISTIC2 identifies genomic regions recurrently gained or lost across the cohort by comparing alteration frequency and amplitude against a background noise model, outputting FDR-corrected q-values. Input segment files use ichorCNA Corrected_Copy_Number values converted to log2(CN/2). Run via Apptainer container on ALICE. Full parameters are listed in Table 2 of the thesis.

## Execution order

1. `prepare_GISTIC_input.R` — prepare combined seg file from corrected ichorCNA output (n=20)
2. `run_GISTIC.slurm` — GISTIC2 execution via Apptainer
3. `plot_GISTIC.R` — significance landscape and frequency plots
4. `per_sample_GISTIC.R` — per-sample peak alteration tables

## Inputs
Corrected seg files from `results/ichorCNA/seg_files_corrected/` and `results/ichorCNA/seg_files/`; hg38 RefGene reference

## Outputs
`results/GISTIC/`: all_lesions, del_genes, amp_genes, frequency plots, input seg file
