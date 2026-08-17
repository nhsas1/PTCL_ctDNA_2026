# RASCAL

Absolute copy number scaling used as secondary cross-validation of ichorCNA tumour fraction estimates (RASCAL v0.7.0).

RASCAL fits absolute copy number by searching across ploidy and cellularity combinations to minimise deviation from integer copy number values. It was applied to high tumour fraction samples (TF >30%) where sufficient segment-level signal exists at this sequencing depth. An unconstrained search was run first; where results diverged from ichorCNA, a constrained search was performed with cellularity restricted to ichorCNA TF ±0.15 (thesis Section 2.11).

## Execution order

1. `export_rascal_input.R` / `run_export_rascal.slurm` — export ichorCNA segments in RASCAL-compatible format, Batches 1+2
2. `run_rascal_batch3.R` / `run_rascal_batch3.slurm` — unconstrained RASCAL search, Batch 3
3. `run_rascal_constrained_batch3.R` / `run_rascal_constrained_batch3.slurm` — TF-constrained search for divergent samples

## Inputs
Segmented ichorCNA output (`.seg.txt` files from `results/ichorCNA/seg_files/`)

## Outputs
Per-sample ploidy and cellularity estimates; expected TF comparison against ichorCNA
