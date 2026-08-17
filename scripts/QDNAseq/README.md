# QDNAseq

Genome-wide relative copy number profiling using QDNAseq (v1.36.0, Bioconductor 3.17).

All 34 plasma samples were profiled at 15 kb bin resolution with GC-content correction and CBS segmentation. Pre-built hg38 bin annotations were not available for this Bioconductor release; bins were constructed from BSgenome.Hsapiens.UCSC.hg38 (v1.4.5). QDNAseq output was retained as an independent cross-validation dataset against ichorCNA calls (Section 2.7 of thesis).

## Execution order

1. `install_packages.R` — install QDNAseq and dependencies (run once)
2. `QDNAseq_ALICE_automated.R` / `run_QDNAseq.slurm` — Batches 1+2 (n=21)
3. `QDNAseq_batch3.R` / `run_QDNAseq_batch3.slurm` — Batch 3 (n=13)
4. `generate_plots.R` / `run_generate_plots.slurm` — per-sample profile plots, Batches 1+2
5. `generate_plots_batch3.R` / `run_generate_plots_batch3.slurm` — Batch 3 plots
6. `export_rascal_input_batch3.R` / `run_export_rascal_batch3.slurm` — export segmented output for RASCAL input

## Inputs
BAM files (not included; patient-level data)

## Outputs
Per-sample RDS files and segmented copy number profiles (not included in repository; used internally for cross-validation against ichorCNA calls
