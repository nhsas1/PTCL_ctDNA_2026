# ichorCNA

Tumour fraction estimation and genome-wide HMM copy number calling at 1 Mb bin resolution (hg38).

ichorCNA was selected for its validated performance at ultra-low-pass cfDNA coverage depths (0.1–1×) and its established 3% tumour fraction detection floor (Adalsteinsson et al., 2017). All 34 plasma samples and 16 matched FFPE tissue samples were processed. Three samples required manual curation after automated anomaly screening; one was excluded for unresolvable ploidy ambiguity (see thesis Section 2.10).

## Installation
Run `install_ichorCNA.R` interactively before first use. See script for full dependency list and HMMcopy readCounter build instructions.

## Execution order

1. `run_readCounter_all.slurm` — bin BAM reads into 1 Mb WIG files, Batches 1+2
2. `run_readCounter_batch3.slurm` — same for Batch 3
3. `run_ichorCNA_all.slurm` — HMM copy number calling, Batches 1+2
4. `run_ichorCNA_batch3.slurm` — same for Batch 3
5. `ichorCNA_anomaly_scan.R` / `run_anomaly_scan.slurm` — automated detection of sign-discordant segments
6. Manual curation — B2_S08, B2_S15, B3_S09 corrected; B3_S12 excluded (see `results/ichorCNA/PTCL_anomaly_corrections_evidence.tsv`)
7. `plot_QDNAseq_ichorCNA_concordance.R` / `run_QDNAseq_ichorCNA_concordance.slurm` — cross-tool concordance validation
8. `plot_bin_size_comparison.R` — 1 Mb vs 500 kb bin-size sensitivity analysis

## Inputs
WIG files from readCounter; GC and mappability reference WIG files (hg38, 1 Mb; bundled with ichorCNA package)

## Outputs
Per-sample: `.seg.txt`, `.params.txt`, genome-wide plots ( not provided here ) 
Corrected seg files: `results/ichorCNA/seg_files_corrected/`
