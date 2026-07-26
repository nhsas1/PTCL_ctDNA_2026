#!/bin/bash
# Regenerates the dedup sensitivity and batch QC figures from the step 2 outputs.
# Does not recompute any metric.
#SBATCH --job-name=frag_figures
#SBATCH --partition=short
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:15:00
#SBATCH --output=/scratch/alice/n/nhsas1/PTCL/fragmentomics/logs/frag_figures_%j.out
#SBATCH --error=/scratch/alice/n/nhsas1/PTCL/fragmentomics/logs/frag_figures_%j.err

module load R/4.3.1

Rscript /scratch/alice/n/nhsas1/ptcl_ctdna_thesis/scripts/fragmentomics/generate_publication_figures.R
