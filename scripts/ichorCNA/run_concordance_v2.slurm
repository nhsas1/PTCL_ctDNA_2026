#!/bin/bash
# Runs the current QDNAseq-vs-ichorCNA concordance analysis, including the Bland-Altman
# agreement panel. This is the version that produced the committed figures and statistics.
#SBATCH --job-name=concordance_v2
#SBATCH --partition=short
#SBATCH --time=00:30:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --output=/scratch/alice/n/nhsas1/ptcl_ctdna_thesis/scripts/ichorCNA/logs/concordance_v2_%j.out
#SBATCH --error=/scratch/alice/n/nhsas1/ptcl_ctdna_thesis/scripts/ichorCNA/logs/concordance_v2_%j.err

source /etc/profile.d/modules.sh
module purge
module load R/4.3.1
export R_LIBS_USER=/home/n/nhsas1/R/library

SCRIPT=/scratch/alice/n/nhsas1/ptcl_ctdna_thesis/scripts/ichorCNA/plot_qdnaseq_ichorcna_concordance_v2.R

echo "Job started: $(date)"

if [ ! -f "$SCRIPT" ]; then
    echo "ERROR: R script not found at $SCRIPT"
    exit 1
fi

Rscript "$SCRIPT"

echo "Exit code: $?"
echo "Job finished: $(date)"
