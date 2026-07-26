#!/bin/bash
# ARCHIVED - do not run. Kept for provenance only.
#
# This submits the original QDNAseq-vs-ichorCNA concordance script, which writes the same
# output filenames as version 2 and would replace the committed figures with an older
# two-panel version lacking the Bland-Altman agreement analysis.
#
# Use scripts/ichorCNA/run_concordance_v2.slurm instead. The SCRIPT path below still points
# at the original deployed location and is deliberately not updated, since running this is
# not intended.
#SBATCH --job-name=qdna_ichor_concordance
#SBATCH --partition=short
#SBATCH --time=00:30:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --output=/scratch/alice/n/nhsas1/ptcl_ctdna_thesis/scripts/ichorCNA/logs/concordance_%j.out
#SBATCH --error=/scratch/alice/n/nhsas1/ptcl_ctdna_thesis/scripts/ichorCNA/logs/concordance_%j.err

source /etc/profile.d/modules.sh
module purge
module load R/4.3.1
export R_LIBS_USER=/home/n/nhsas1/R/library

SCRIPT=/scratch/alice/n/nhsas1/ptcl_ctdna_thesis/scripts/ichorCNA/plot_qdnaseq_ichorcna_concordance.R

echo "Job started: $(date)"
echo "Running concordance analysis"

if [ ! -f "$SCRIPT" ]; then
    echo "ERROR: R script not found at $SCRIPT"
    exit 1
fi

Rscript "$SCRIPT"

echo "Exit code: $?"
echo "Job finished: $(date)"
