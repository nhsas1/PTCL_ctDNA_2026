#!/bin/bash
# Submits generate_plots.R, producing four diagnostic QC plots per Batch 1-2 sample.
#
# 4 hours is shorter than the main pipeline because plotting stops before the seg export,
# but memory stays at 32G since this script re-runs the full pipeline from BAM and holds
# the same objects. cpus-per-task=4 is unused.
#
# No exit code is captured here, so a failed R run is only visible in the .err log.
#SBATCH --job-name=QDNAseq_plots
#SBATCH --output=/scratch/alice/n/nhsas1/PTCL/scripts/logs/plots_%j.out
#SBATCH --error=/scratch/alice/n/nhsas1/PTCL/scripts/logs/plots_%j.err
#SBATCH --time=04:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --partition=short
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=nhsas1@le.ac.uk

source /etc/profile.d/modules.sh
module purge
module load R/4.3.1

echo "Starting plot generation: $(date)"
Rscript /scratch/alice/n/nhsas1/PTCL/scripts/generate_plots.R
echo "Finished: $(date)"
