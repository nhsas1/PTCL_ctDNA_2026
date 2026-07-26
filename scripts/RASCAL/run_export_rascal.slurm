#!/bin/bash
# Runs export_rascal_input.R for Batches 1-2. Despite the job name, the output feeds the
# QDNAseq-vs-ichorCNA concordance analysis, not RASCAL; see the R script header.
#
# 4 hours and 32G for a pass over 21 BAMs that stops before segmentation.
# cpus-per-task is requested but unused, as the R script is single-threaded.
#SBATCH --job-name=RASCAL_export
#SBATCH --output=/scratch/alice/n/nhsas1/PTCL/scripts/logs/rascal_export_%j.out
#SBATCH --error=/scratch/alice/n/nhsas1/PTCL/scripts/logs/rascal_export_%j.err
#SBATCH --time=04:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --partition=short
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=nhsas1@le.ac.uk

source /etc/profile.d/modules.sh
module purge
module load R/4.3.1

echo "RASCAL input export starting: $(date)"
Rscript /scratch/alice/n/nhsas1/PTCL/scripts/export_rascal_input.R
echo "Finished: $(date)"
