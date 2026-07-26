#!/bin/bash
# Submits export_rascal_input_batch3.R, writing bin-level linear copy number for Batch 3.
#
# The job name says RASCAL, but these files are read by the QDNAseq-vs-ichorCNA
# concordance scripts, not by RASCAL. See the header of export_rascal_input_batch3.R.
#
# 4 hours and 32G for a pass over 13 BAMs that stops before segmentation.
# cpus-per-task=4 is unused.
#SBATCH --job-name=RASCAL_export_B3
#SBATCH --output=/scratch/alice/n/nhsas1/PTCL/scripts/logs/rascal_export_batch3_%j.out
#SBATCH --error=/scratch/alice/n/nhsas1/PTCL/scripts/logs/rascal_export_batch3_%j.err
#SBATCH --time=04:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --partition=short
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=nhsas1@le.ac.uk

source /etc/profile.d/modules.sh
module purge
module load R/4.3.1

echo "RASCAL Batch 3 input export starting: $(date)"
Rscript /scratch/alice/n/nhsas1/PTCL/scripts/export_rascal_input_batch3.R
echo "Finished: $(date)"
