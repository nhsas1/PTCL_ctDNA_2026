#!/bin/bash
# Runs the unconstrained RASCAL fit for all 13 Batch 3 samples.
#
# 3 hours and 32G cover a ploidy-cellularity grid search per sample; the grid is the cost,
# not the input size. cpus-per-task=2 is unused because the search is single-threaded.
#
# Unlike the QDNAseq wrappers this writes its logs to scripts/RASCAL/logs/ and sets no
# mail notification, so a failure here is silent. The directory must already exist.
#
# The R script catches every per-sample error and never calls stop, so the job exits zero
# even if every sample failed. Check the summary CSV rather than the exit status.
#SBATCH --job-name=RASCAL_B3
#SBATCH --partition=short
#SBATCH --time=03:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=2
#SBATCH --output=/scratch/alice/n/nhsas1/PTCL/scripts/RASCAL/logs/rascal_batch3_%j.out
#SBATCH --error=/scratch/alice/n/nhsas1/PTCL/scripts/RASCAL/logs/rascal_batch3_%j.err

echo "Job started: $(date)"
echo "Node: $(hostname)"
echo "Job ID: $SLURM_JOB_ID"

source /etc/profile.d/modules.sh
module purge
module load R/4.3.1

Rscript /scratch/alice/n/nhsas1/PTCL/scripts/RASCAL/run_rascal_batch3.R

echo "Job finished: $(date)"
