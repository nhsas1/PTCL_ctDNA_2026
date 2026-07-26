#!/bin/bash
# Submits the Batch 3 QDNAseq pipeline (QDNAseq_batch3.R) for 13 samples.
#
# Resources match the Batch 1-2 job despite the smaller sample count, because Batch 3 BAMs
# are not deduplicated and so carry more reads per sample. cpus-per-task=4 is unused; the
# R script is single-threaded.
#
# See run_QDNAseq.slurm for the shared caveats: scripts/logs/ must pre-exist, the deployed
# script path is flat rather than nested, and the SUCCESS banner reflects only that the R
# script ran to completion, not that any sample was processed successfully.
#SBATCH --job-name=QDNAseq_B3
#SBATCH --output=/scratch/alice/n/nhsas1/PTCL/scripts/logs/QDNAseq_batch3_%j.out
#SBATCH --error=/scratch/alice/n/nhsas1/PTCL/scripts/logs/QDNAseq_batch3_%j.err
#SBATCH --time=08:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --partition=short
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=nhsas1@le.ac.uk

echo "========================================"
echo "QDNAseq PTCL Pipeline — Batch 3"
echo "Job ID:    $SLURM_JOB_ID"
echo "Node:      $SLURMD_NODENAME"
echo "Start:     $(date)"
echo "========================================"

source /etc/profile.d/modules.sh
module purge
module load R/4.3.1

echo "R version:"
R --version | head -1
echo "Which Rscript:"
which Rscript

Rscript /scratch/alice/n/nhsas1/PTCL/scripts/QDNAseq_batch3.R

EXIT_CODE=$?
echo "========================================"
echo "Finished: $(date)"
echo "Exit code: $EXIT_CODE"
if [ $EXIT_CODE -eq 0 ]; then
    echo "STATUS: SUCCESS"
else
    echo "STATUS: FAILED - check .err log"
fi
echo "========================================"
