#!/bin/bash
# Submits the Batch 1-2 QDNAseq pipeline (QDNAseq_ALICE_automated.R) for all 21 samples.
#
# 8 hours and 32G cover a serial pass over 21 BAMs; the memory is driven by holding the
# 15kb genome-wide bin object alongside per-sample read counts. cpus-per-task=4 is
# requested but not used - no BiocParallel backend is registered and the R script is
# single-threaded, so this only lengthens the queue wait.
#
# The log directory must already exist; SLURM will reject the job at submission if
# scripts/logs/ is missing. Note also that the Rscript path below is flat
# (scripts/QDNAseq_ALICE_automated.R) whereas the repository nests these under
# scripts/QDNAseq/, so the deployed layout on ALICE differs from the repository layout.
#
# CAVEAT on the SUCCESS banner: the R script wraps every sample in tryCatch and never
# calls stop(), so it exits 0 even if every sample failed. A green STATUS here means the
# script ran to completion, not that any sample succeeded - check the .out log for the
# per-sample SUCCESS lines.
#SBATCH --job-name=QDNAseq_PTCL
#SBATCH --output=/scratch/alice/n/nhsas1/PTCL/scripts/logs/QDNAseq_%j.out
#SBATCH --error=/scratch/alice/n/nhsas1/PTCL/scripts/logs/QDNAseq_%j.err
#SBATCH --time=08:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --partition=short
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=nhsas1@le.ac.uk

echo "========================================"
echo "QDNAseq PTCL Pipeline"
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

Rscript /scratch/alice/n/nhsas1/PTCL/scripts/QDNAseq_ALICE_automated.R

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
