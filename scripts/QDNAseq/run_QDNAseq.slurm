#!/bin/bash
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
