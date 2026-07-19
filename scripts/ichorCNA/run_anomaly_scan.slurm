#!/bin/bash
#SBATCH --job-name=ichor_anomaly
#SBATCH --partition=short
#SBATCH --time=00:20:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --output=/scratch/alice/n/nhsas1/ptcl_ctdna_thesis/scripts/ichorCNA/logs/anomaly_scan_%j.out
#SBATCH --error=/scratch/alice/n/nhsas1/ptcl_ctdna_thesis/scripts/ichorCNA/logs/anomaly_scan_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=nhsas1@le.ac.uk

echo "Job started: $(date)"
echo "Node: $(hostname)"
echo "Job ID: $SLURM_JOB_ID"
echo ""

source /etc/profile.d/modules.sh
module purge
module load R/4.3.1

Rscript /scratch/alice/n/nhsas1/ptcl_ctdna_thesis/scripts/ichorCNA/ichorCNA_anomaly_scan.R

EXIT_CODE=$?
echo ""
echo "Job finished: $(date)"
echo "Exit code: $EXIT_CODE"
if [ $EXIT_CODE -eq 0 ]; then
    echo "STATUS: SUCCESS"
else
    echo "STATUS: FAILED — check .err log"
fi
