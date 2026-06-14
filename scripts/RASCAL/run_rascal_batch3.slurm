#!/bin/bash
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
