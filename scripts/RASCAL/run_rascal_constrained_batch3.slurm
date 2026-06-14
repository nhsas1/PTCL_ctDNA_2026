#!/bin/bash
#SBATCH --job-name=RASCAL_con_B3
#SBATCH --partition=short
#SBATCH --time=01:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=2
#SBATCH --output=/scratch/alice/n/nhsas1/PTCL/scripts/RASCAL/logs/rascal_con_batch3_%j.out
#SBATCH --error=/scratch/alice/n/nhsas1/PTCL/scripts/RASCAL/logs/rascal_con_batch3_%j.err

echo "Job started: $(date)"
echo "Node: $(hostname)"
echo "Job ID: $SLURM_JOB_ID"

source /etc/profile.d/modules.sh
module purge
module load R/4.3.1

Rscript /scratch/alice/n/nhsas1/PTCL/scripts/RASCAL/run_rascal_constrained_batch3.R

echo "Job finished: $(date)"
