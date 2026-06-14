#!/bin/bash
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
