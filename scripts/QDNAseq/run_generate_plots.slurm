#!/bin/bash
#SBATCH --job-name=QDNAseq_plots
#SBATCH --output=/scratch/alice/n/nhsas1/PTCL/scripts/logs/plots_%j.out
#SBATCH --error=/scratch/alice/n/nhsas1/PTCL/scripts/logs/plots_%j.err
#SBATCH --time=04:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --partition=short
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=nhsas1@le.ac.uk

source /etc/profile.d/modules.sh
module purge
module load R/4.3.1

echo "Starting plot generation: $(date)"
Rscript /scratch/alice/n/nhsas1/PTCL/scripts/generate_plots.R
echo "Finished: $(date)"
