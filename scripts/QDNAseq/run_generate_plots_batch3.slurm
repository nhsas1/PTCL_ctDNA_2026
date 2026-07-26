#!/bin/bash
# Submits generate_plots_batch3.R, producing four diagnostic QC plots per Batch 3 sample.
# Resources and caveats match run_generate_plots.slurm.
#SBATCH --job-name=plots_B3
#SBATCH --output=/scratch/alice/n/nhsas1/PTCL/scripts/logs/plots_batch3_%j.out
#SBATCH --error=/scratch/alice/n/nhsas1/PTCL/scripts/logs/plots_batch3_%j.err
#SBATCH --time=04:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --partition=short
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=nhsas1@le.ac.uk

source /etc/profile.d/modules.sh
module purge
module load R/4.3.1

echo "Starting Batch 3 plot generation: $(date)"
Rscript /scratch/alice/n/nhsas1/PTCL/scripts/generate_plots_batch3.R
echo "Finished: $(date)"
