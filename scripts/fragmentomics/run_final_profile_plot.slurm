#!/bin/bash
# Step 5d: excludes centromeric bins, recomputes within-sample z-scores over the retained
# bins, and writes the final profile figures and table.
#SBATCH --job-name=frag_finalplot
#SBATCH --partition=short
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:20:00
#SBATCH --output=/scratch/alice/n/nhsas1/PTCL/fragmentomics/logs/frag_finalplot_%j.out
#SBATCH --error=/scratch/alice/n/nhsas1/PTCL/fragmentomics/logs/frag_finalplot_%j.err

module load R/4.3.1

Rscript /scratch/alice/n/nhsas1/ptcl_ctdna_thesis/scripts/fragmentomics/step5d_final_profile_plot.R
