#!/bin/bash
# Step 5c: diagnostic report on the most extreme bins, the evidence behind the centromere
# exclusion in step 5d. Produces no figures - the chromosome ordering and palette fixes it
# used to apply now live in step 5 itself, so there is no longer a filename collision.
#SBATCH --job-name=frag_fixplot
#SBATCH --partition=short
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:20:00
#SBATCH --output=/scratch/alice/n/nhsas1/PTCL/fragmentomics/logs/frag_fixplot_%j.out
#SBATCH --error=/scratch/alice/n/nhsas1/PTCL/fragmentomics/logs/frag_fixplot_%j.err

module load R/4.3.1

Rscript /scratch/alice/n/nhsas1/ptcl_ctdna_thesis/scripts/fragmentomics/step5c_fix_profile_plot.R
