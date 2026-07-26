#!/bin/bash
# Step 5: DELFI-style genome-wide fragmentation profile, 5Mb bins with LOESS GC correction.
#SBATCH --job-name=frag_gwprofile
#SBATCH --partition=short
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --output=/scratch/alice/n/nhsas1/PTCL/fragmentomics/logs/frag_gwprofile_%j.out
#SBATCH --error=/scratch/alice/n/nhsas1/PTCL/fragmentomics/logs/frag_gwprofile_%j.err

module load R/4.3.1

Rscript /scratch/alice/n/nhsas1/ptcl_ctdna_thesis/scripts/fragmentomics/step5_genomewide_profile.R
