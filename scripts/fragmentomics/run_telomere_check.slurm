#!/bin/bash
# Step 5f: flags terminal bins in the per-bin results. See the R script header - this does
# not measure telomeric sequence content and cannot support a telomeric artefact argument.
#SBATCH --job-name=frag_telocheck
#SBATCH --partition=short
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:15:00
#SBATCH --output=/scratch/alice/n/nhsas1/PTCL/fragmentomics/logs/frag_telocheck_%j.out
#SBATCH --error=/scratch/alice/n/nhsas1/PTCL/fragmentomics/logs/frag_telocheck_%j.err

module load R/4.3.1
Rscript /scratch/alice/n/nhsas1/ptcl_ctdna_thesis/scripts/fragmentomics/step5f_telomere_check.R
