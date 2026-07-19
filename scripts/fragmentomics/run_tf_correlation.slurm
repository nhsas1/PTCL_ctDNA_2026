#!/bin/bash
#SBATCH --job-name=frag_tfcorr
#SBATCH --partition=short
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:20:00
#SBATCH --output=/scratch/alice/n/nhsas1/PTCL/fragmentomics/logs/frag_tfcorr_%j.out
#SBATCH --error=/scratch/alice/n/nhsas1/PTCL/fragmentomics/logs/frag_tfcorr_%j.err

module load R/4.3.1

Rscript /scratch/alice/n/nhsas1/ptcl_ctdna_thesis/scripts/fragmentomics/step4_tf_correlation.R
