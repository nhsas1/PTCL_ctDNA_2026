#!/bin/bash
#SBATCH --job-name=ptcl_bam_qc
#SBATCH --output=/scratch/alice/n/nhsas1/PTCL/QC_output/slurm_%j.log
#SBATCH --error=/scratch/alice/n/nhsas1/PTCL/QC_output/slurm_%j.err
#SBATCH --time=08:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --partition=short

module load samtools/1.18

echo "Samtools version: $(samtools --version | head -1)"
echo "Job started: $(date)"
echo "Running on node: $(hostname)"
echo ""

bash /scratch/alice/n/nhsas1/PTCL/scripts/bam_qc_final.sh

echo ""
echo "Job completed: $(date)"
