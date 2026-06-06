#!/bin/bash
#SBATCH --job-name=ptcl_bam_qc_batch3
#SBATCH --output=/scratch/alice/n/nhsas1/PTCL/QC_output_batch3/slurm_%j.log
#SBATCH --error=/scratch/alice/n/nhsas1/PTCL/QC_output_batch3/slurm_%j.err
#SBATCH --time=08:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --partition=short

# Create output directory before SLURM tries to write logs there
mkdir -p /scratch/alice/n/nhsas1/PTCL/QC_output_batch3

source /etc/profile.d/modules.sh
module purge
module load samtools/1.18

echo "Samtools version: $(samtools --version | head -1)"
echo "Job started: $(date)"
echo "Running on node: $(hostname)"
echo ""

bash /scratch/alice/n/nhsas1/PTCL/scripts/bam_qc_batch3.sh

echo ""
echo "Job completed: $(date)"
