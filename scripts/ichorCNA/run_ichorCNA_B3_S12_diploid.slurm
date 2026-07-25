#!/bin/bash
#SBATCH --job-name=ichor_B3S12_diploid
#SBATCH --partition=short
#SBATCH --time=00:30:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --output=/scratch/alice/n/nhsas1/ptcl_ctdna_thesis/scripts/ichorCNA/logs/B3S12_diploid_%j.out
#SBATCH --error=/scratch/alice/n/nhsas1/ptcl_ctdna_thesis/scripts/ichorCNA/logs/B3S12_diploid_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=nhsas1@le.ac.uk

source /etc/profile.d/modules.sh
module purge
module load R/4.3.1

echo "Job started: $(date)"
echo "Node: $(hostname)"
echo "Job ID: $SLURM_JOB_ID"

mkdir -p /scratch/alice/n/nhsas1/PTCL/ichorCNA/output_batch3/B3_S12_diploid

Rscript /scratch/alice/n/nhsas1/PTCL/scripts/ichorCNA/scripts/runIchorCNA.R \
    --id "B3_S12_diploid" \
    --WIG "/scratch/alice/n/nhsas1/PTCL/ichorCNA/wig_batch3/B3_S12.wig" \
    --gcWig "/scratch/alice/n/nhsas1/PTCL/scripts/ichorCNA/inst/extdata/gc_hg38_1000kb.wig" \
    --mapWig "/scratch/alice/n/nhsas1/PTCL/scripts/ichorCNA/inst/extdata/map_hg38_1000kb.wig" \
    --centromere "/scratch/alice/n/nhsas1/PTCL/scripts/ichorCNA/inst/extdata/GRCh38.GCA_000001405.2_centromere_acen.txt" \
    --ploidy "c(2)" \
    --normal "c(0.5,0.6,0.7,0.8,0.9,0.95,0.99)" \
    --maxCN 5 \
    --includeHOMD False \
    --chrs "c(1:22)" \
    --chrTrain "c(1:22)" \
    --normalizeMaleX True \
    --estimateNormal True \
    --estimatePloidy True \
    --estimateScPrevalence True \
    --scStates "c(1,3)" \
    --txnE 0.9999 \
    --txnStrength 10000 \
    --plotFileType pdf \
    --plotYLim "c(-2,2)" \
    --outDir "/scratch/alice/n/nhsas1/PTCL/ichorCNA/output_batch3/B3_S12_diploid"

echo "Job finished: $(date)"
