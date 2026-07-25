#!/bin/bash
#SBATCH --job-name=ichor_B2S08_noSub
#SBATCH --partition=short
#SBATCH --time=00:30:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --output=/scratch/alice/n/nhsas1/ptcl_ctdna_thesis/scripts/ichorCNA/logs/B2S08_noSub_%j.out
#SBATCH --error=/scratch/alice/n/nhsas1/ptcl_ctdna_thesis/scripts/ichorCNA/logs/B2S08_noSub_%j.err

source /etc/profile.d/modules.sh
module purge
module load R/4.3.1

mkdir -p /scratch/alice/n/nhsas1/PTCL/ichorCNA/output/B2_S08_noSubclone

Rscript /scratch/alice/n/nhsas1/PTCL/scripts/ichorCNA/scripts/runIchorCNA.R \
    --id "B2_S08_noSubclone" \
    --WIG "/scratch/alice/n/nhsas1/PTCL/ichorCNA/wig/B2_S08.wig" \
    --gcWig "/scratch/alice/n/nhsas1/PTCL/scripts/ichorCNA/inst/extdata/gc_hg38_1000kb.wig" \
    --mapWig "/scratch/alice/n/nhsas1/PTCL/scripts/ichorCNA/inst/extdata/map_hg38_1000kb.wig" \
    --centromere "/scratch/alice/n/nhsas1/PTCL/scripts/ichorCNA/inst/extdata/GRCh38.GCA_000001405.2_centromere_acen.txt" \
    --ploidy "c(2,3)" \
    --normal "c(0.5,0.6,0.7,0.8,0.9,0.95,0.99)" \
    --maxCN 5 \
    --includeHOMD False \
    --chrs "c(1:22)" \
    --chrTrain "c(1:22)" \
    --normalizeMaleX True \
    --estimateNormal True \
    --estimatePloidy True \
    --estimateScPrevalence False \
    --scStates "c()" \
    --txnE 0.9999 \
    --txnStrength 10000 \
    --plotFileType pdf \
    --plotYLim "c(-2,2)" \
    --outDir "/scratch/alice/n/nhsas1/PTCL/ichorCNA/output/B2_S08_noSubclone"

echo "Done: $(date)"
