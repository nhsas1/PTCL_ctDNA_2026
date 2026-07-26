#!/bin/bash
# Diagnostic re-run of B3_S12 with subclonal states disabled, testing whether the high
# fitted ploidy is an artefact of the subclonal component absorbing signal.
#
# NOTE on --ploidy: this sets the starting values for ichorCNA's ploidy search, not a hard
# constraint. A run labelled "diploid" begins the search at 2 but the fitted ploidy is still
# estimated and can move away from it. Any argument resting on these runs needs the fitted
# output, not the requested parameter.
#
# EVIDENCE GAP: the params and seg files these alternative runs produced are not committed
# to the repository, and the captured logs echo the input parameters only. Their outcomes
# therefore cannot be checked from the repository as it stands.
#SBATCH --job-name=ichor_B3S12_noSub
#SBATCH --partition=short
#SBATCH --time=00:30:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --output=/scratch/alice/n/nhsas1/ptcl_ctdna_thesis/scripts/ichorCNA/logs/B3S12_noSub_%j.out
#SBATCH --error=/scratch/alice/n/nhsas1/ptcl_ctdna_thesis/scripts/ichorCNA/logs/B3S12_noSub_%j.err

source /etc/profile.d/modules.sh
module purge
module load R/4.3.1

mkdir -p /scratch/alice/n/nhsas1/PTCL/ichorCNA/output_batch3/B3_S12_noSubclone

Rscript /scratch/alice/n/nhsas1/PTCL/scripts/ichorCNA/scripts/runIchorCNA.R \
    --id "B3_S12_noSubclone" \
    --WIG "/scratch/alice/n/nhsas1/PTCL/ichorCNA/wig_batch3/B3_S12.wig" \
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
    --outDir "/scratch/alice/n/nhsas1/PTCL/ichorCNA/output_batch3/B3_S12_noSubclone"

echo "Done: $(date)"
