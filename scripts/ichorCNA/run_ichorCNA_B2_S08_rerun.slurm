#!/bin/bash
# Diagnostic re-run of B2_S08 with the ploidy search started at 2.
#
# NOTE on --ploidy: this sets the starting values for ichorCNA's ploidy search, not a hard
# constraint. A run labelled "diploid" begins the search at 2 but the fitted ploidy is still
# estimated and can move away from it. Any argument resting on these runs needs the fitted
# output, not the requested parameter.
#
# EVIDENCE GAP: the params and seg files these alternative runs produced are not committed
# to the repository, and the captured logs echo the input parameters only. Their outcomes
# therefore cannot be checked from the repository as it stands.
#SBATCH --job-name=ichor_B2S08_rerun
#SBATCH --partition=short
#SBATCH --time=00:30:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --output=/scratch/alice/n/nhsas1/ptcl_ctdna_thesis/scripts/ichorCNA/logs/B2S08_rerun_%j.out
#SBATCH --error=/scratch/alice/n/nhsas1/ptcl_ctdna_thesis/scripts/ichorCNA/logs/B2S08_rerun_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=nhsas1@le.ac.uk

source /etc/profile.d/modules.sh
module purge
module load R/4.3.1

echo "Re-running ichorCNA for B2_S08 with constrained ploidy=2"
echo "Started: $(date)"

Rscript /scratch/alice/n/nhsas1/PTCL/scripts/ichorCNA/scripts/runIchorCNA.R \
    --id "B2_S08_ploidy2only" \
    --WIG "/scratch/alice/n/nhsas1/PTCL/ichorCNA/wig/B2_S08.wig" \
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
    --outDir "/scratch/alice/n/nhsas1/PTCL/ichorCNA/output/B2_S08_ploidy2only"

echo "Finished: $(date)"
