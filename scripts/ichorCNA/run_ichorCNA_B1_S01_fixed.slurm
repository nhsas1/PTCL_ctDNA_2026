#!/bin/bash
#SBATCH --job-name=ichorCNA_B1_S01
#SBATCH --partition=short
#SBATCH --time=02:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=2
#SBATCH --output=/scratch/alice/n/nhsas1/PTCL/ichorCNA/logs/ichorCNA_B1_S01_%j.out
#SBATCH --error=/scratch/alice/n/nhsas1/PTCL/ichorCNA/logs/ichorCNA_B1_S01_%j.err

# ============================================================
# ichorCNA — Single Sample Test: B1_S01
# Calls runIchorCNA.R directly (no intermediate R wrapper)
# ============================================================

echo "Job started:  $(date)"
echo "Node:         $HOSTNAME"
echo "Job ID:       $SLURM_JOB_ID"
echo ""

# --- Load R ---
source /etc/profile.d/modules.sh
module purge
module load R/4.3.1

echo "R version: $(R --version | head -1)"
echo ""

# --- Set library path so R finds ichorCNA ---
export R_LIBS_USER=/home/n/nhsas1/R/library

# --- Paths ---
ICHORCNA_SCRIPT=/scratch/alice/n/nhsas1/PTCL/scripts/ichorCNA/scripts/runIchorCNA.R
WIG=/scratch/alice/n/nhsas1/PTCL/ichorCNA/wig/B1_S01.wig
GC_WIG=/scratch/alice/n/nhsas1/PTCL/scripts/ichorCNA/inst/extdata/gc_hg38_1000kb.wig
MAP_WIG=/scratch/alice/n/nhsas1/PTCL/scripts/ichorCNA/inst/extdata/map_hg38_1000kb.wig
CENTROMERE=/scratch/alice/n/nhsas1/PTCL/scripts/ichorCNA/inst/extdata/GRCh38.GCA_000001405.2_centromere_acen.txt
OUT_DIR=/scratch/alice/n/nhsas1/PTCL/ichorCNA/output/B1_S01

# --- Pre-flight checks ---
echo "=== Pre-flight checks ==="
for f in "$WIG" "$GC_WIG" "$MAP_WIG" "$CENTROMERE" "$ICHORCNA_SCRIPT"; do
    if [ -f "$f" ]; then
        echo "  FOUND:   $f"
    else
        echo "  MISSING: $f"
        echo "Aborting — required file missing"
        exit 1
    fi
done
echo ""

# --- Create output directory ---
mkdir -p "$OUT_DIR"
echo "Output directory: $OUT_DIR"
echo ""

# --- Run ichorCNA directly ---
echo "=== Running ichorCNA ==="
echo "Started: $(date)"
echo ""

Rscript "$ICHORCNA_SCRIPT" \
    --id "B1_S01" \
    --WIG "$WIG" \
    --gcWig "$GC_WIG" \
    --mapWig "$MAP_WIG" \
    --centromere "$CENTROMERE" \
    --ploidy "c(2,3)" \
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
    --outDir "$OUT_DIR"

EXIT_CODE=$?
echo ""
echo "ichorCNA exit code: $EXIT_CODE"
echo "Finished: $(date)"
echo ""

# --- Check outputs ---
echo "=== Output files ==="
ls -lh "$OUT_DIR/" 2>/dev/null || echo "Output directory is empty"
echo ""

# --- Show tumour fraction result ---
PARAMS="$OUT_DIR/B1_S01.params.txt"
if [ -f "$PARAMS" ]; then
    echo "=== TUMOUR FRACTION RESULT ==="
    cat "$PARAMS"
    echo ""
    echo "SUCCESS: ichorCNA completed for B1_S01"
else
    echo "FAILED: params.txt not found"
    echo "Check the .err log for the actual error message"
    exit 1
fi

echo ""
echo "Job finished: $(date)"
