#!/bin/bash
#SBATCH --job-name=ichorCNA_all
#SBATCH --partition=short
#SBATCH --time=02:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=2
#SBATCH --array=11-21%4
#SBATCH --output=/scratch/alice/n/nhsas1/PTCL/ichorCNA/logs/ichorCNA_%a_%j.out
#SBATCH --error=/scratch/alice/n/nhsas1/PTCL/ichorCNA/logs/ichorCNA_%a_%j.err

# ============================================================
# ichorCNA — All 21 samples (SLURM array job)
# Includes B1_S01 — will skip if already done
# Each array task processes ONE sample independently
# ============================================================

echo "Array task:  $SLURM_ARRAY_TASK_ID"
echo "Job started: $(date)"
echo "Node:        $HOSTNAME"
echo ""

# --- Load R ---
source /etc/profile.d/modules.sh
module purge
module load R/4.3.1
export R_LIBS_USER=/home/n/nhsas1/R/library

# --- Shared paths ---
ICHORCNA_SCRIPT=/scratch/alice/n/nhsas1/PTCL/scripts/ichorCNA/scripts/runIchorCNA.R
WIG_DIR=/scratch/alice/n/nhsas1/PTCL/ichorCNA/wig
OUT_BASE=/scratch/alice/n/nhsas1/PTCL/ichorCNA/output
GC_WIG=/scratch/alice/n/nhsas1/PTCL/scripts/ichorCNA/inst/extdata/gc_hg38_1000kb.wig
MAP_WIG=/scratch/alice/n/nhsas1/PTCL/scripts/ichorCNA/inst/extdata/map_hg38_1000kb.wig
CENTROMERE=/scratch/alice/n/nhsas1/PTCL/scripts/ichorCNA/inst/extdata/GRCh38.GCA_000001405.2_centromere_acen.txt

# --- All 21 samples ---
declare -A SAMPLES
SAMPLES[1]="B1_S01"
SAMPLES[2]="B1_S02"
SAMPLES[3]="B1_S03"
SAMPLES[4]="B1_S04"
SAMPLES[5]="B1_S05"
SAMPLES[6]="B2_S01"
SAMPLES[7]="B2_S02"
SAMPLES[8]="B2_S03"
SAMPLES[9]="B2_S04"
SAMPLES[10]="B2_S05"
SAMPLES[11]="B2_S06"
SAMPLES[12]="B2_S07"
SAMPLES[13]="B2_S08"
SAMPLES[14]="B2_S09"
SAMPLES[15]="B2_S10"
SAMPLES[16]="B2_S11"
SAMPLES[17]="B2_S12"
SAMPLES[18]="B2_S13"
SAMPLES[19]="B2_S14"
SAMPLES[20]="B2_S15"
SAMPLES[21]="B2_S16"

# --- Get this task's sample ---
SAMPLE_ID=${SAMPLES[$SLURM_ARRAY_TASK_ID]}
WIG_FILE="$WIG_DIR/${SAMPLE_ID}.wig"
OUT_DIR="$OUT_BASE/${SAMPLE_ID}"

echo "Sample ID:  $SAMPLE_ID"
echo "WIG file:   $WIG_FILE"
echo "Output dir: $OUT_DIR"
echo ""

# --- Skip if already completed successfully ---
PARAMS_FILE="$OUT_DIR/${SAMPLE_ID}.params.txt"
if [ -f "$PARAMS_FILE" ]; then
    echo "params.txt already exists — $SAMPLE_ID already completed"
    echo "Delete $OUT_DIR to force rerun"
    echo "Tumour fraction from existing run:"
    grep "Tumor Fraction" "$PARAMS_FILE" | head -1
    exit 0
fi

# --- Check WIG file exists and is not empty ---
if [ ! -s "$WIG_FILE" ]; then
    echo "ERROR: WIG file missing or empty: $WIG_FILE"
    echo "Run readCounter first before ichorCNA"
    exit 1
fi

WIG_LINES=$(wc -l < "$WIG_FILE")
echo "WIG file: $WIG_LINES lines — OK"

# --- Create output directory ---
mkdir -p "$OUT_DIR"

# --- Run ichorCNA ---
echo ""
echo "Running ichorCNA for $SAMPLE_ID..."
echo "Started: $(date)"

Rscript "$ICHORCNA_SCRIPT" \
    --id "$SAMPLE_ID" \
    --WIG "$WIG_FILE" \
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
echo "ichorCNA finished: $(date)"
echo "Exit code: $EXIT_CODE"

# --- Verify output ---
echo ""
echo "=== Output files ==="
ls -lh "$OUT_DIR/" 2>/dev/null

if [ -f "$PARAMS_FILE" ]; then
    echo ""
    echo "=== RESULT: $SAMPLE_ID ==="
    grep -E "Tumor Fraction|Ploidy|GC-Map" "$PARAMS_FILE" | head -5
    echo "SUCCESS"
else
    echo "FAILED: params.txt not created — check .err log"
    exit 1
fi

echo ""
echo "Job finished: $(date)"
