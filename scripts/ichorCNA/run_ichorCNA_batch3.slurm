#!/bin/bash
#SBATCH --job-name=ichorCNA_B3
#SBATCH --partition=short
#SBATCH --time=02:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=2
#SBATCH --array=1-13%4
#SBATCH --output=/scratch/alice/n/nhsas1/PTCL/ichorCNA/logs/ichorCNA_B3_%a_%j.out
#SBATCH --error=/scratch/alice/n/nhsas1/PTCL/ichorCNA/logs/ichorCNA_B3_%a_%j.err

# ============================================================
# ichorCNA — Batch 3 (13 new samples)
# Each array task processes ONE sample independently
# Run AFTER run_readCounter_batch3.slurm completes
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
WIG_DIR=/scratch/alice/n/nhsas1/PTCL/ichorCNA/wig_batch3
OUT_BASE=/scratch/alice/n/nhsas1/PTCL/ichorCNA/output_batch3
GC_WIG=/scratch/alice/n/nhsas1/PTCL/scripts/ichorCNA/inst/extdata/gc_hg38_1000kb.wig
MAP_WIG=/scratch/alice/n/nhsas1/PTCL/scripts/ichorCNA/inst/extdata/map_hg38_1000kb.wig
CENTROMERE=/scratch/alice/n/nhsas1/PTCL/scripts/ichorCNA/inst/extdata/GRCh38.GCA_000001405.2_centromere_acen.txt

# --- Sample list ---
declare -A SAMPLES
SAMPLES[1]="B3_S01"
SAMPLES[2]="B3_S02"
SAMPLES[3]="B3_S03"
SAMPLES[4]="B3_S04"
SAMPLES[5]="B3_S05"
SAMPLES[6]="B3_S06"
SAMPLES[7]="B3_S07"
SAMPLES[8]="B3_S08"
SAMPLES[9]="B3_S09"
SAMPLES[10]="B3_S10"
SAMPLES[11]="B3_S11"
SAMPLES[12]="B3_S12"
SAMPLES[13]="B3_S13"

# --- Get this task's sample ---
SAMPLE_ID=${SAMPLES[$SLURM_ARRAY_TASK_ID]}
WIG_FILE="$WIG_DIR/${SAMPLE_ID}.wig"
OUT_DIR="$OUT_BASE/${SAMPLE_ID}"

echo "Sample ID:  $SAMPLE_ID"
echo "WIG file:   $WIG_FILE"
echo "Output dir: $OUT_DIR"
echo ""

# --- Skip if already completed ---
PARAMS_FILE="$OUT_DIR/${SAMPLE_ID}.params.txt"
if [ -f "$PARAMS_FILE" ]; then
    echo "params.txt already exists — $SAMPLE_ID already completed"
    echo "Delete $OUT_DIR to force rerun"
    echo "Tumour fraction from existing run:"
    grep "Tumor Fraction" "$PARAMS_FILE" | head -1
    exit 0
fi

# --- Check WIG file exists ---
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
