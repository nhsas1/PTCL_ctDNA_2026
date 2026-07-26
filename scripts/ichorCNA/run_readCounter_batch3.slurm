#!/bin/bash
# Batch 3 counterpart of run_readCounter_all.slurm, producing 1Mb WIG files for the 13
# Batch 3 samples. Same bin size and quality settings, so counts remain comparable across
# batches, with the caveat that Batch 3 BAMs retain their duplicates.
#SBATCH --job-name=readCounter_B3
#SBATCH --partition=short
#SBATCH --time=01:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --array=1-13%4
#SBATCH --output=/scratch/alice/n/nhsas1/PTCL/ichorCNA/logs/readCounter_B3_%a_%j.out
#SBATCH --error=/scratch/alice/n/nhsas1/PTCL/ichorCNA/logs/readCounter_B3_%a_%j.err

# ============================================================
# readCounter — Batch 3 (13 new samples)
# Each array task processes ONE sample independently
# ============================================================

echo "Array task:  $SLURM_ARRAY_TASK_ID"
echo "Job started: $(date)"
echo "Node:        $HOSTNAME"
echo ""

# --- Tool and shared paths ---
READCOUNTER=/scratch/alice/n/nhsas1/PTCL/hmmcopy_utils/bin/readCounter
BAM_DIR=/scratch/alice/n/nhsas1/PTCL/BAMs_batch3
WIG_DIR=/scratch/alice/n/nhsas1/PTCL/ichorCNA/wig_batch3
CHROMOSOMES="chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22"
WINDOW=1000000
QUALITY=20

# --- Create WIG output directory if not present ---
mkdir -p "$WIG_DIR"

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

declare -A BAMS
BAMS[1]="LB_20241212_1_SI.bam"
BAMS[2]="LB_20241212_2_SI.bam"
BAMS[3]="LB_20241212_3_SI.bam"
BAMS[4]="LB_20241212_4_SI.bam"
BAMS[5]="LB_20241212_6_SI.bam"
BAMS[6]="LB_241114_1_SI.bam"
BAMS[7]="LB_241114_2_SI.bam"
BAMS[8]="LB_241114_3_SI.bam"
BAMS[9]="LB_241114_4_SI.bam"
BAMS[10]="LB_241114_5_SI.bam"
BAMS[11]="LB_241114_6_SI.bam"
BAMS[12]="LB_241114_7_SI.bam"
BAMS[13]="LB_241114_8_SI.bam"

# --- Get this task's sample ---
SAMPLE_ID=${SAMPLES[$SLURM_ARRAY_TASK_ID]}
BAM_FILE="$BAM_DIR/${BAMS[$SLURM_ARRAY_TASK_ID]}"
WIG_OUT="$WIG_DIR/${SAMPLE_ID}.wig"

echo "Sample ID:  $SAMPLE_ID"
echo "BAM file:   $BAM_FILE"
echo "WIG output: $WIG_OUT"
echo ""

# --- Check BAM exists ---
if [ ! -f "$BAM_FILE" ]; then
    echo "ERROR: BAM file not found: $BAM_FILE"
    exit 1
fi

if [ ! -f "${BAM_FILE}.bai" ]; then
    echo "ERROR: BAM index not found: ${BAM_FILE}.bai"
    exit 1
fi

# --- Skip if WIG already exists and is not empty ---
if [ -s "$WIG_OUT" ]; then
    echo "WIG file already exists and is non-empty — skipping"
    echo "Delete $WIG_OUT to force rerun"
    exit 0
fi

# --- Run readCounter ---
echo "Running readCounter..."
$READCOUNTER \
    --window $WINDOW \
    --quality $QUALITY \
    --chromosome "$CHROMOSOMES" \
    "$BAM_FILE" > "$WIG_OUT"

# --- Verify output ---
if [ -s "$WIG_OUT" ]; then
    LINES=$(wc -l < "$WIG_OUT")
    echo "SUCCESS: $SAMPLE_ID — WIG file created ($LINES lines)"
    echo "First 3 data lines:"
    head -4 "$WIG_OUT" | tail -3
else
    echo "ERROR: WIG file is empty for $SAMPLE_ID"
    rm -f "$WIG_OUT"
    exit 1
fi

echo ""
echo "Finished: $(date)"
