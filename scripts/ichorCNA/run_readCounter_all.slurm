#!/bin/bash
#SBATCH --job-name=readCounter_all
#SBATCH --partition=short
#SBATCH --time=01:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --array=11-20%4
#SBATCH --output=/scratch/alice/n/nhsas1/PTCL/ichorCNA/logs/readCounter_%a_%j.out
#SBATCH --error=/scratch/alice/n/nhsas1/PTCL/ichorCNA/logs/readCounter_%a_%j.err

# readCounter — All remaining 20 samples (SLURM array job)

echo "Array task:  $SLURM_ARRAY_TASK_ID"
echo "Job started: $(date)"
echo "Node:        $HOSTNAME"
echo ""

# Tool and shared paths
READCOUNTER=/scratch/alice/n/nhsas1/PTCL/hmmcopy_utils/bin/readCounter
BAM_DIR=/scratch/alice/n/nhsas1/PTCL/BAMs
WIG_DIR=/scratch/alice/n/nhsas1/PTCL/ichorCNA/wig
CHROMOSOMES="chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22"
WINDOW=1000000
QUALITY=20

# Sample list (array index 1-20 maps to samples B1_S02 to B2_S16)
# B1_S01 is already done so it is  excluded from this array
declare -A SAMPLES
SAMPLES[1]="B1_S02"
SAMPLES[2]="B1_S03"
SAMPLES[3]="B1_S04"
SAMPLES[4]="B1_S05"
SAMPLES[5]="B2_S01"
SAMPLES[6]="B2_S02"
SAMPLES[7]="B2_S03"
SAMPLES[8]="B2_S04"
SAMPLES[9]="B2_S05"
SAMPLES[10]="B2_S06"
SAMPLES[11]="B2_S07"
SAMPLES[12]="B2_S08"
SAMPLES[13]="B2_S09"
SAMPLES[14]="B2_S10"
SAMPLES[15]="B2_S11"
SAMPLES[16]="B2_S12"
SAMPLES[17]="B2_S13"
SAMPLES[18]="B2_S14"
SAMPLES[19]="B2_S15"
SAMPLES[20]="B2_S16"

declare -A BAMS
BAMS[1]="LB20231211_2_EKDL240006163-1A_227NJFLT4_L7_dD-RG-BQSR_FINAL.bam"
BAMS[2]="LB20231211_3_EKDL240006164-1A_227NJFLT4_L7_dD-RG-BQSR_FINAL.bam"
BAMS[3]="LB20231211_4_EKDL240006165-1A_227NJFLT4_L7_dD-RG-BQSR_FINAL.bam"
BAMS[4]="LB20231211_6_EKDL240006167-1A_227NJFLT4_L7_dD-RG-BQSR_FINAL.bam"
BAMS[5]="LB_20251209_PRE_1_EKDL250035517-1A_23_merge_dD_SI.bam"
BAMS[6]="LB_20251209_PRE_2_EKDL250035518-1A_237N3GLT4_L2_dD_SI.bam"
BAMS[7]="LB_20251209_PRE_3_EKDL250035519-1A_237N3GLT4_L2_dD_SI.bam"
BAMS[8]="LB_20251209_PRE_4_EKDL250035520-1A_237N3GLT4_L2_dD_SI.bam"
BAMS[9]="LB_20251209_PRE_5_EKDL250035521-1A_237N3GLT4_L1_dD_SI.bam"
BAMS[10]="LB_20251209_PRE_6_EKDL250035522-1A_23_merge_dD_SI.bam"
BAMS[11]="LB_20251209_PRE_7_EKDL250035523-1A_23_merge_dD_SI.bam"
BAMS[12]="LB_20251209_PRE_8_EKDL250035524-1A_23_merge_dD_SI.bam"
BAMS[13]="LB_20251210_PRE_1_EKDL250035525-1A_23_merge_dD_SI.bam"
BAMS[14]="LB_20251210_PRE_2_EKDL250035526-1A_237N3GLT4_L8_dD_SI.bam"
BAMS[15]="LB_20251210_PRE_3_EKDL250035527-1A_237N3GLT4_L8_dD_SI.bam"
BAMS[16]="LB_20251210_PRE_4_EKDL250035528-1A_237N3GLT4_L8_dD_SI.bam"
BAMS[17]="LB_20251210_PRE_5_EKDL250035529-1A_23_merge_dD_SI.bam"
BAMS[18]="LB_20251210_PRE_6_EKDL250035530-1A_237N3GLT4_L5_dD_SI.bam"
BAMS[19]="LB_20251210_PRE_7_EKDL250035531-1A_237N3GLT4_L5_dD_SI.bam"
BAMS[20]="LB_20251210_PRE_8_EKDL250035532-1A_237N3GLT4_L5_dD_SI.bam"

# Get this task's sample
SAMPLE_ID=${SAMPLES[$SLURM_ARRAY_TASK_ID]}
BAM_FILE="$BAM_DIR/${BAMS[$SLURM_ARRAY_TASK_ID]}"
WIG_OUT="$WIG_DIR/${SAMPLE_ID}.wig"

echo "Sample ID:  $SAMPLE_ID"
echo "BAM file:   $BAM_FILE"
echo "WIG output: $WIG_OUT"
echo ""

# Check BAM exists 
if [ ! -f "$BAM_FILE" ]; then
    echo "ERROR: BAM file not found: $BAM_FILE"
    exit 1
fi

if [ ! -f "${BAM_FILE}.bai" ]; then
    echo "ERROR: BAM index not found: ${BAM_FILE}.bai"
    exit 1
fi

#  Skip if WIG already exists and is not empty 
if [ -s "$WIG_OUT" ]; then
    echo "WIG file already exists and is non-empty — skipping"
    echo "Delete $WIG_OUT to force rerun"
    exit 0
fi

#  Run readCounter 
echo "Running readCounter..."
$READCOUNTER \
    --window $WINDOW \
    --quality $QUALITY \
    --chromosome "$CHROMOSOMES" \
    "$BAM_FILE" > "$WIG_OUT"

#  Verify output 
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
