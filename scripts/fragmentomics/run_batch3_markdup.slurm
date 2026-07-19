#!/bin/bash
#SBATCH --job-name=b3_markdup
#SBATCH --partition=short
#SBATCH --array=22-34%4
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=/scratch/alice/n/nhsas1/PTCL/fragmentomics/logs/markdup_%A_%a.out
#SBATCH --error=/scratch/alice/n/nhsas1/PTCL/fragmentomics/logs/markdup_%A_%a.err

module load gcc/12.3.0-yxgv2bl
module load picard/3.0.0-56ewicy
module load samtools/1.18

SHEET=/scratch/alice/n/nhsas1/ptcl_ctdna_thesis/scripts/fragmentomics/fragment_sample_sheet.tsv
OUTDIR=/scratch/alice/n/nhsas1/PTCL/fragmentomics/BAMs_batch3_dedup
METRICSDIR=/scratch/alice/n/nhsas1/PTCL/fragmentomics/markdup_metrics

mkdir -p ${OUTDIR}
mkdir -p ${METRICSDIR}
mkdir -p /scratch/alice/n/nhsas1/PTCL/fragmentomics/logs

LINE=$(awk -v n=${SLURM_ARRAY_TASK_ID} 'NR==n+1' ${SHEET})
SAMPLE=$(echo "${LINE}" | cut -f2)
BATCH=$(echo "${LINE}" | cut -f3)
BAMPATH=$(echo "${LINE}" | cut -f4)

if [ "${BATCH}" != "Batch3" ]; then
    echo "ERROR: array index ${SLURM_ARRAY_TASK_ID} is not a Batch3 sample (got ${SAMPLE}, ${BATCH}). Aborting."
    exit 1
fi

if [ ! -f "${BAMPATH}" ]; then
    echo "ERROR: BAM not found for ${SAMPLE}: ${BAMPATH}"
    exit 1
fi

echo "Marking duplicates for ${SAMPLE} (${BATCH})"

OUTBAM=${OUTDIR}/${SAMPLE}_dedup.bam
METRICS=${METRICSDIR}/${SAMPLE}_markdup_metrics.txt

export JAVA_TOOL_OPTIONS="-Xmx6g"
picard MarkDuplicates \
    --INPUT "${BAMPATH}" \
    --OUTPUT "${OUTBAM}" \
    --METRICS_FILE "${METRICS}" \
    --REMOVE_DUPLICATES false \
    --READ_NAME_REGEX null

if [ $? -ne 0 ]; then
    echo "ERROR: MarkDuplicates failed for ${SAMPLE}"
    exit 1
fi

samtools index -@ ${SLURM_CPUS_PER_TASK} "${OUTBAM}"

DUP_PCT=$(grep -A1 "^LIBRARY" "${METRICS}" | tail -1 | awk -F'\t' '{print $9}')
echo "${SAMPLE}: MarkDuplicates complete, duplicate fraction = ${DUP_PCT}"
