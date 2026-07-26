#!/bin/bash
# Step 1 for the deduplicated Batch 3 BAMs, feeding the dedup sensitivity analysis.
#
# Uses -F 3328 rather than 2304 because Picard has now marked duplicates in these files, so
# the 1024 bit is present and must be filtered. See run_fragment_extraction.slurm.
#
# The 13 sample IDs are hardcoded here rather than read from fragment_sample_sheet.tsv,
# so this is a second place the Batch 3 sample list must be kept in step.
#SBATCH --job-name=frag_extract_b3dedup
#SBATCH --partition=short
#SBATCH --array=1-13%4
#SBATCH --cpus-per-task=4
#SBATCH --mem=4G
#SBATCH --time=02:00:00
#SBATCH --output=/scratch/alice/n/nhsas1/PTCL/fragmentomics/logs/frag_extract_b3dedup_%A_%a.out
#SBATCH --error=/scratch/alice/n/nhsas1/PTCL/fragmentomics/logs/frag_extract_b3dedup_%A_%a.err

module load samtools/1.18

SAMPLES=(B3_S01 B3_S02 B3_S03 B3_S04 B3_S05 B3_S06 B3_S07 B3_S08 B3_S09 B3_S10 B3_S11 B3_S12 B3_S13)
SAMPLE=${SAMPLES[$((SLURM_ARRAY_TASK_ID - 1))]}
BAMPATH=/scratch/alice/n/nhsas1/PTCL/fragmentomics/BAMs_batch3_dedup/${SAMPLE}_dedup.bam

OUTDIR=/scratch/alice/n/nhsas1/PTCL/fragmentomics/fragment_histograms_batch3_dedup
mkdir -p ${OUTDIR}
mkdir -p /scratch/alice/n/nhsas1/PTCL/fragmentomics/logs

if [ ! -f "${BAMPATH}" ]; then
    echo "ERROR: BAM not found for ${SAMPLE}: ${BAMPATH}"
    exit 1
fi

echo "Processing ${SAMPLE} (Batch3, deduplicated), BAM: ${BAMPATH}, exclude flag 3328"

OUTFILE=${OUTDIR}/${SAMPLE}_dedup_fragment_length_hist.tsv

samtools view -@ ${SLURM_CPUS_PER_TASK} -q 20 -f 2 -F 3328 "${BAMPATH}" \
    chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 | \
    awk '{if ($9 >= 50 && $9 <= 1000) print $9}' | \
    sort -n | \
    uniq -c | \
    awk 'BEGIN{OFS="\t"; print "fragment_length", "count"} {print $2, $1}' \
    > ${OUTFILE}

N_FRAGMENTS=$(awk 'NR>1{sum+=$2} END{print sum+0}' ${OUTFILE})
echo "${SAMPLE}: ${N_FRAGMENTS} fragments extracted, written to ${OUTFILE}"
