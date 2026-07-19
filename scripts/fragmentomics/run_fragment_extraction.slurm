#!/bin/bash
#SBATCH --job-name=frag_extract
#SBATCH --partition=short
#SBATCH --array=1-34
#SBATCH --cpus-per-task=4
#SBATCH --mem=4G
#SBATCH --time=02:00:00
#SBATCH --output=/scratch/alice/n/nhsas1/PTCL/fragmentomics/logs/frag_extract_%A_%a.out
#SBATCH --error=/scratch/alice/n/nhsas1/PTCL/fragmentomics/logs/frag_extract_%A_%a.err

module load samtools/1.18

SHEET=/scratch/alice/n/nhsas1/ptcl_ctdna_thesis/scripts/fragmentomics/fragment_sample_sheet.tsv
OUTDIR=/scratch/alice/n/nhsas1/PTCL/fragmentomics/fragment_histograms

mkdir -p ${OUTDIR}
mkdir -p /scratch/alice/n/nhsas1/PTCL/fragmentomics/logs

LINE=$(awk -v n=${SLURM_ARRAY_TASK_ID} 'NR==n+1' ${SHEET})
SAMPLE=$(echo "${LINE}" | cut -f2)
BATCH=$(echo "${LINE}" | cut -f3)
BAMPATH=$(echo "${LINE}" | cut -f4)

if [ -z "${SAMPLE}" ]; then
    echo "ERROR: no sample sheet entry for array task ${SLURM_ARRAY_TASK_ID}"
    exit 1
fi

if [ ! -f "${BAMPATH}" ]; then
    echo "ERROR: BAM not found for ${SAMPLE}: ${BAMPATH}"
    exit 1
fi

if [ "${BATCH}" == "Batch3" ]; then
    EXCLUDE_FLAG=2304
else
    EXCLUDE_FLAG=3328
fi

echo "Processing ${SAMPLE} (${BATCH}), BAM: ${BAMPATH}, exclude flag ${EXCLUDE_FLAG}"

OUTFILE=${OUTDIR}/${SAMPLE}_fragment_length_hist.tsv

samtools view -@ ${SLURM_CPUS_PER_TASK} -q 20 -f 2 -F ${EXCLUDE_FLAG} "${BAMPATH}" \
    chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 | \
    awk '{if ($9 >= 50 && $9 <= 1000) print $9}' | \
    sort -n | \
    uniq -c | \
    awk 'BEGIN{OFS="\t"; print "fragment_length", "count"} {print $2, $1}' \
    > ${OUTFILE}

N_FRAGMENTS=$(awk 'NR>1{sum+=$2} END{print sum+0}' ${OUTFILE})
echo "${SAMPLE}: ${N_FRAGMENTS} fragments extracted, written to ${OUTFILE}"
