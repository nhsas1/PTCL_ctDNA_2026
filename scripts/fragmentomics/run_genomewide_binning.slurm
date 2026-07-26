#!/bin/bash
# Step 5 input: bins fragments into 5Mb windows genome-wide, counting short (100-150bp) and
# long (151-220bp) fragments per bin. Cut-offs match calculate_fragment_metrics.R exactly.
#
# Fragments are assigned by midpoint, so each is counted once. Note SAM POS is 1-based
# while compute_gc_content_bins.sh emits 0-based bin starts, so a midpoint landing exactly
# on a 5Mb boundary is assigned one bin high. That is one base in five million and does not
# affect any result, but the two coordinate conventions do differ.
#
# No --array directive is set here; the range comes from submit_chunked.sh. Launching this
# with a bare sbatch produces an empty task ID and the guard below catches it.
#SBATCH --job-name=frag_genomebin
#SBATCH --partition=short
#SBATCH --cpus-per-task=4
#SBATCH --mem=4G
#SBATCH --time=02:00:00
#SBATCH --output=/scratch/alice/n/nhsas1/PTCL/fragmentomics/logs/frag_genomebin_%A_%a.out
#SBATCH --error=/scratch/alice/n/nhsas1/PTCL/fragmentomics/logs/frag_genomebin_%A_%a.err

module load samtools/1.18

SHEET=/scratch/alice/n/nhsas1/ptcl_ctdna_thesis/scripts/fragmentomics/fragment_sample_sheet.tsv
OUTDIR=/scratch/alice/n/nhsas1/PTCL/fragmentomics/genomewide_bins
BATCH3_DEDUP_DIR=/scratch/alice/n/nhsas1/PTCL/fragmentomics/BAMs_batch3_dedup

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

if [ "${BATCH}" == "Batch3" ]; then
    BAMPATH=${BATCH3_DEDUP_DIR}/${SAMPLE}_dedup.bam
    EXCLUDE_FLAG=3328
else
    EXCLUDE_FLAG=3328
fi

if [ ! -f "${BAMPATH}" ]; then
    echo "ERROR: BAM not found for ${SAMPLE}: ${BAMPATH}"
    exit 1
fi

echo "Processing ${SAMPLE} (${BATCH}), BAM: ${BAMPATH}"

OUTFILE=${OUTDIR}/${SAMPLE}_5Mb_bins.tsv

samtools view -@ ${SLURM_CPUS_PER_TASK} -q 20 -f 2 -F ${EXCLUDE_FLAG} "${BAMPATH}" \
    chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 | \
    awk -v binsize=5000000 '
    {
        tlen = $9
        if (tlen <= 0 || tlen < 50 || tlen > 1000) next
        chr = $3
        pos = $4
        mid = pos + int(tlen/2)
        bin = int(mid / binsize) * binsize
        key = chr "\t" bin
        total[key]++
        if (tlen >= 100 && tlen <= 150) short[key]++
        else if (tlen >= 151 && tlen <= 220) long[key]++
    }
    END {
        print "chr\tbin_start\tshort_count\tlong_count\ttotal_count"
        for (k in total) {
            print k "\t" (short[k]+0) "\t" (long[k]+0) "\t" total[k]
        }
    }' | sort -k1,1 -k2,2n > ${OUTFILE}

N_BINS=$(tail -n +2 ${OUTFILE} | wc -l)
echo "${SAMPLE}: ${N_BINS} genomic bins written to ${OUTFILE}"
