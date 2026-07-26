#!/bin/bash
#SBATCH --job-name=frag_gccontent
#SBATCH --partition=short
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH --output=/scratch/alice/n/nhsas1/PTCL/fragmentomics/logs/frag_gccontent_%j.out
#SBATCH --error=/scratch/alice/n/nhsas1/PTCL/fragmentomics/logs/frag_gccontent_%j.err

module load samtools/1.18

REF=/scratch/alice/n/nhsas1/PTCL/reference/hg38.fa
OUTDIR=/scratch/alice/n/nhsas1/PTCL/fragmentomics/genomewide_bins
BINSIZE=5000000

mkdir -p ${OUTDIR}
OUTFILE=${OUTDIR}/gc_content_5mb_bins.tsv

echo -e "chr\tbin_start\tgc_fraction\tvalid_base_fraction" > ${OUTFILE}

for chr in chr{1..22}; do
    CHRLEN=$(awk -v c=${chr} '$1==c {print $2}' ${REF}.fai)
    echo "Processing ${chr} (length ${CHRLEN})"

    start=0
    while [ ${start} -lt ${CHRLEN} ]; do
        end=$((start + BINSIZE))
        if [ ${end} -gt ${CHRLEN} ]; then
            end=${CHRLEN}
        fi
        region_start=$((start + 1))

        seq=$(samtools faidx ${REF} ${chr}:${region_start}-${end} | tail -n +2 | tr -d '\n')

        echo "${seq}" | awk -v chr=${chr} -v bstart=${start} '
        {
            seq = toupper($0)
            n = length(seq)
            # Count bases by deleting everything else and measuring the remaining length.
            # This replaces a per-character substr() loop which, across the ~2.9Gb genome,
            # ran roughly 3e9 iterations and was very likely to exceed the one hour
            # walltime. A timeout mid-run would have been quietly damaging: the last
            # chromosomes would simply be missing from the GC table, and the merge in
            # step5_genomewide_profile.R would silently drop their bins rather than error.
            t = seq; gsub(/[^GC]/, "", t); gc = length(t)
            t = seq; gsub(/[^AT]/, "", t); at = length(t)
            nn = n - gc - at
            valid = gc + at
            total = valid + nn
            gc_frac = (valid > 0) ? gc / valid : "NA"
            valid_frac = (total > 0) ? valid / total : 0
            print chr "\t" bstart "\t" gc_frac "\t" valid_frac
        }' >> ${OUTFILE}

        start=${end}
    done
done

N_BINS=$(tail -n +2 ${OUTFILE} | wc -l)
echo "Done. ${N_BINS} bins written to ${OUTFILE}"
