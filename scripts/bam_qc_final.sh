#!/bin/bash
# ============================================================
# PTCL ctDNA — Final BAM QC Script
# Author: [Noor Shaban]
# Date: May 2026
# Purpose: Comprehensive QC of all BAM files before analysis
# Run via: sbatch submit_bam_qc.sh
# ============================================================

# ── Directory paths ─────────────────────────────────────────
BAM_DIR="/scratch/alice/n/nhsas1/PTCL/BAMs"
OUT_DIR="/scratch/alice/n/nhsas1/PTCL/QC_output"
mkdir -p "$OUT_DIR"

# ── Output files ────────────────────────────────────────────
REPORT="${OUT_DIR}/bam_qc_report.txt"
TSV="${OUT_DIR}/bam_qc_summary.tsv"
COVERAGE="${OUT_DIR}/samtools_coverage_all.txt"
GENOME_CHECK="${OUT_DIR}/genome_check.txt"
DUP_REPORT="${OUT_DIR}/duplicate_check.txt"
INDEX_LOG="${OUT_DIR}/indexing_log.txt"

# ── Initialise output files ─────────────────────────────────
echo "=================================================" >  "$REPORT"
echo "PTCL ctDNA BAM QC Report"                         >> "$REPORT"
echo "Generated: $(date)"                               >> "$REPORT"
echo "User: ${USER}"                                    >> "$REPORT"
echo "BAM directory: ${BAM_DIR}"                        >> "$REPORT"
echo "Samtools version: $(samtools --version | head -1)" >> "$REPORT"
echo "=================================================" >> "$REPORT"
echo ""                                                 >> "$REPORT"

# TSV header
echo -e "sample\tfilename\ttotal_reads\tmapped_reads\tmapping_pct\t\
duplicate_reads\tduplicate_pct\tread_length\tcoverage_est\t\
reference_genome\thas_RG\tBQSR_confirmed\tdup_tool\tdup_status\t\
index_status\tbatch\tnotes" > "$TSV"

# ── Count BAM files ─────────────────────────────────────────
cd "$BAM_DIR"
BAM_COUNT=$(ls *.bam 2>/dev/null | wc -l)
echo "BAM files found: ${BAM_COUNT}" >> "$REPORT"
echo "" >> "$REPORT"

if [ "$BAM_COUNT" -eq 0 ]; then
    echo "ERROR: No BAM files found in ${BAM_DIR}" >> "$REPORT"
    exit 1
fi

# ============================================================
# SECTION 1 — BAM INDEXING
# Creates .bai index for any BAM missing one
# Required by IGV and all downstream tools
# ============================================================
echo "=================================================" >> "$REPORT"
echo "SECTION 1: BAM INDEXING"                          >> "$REPORT"
echo "=================================================" >> "$REPORT"
echo "Indexing started: $(date)"                        >> "$INDEX_LOG"

for BAM in *.bam; do
    BAI="${BAM}.bai"
    BAI_ALT="${BAM%.bam}.bai"

    if [ -f "$BAI" ] || [ -f "$BAI_ALT" ]; then
        echo "Already indexed: ${BAM}" >> "$REPORT"
        echo "SKIP  ${BAM}" >> "$INDEX_LOG"
    else
        echo "Indexing: ${BAM}" >> "$REPORT"
        samtools index "$BAM"
        if [ $? -eq 0 ]; then
            echo "OK    ${BAM}" >> "$INDEX_LOG"
        else
            echo "FAIL  ${BAM}" >> "$INDEX_LOG"
            echo "WARNING: Indexing failed for ${BAM}" >> "$REPORT"
        fi
    fi
done

INDEX_COUNT=$(ls *.bai 2>/dev/null | wc -l)
echo "Index files present after indexing: ${INDEX_COUNT}" >> "$REPORT"
echo "Indexing completed: $(date)" >> "$INDEX_LOG"
echo "" >> "$REPORT"

# ============================================================
# SECTION 2 — REFERENCE GENOME VERIFICATION
# Checks chr1 length from BAM header @SQ lines
# hg38: chr1 = 248,956,422 bp
# hg19: chr1 = 249,250,621 bp
# ============================================================
echo "=================================================" >> "$REPORT"
echo "SECTION 2: REFERENCE GENOME CHECK"               >> "$REPORT"
echo "=================================================" >> "$REPORT"
echo "Reference genome check — $(date)" > "$GENOME_CHECK"
echo "" >> "$GENOME_CHECK"

ALL_HG38=true

for BAM in *.bam; do
    # Extract chr1 length using tab-separated grep to avoid partial matches
    CHR1_LEN=$(samtools view -H "$BAM" \
        | awk '/^@SQ/ && /SN:chr1\t/' \
        | grep -o 'LN:[0-9]*' \
        | head -1 \
        | cut -d: -f2)

    if [ "$CHR1_LEN" = "248956422" ]; then
        REF="hg38"
        echo "hg38 CONFIRMED   ${BAM}" >> "$GENOME_CHECK"
        echo "hg38 CONFIRMED   ${BAM}" >> "$REPORT"
    elif [ "$CHR1_LEN" = "249250621" ]; then
        REF="hg19"
        ALL_HG38=false
        echo "hg19 DETECTED    ${BAM} — WARNING: inconsistent with cohort" >> "$GENOME_CHECK"
        echo "WARNING: hg19 detected in ${BAM}" >> "$REPORT"
    elif [ -z "$CHR1_LEN" ]; then
        REF="UNKNOWN"
        ALL_HG38=false
        echo "UNKNOWN          ${BAM} — chr1 @SQ line not found" >> "$GENOME_CHECK"
        echo "WARNING: Reference genome unknown for ${BAM}" >> "$REPORT"
    else
        REF="UNKNOWN_chr1=${CHR1_LEN}"
        ALL_HG38=false
        echo "UNKNOWN          ${BAM} — chr1 length=${CHR1_LEN}" >> "$GENOME_CHECK"
    fi
done

if [ "$ALL_HG38" = true ]; then
    echo "" >> "$REPORT"
    echo "RESULT: All BAM files confirmed hg38" >> "$REPORT"
    echo "" >> "$REPORT"
    echo "RESULT: All files confirmed hg38" >> "$GENOME_CHECK"
else
    echo "WARNING: Not all files confirmed hg38 — check genome_check.txt" >> "$REPORT"
fi
echo "" >> "$REPORT"

# ============================================================
# SECTION 3 — DUPLICATE HANDLING CHECK
# Reads Picard MarkDuplicates parameters from @PG header
# Reports REMOVE_DUPLICATES setting and duplicate read counts
# ============================================================
echo "=================================================" >> "$REPORT"
echo "SECTION 3: DUPLICATE HANDLING"                   >> "$REPORT"
echo "=================================================" >> "$REPORT"
echo "Duplicate check — $(date)" > "$DUP_REPORT"
echo "" >> "$DUP_REPORT"
printf "%-55s %8s %8s %8s %s\n" \
    "File" "Total" "Dups" "Dup%" "REMOVE_DUPS" >> "$DUP_REPORT"
echo "--------------------------------------------------------------------" >> "$DUP_REPORT"

for BAM in *.bam; do
    # Get duplicate handling setting from header
    REMOVE_DUPS=$(samtools view -H "$BAM" \
        | grep -i "MarkDuplicates" \
        | grep -o "REMOVE_DUPLICATES=[a-zA-Z]*" \
        | head -1)

    DUP_TOOL=$(samtools view -H "$BAM" \
        | grep "^@PG" \
        | grep -i "markdup\|MarkDuplicates\|sambamba\|picard" \
        | grep -o "PN:[^ \t]*" \
        | head -1 \
        | cut -d: -f2)

    # Get read counts from flagstat
    FLAGSTAT=$(samtools flagstat "$BAM" 2>/dev/null)
    TOTAL=$(echo "$FLAGSTAT" \
        | grep "in total" | awk '{print $1}')
    MAPPED=$(echo "$FLAGSTAT" \
        | grep " mapped (" | head -1 | awk '{print $1}')
    MAP_PCT=$(echo "$FLAGSTAT" \
        | grep " mapped (" | head -1 \
        | grep -o '([0-9.]*%' | tr -d '(')
    DUPS=$(echo "$FLAGSTAT" \
        | grep "duplicate" | head -1 | awk '{print $1}')

    # Calculate duplicate percentage
    if [ -n "$TOTAL" ] && [ "$TOTAL" -gt 0 ] && [ -n "$DUPS" ]; then
        DUP_PCT=$(echo "scale=2; $DUPS * 100 / $TOTAL" | bc)
    else
        DUP_PCT="N/A"
    fi

    # Get read length
    READ_LEN=$(samtools stats "$BAM" 2>/dev/null \
        | grep "^RL" \
        | awk '{print $2}' \
        | sort -n | tail -1)

    # Estimated coverage
    if [ -n "$TOTAL" ] && [ -n "$READ_LEN" ]; then
        COVERAGE_EST=$(echo "scale=2; \
            $TOTAL * $READ_LEN / 3200000000" | bc)
    else
        COVERAGE_EST="N/A"
    fi

    # Check for RG tag
    RG=$(samtools view -H "$BAM" \
        | grep "^@RG" | head -1)
    HAS_RG=$( [ -n "$RG" ] && echo "YES" || echo "NO" )

    # Check BQSR
    BQSR=$(samtools view -H "$BAM" \
        | grep -i "BQSR\|BaseRecalibrator\|ApplyBQSR" \
        | head -1)
    if echo "$BAM" | grep -qi "BQSR" || [ -n "$BQSR" ]; then
        BQSR_STATUS="YES"
    else
        BQSR_STATUS="NOT_CONFIRMED"
    fi

    # Determine batch from filename
    if echo "$BAM" | grep -q "LB20231211"; then
        BATCH="Batch1"
    else
        BATCH="Batch2"
    fi

    # Determine duplicate status label
    if [ "$REMOVE_DUPS" = "REMOVE_DUPLICATES=false" ]; then
        DUP_STATUS="Marked_only"
    elif [ "$REMOVE_DUPS" = "REMOVE_DUPLICATES=true" ]; then
        DUP_STATUS="Removed"
    else
        DUP_STATUS="Unknown"
    fi

    # Flag outlier
    NOTES=""
    if [ -n "$DUP_PCT" ] && [ "$DUP_PCT" != "N/A" ]; then
        OUTLIER=$(echo "$DUP_PCT > 14" | bc 2>/dev/null)
        if [ "$OUTLIER" = "1" ]; then
            NOTES="OUTLIER: high duplicate rate"
        fi
    fi

    # Check index status
    if [ -f "${BAM}.bai" ] || [ -f "${BAM%.bam}.bai" ]; then
        INDEX_STATUS="YES"
    else
        INDEX_STATUS="NO"
    fi

    # Write to duplicate report
    printf "%-55s %8s %8s %7s%% %s\n" \
        "$BAM" "$TOTAL" "$DUPS" "$DUP_PCT" "$REMOVE_DUPS" >> "$DUP_REPORT"

    # Write to main report
    echo "--- ${BAM} ---"                           >> "$REPORT"
    echo "  Batch:              ${BATCH}"           >> "$REPORT"
    echo "  Total reads:        ${TOTAL}"           >> "$REPORT"
    echo "  Mapped reads:       ${MAPPED} (${MAP_PCT})" >> "$REPORT"
    echo "  Duplicate reads:    ${DUPS} (${DUP_PCT}%)" >> "$REPORT"
    echo "  Duplicate tool:     ${DUP_TOOL:-Not found in header}" >> "$REPORT"
    echo "  Duplicate status:   ${DUP_STATUS}"      >> "$REPORT"
    echo "  Read length:        ${READ_LEN} bp"     >> "$REPORT"
    echo "  Coverage estimate:  ~${COVERAGE_EST}x"  >> "$REPORT"
    echo "  Reference genome:   ${REF}"             >> "$REPORT"
    echo "  RG tag:             ${HAS_RG}"          >> "$REPORT"
    echo "  BQSR:               ${BQSR_STATUS}"     >> "$REPORT"
    echo "  BAM index:          ${INDEX_STATUS}"    >> "$REPORT"
    [ -n "$NOTES" ] && echo "  NOTE: ${NOTES}"      >> "$REPORT"
    echo ""                                         >> "$REPORT"

    # Extract sample ID from filename for TSV
    SAMPLE_ID=$(echo "$BAM" \
        | grep -o "PRE_[0-9]*\|LB20231211_[0-9]*" \
        | head -1)

    # Write TSV row
    echo -e "${SAMPLE_ID}\t${BAM}\t${TOTAL}\t${MAPPED}\t\
${MAP_PCT}\t${DUPS}\t${DUP_PCT}%\t${READ_LEN}bp\t\
~${COVERAGE_EST}x\t${REF}\t${HAS_RG}\t${BQSR_STATUS}\t\
${DUP_TOOL}\t${DUP_STATUS}\t${INDEX_STATUS}\t${BATCH}\t${NOTES}" >> "$TSV"

done

echo "" >> "$REPORT"

# ============================================================
# SECTION 4 — SAMTOOLS COVERAGE
# Per-chromosome coverage statistics for all 21 samples
# Columns: chr, numreads, covbases, coverage%, meandepth,
#          meanbaseq, meanmapq
# ============================================================
echo "=================================================" >> "$REPORT"
echo "SECTION 4: SAMTOOLS COVERAGE (per chromosome)"   >> "$REPORT"
echo "=================================================" >> "$REPORT"
echo "Coverage analysis started: $(date)" >> "$REPORT"
echo "Full results in: ${COVERAGE}" >> "$REPORT"
echo "" >> "$REPORT"

echo "PTCL ctDNA — Samtools Coverage Report" > "$COVERAGE"
echo "Generated: $(date)" >> "$COVERAGE"
echo "" >> "$COVERAGE"

for BAM in *.bam; do
    echo "=================================================" >> "$COVERAGE"
    echo "FILE: ${BAM}"                                      >> "$COVERAGE"
    echo "=================================================" >> "$COVERAGE"

    # Run coverage — filter to standard chromosomes only
    samtools coverage "$BAM" \
        | awk 'NR==1 || $1 ~ /^chr([0-9]+|X|Y)$/' \
        >> "$COVERAGE"

    # Extract genome-wide mean depth for summary
    MEAN_DEPTH=$(samtools coverage "$BAM" \
        | awk '$1 ~ /^chr([0-9]+|X|Y)$/ {
            bases += $3 - $2
            depth += $7 * ($3 - $2)
          }
          END {
            if (bases > 0) printf "%.3f", depth/bases
          }')

    echo "" >> "$COVERAGE"
    echo "Genome-wide mean depth: ${MEAN_DEPTH}x" >> "$COVERAGE"
    echo "" >> "$COVERAGE"

    # Add summary line to main report
    echo "  ${BAM}: mean depth = ${MEAN_DEPTH}x" >> "$REPORT"
done

echo "" >> "$REPORT"
echo "Coverage analysis completed: $(date)" >> "$REPORT"

# ============================================================
# SECTION 5 — FINAL SUMMARY
# ============================================================
echo "=================================================" >> "$REPORT"
echo "FINAL SUMMARY"                                    >> "$REPORT"
echo "=================================================" >> "$REPORT"
echo "Total BAM files processed:  ${BAM_COUNT}"        >> "$REPORT"
echo "Index files present:        $(ls *.bai 2>/dev/null | wc -l)" >> "$REPORT"
echo "All hg38 confirmed:         ${ALL_HG38}"         >> "$REPORT"
echo "" >> "$REPORT"
echo "Output files produced:"                          >> "$REPORT"
echo "  ${REPORT}"                                     >> "$REPORT"
echo "  ${TSV}"                                        >> "$REPORT"
echo "  ${COVERAGE}"                                   >> "$REPORT"
echo "  ${GENOME_CHECK}"                               >> "$REPORT"
echo "  ${DUP_REPORT}"                                 >> "$REPORT"
echo "  ${INDEX_LOG}"                                  >> "$REPORT"
echo "" >> "$REPORT"
echo "QC completed: $(date)"                           >> "$REPORT"
echo "=================================================" >> "$REPORT"

echo "All QC complete. Results in: ${OUT_DIR}"
