#!/usr/bin/env python3
# ============================================================
# PTCL ctDNA — Add original_filename column to ichorCNA_summary.csv
# Author: Noor Shaban
# Date: June 2026
# Purpose: Add original BAM filename column to ichorCNA summary
# ============================================================

import csv
import os

INPUT  = '/scratch/alice/n/nhsas1/ptcl_ctdna_thesis/results/ichorCNA/ichorCNA_summary.csv'
OUTPUT = '/scratch/alice/n/nhsas1/ptcl_ctdna_thesis/results/ichorCNA/ichorCNA_summary.csv'
BACKUP = '/scratch/alice/n/nhsas1/ptcl_ctdna_thesis/results/ichorCNA/ichorCNA_summary_backup.csv'

# ── Complete sample ID to original BAM filename mapping ─────
FILENAME_MAP = {
    # Batch 1
    'B1_S01': 'LB20231211_1_EKDL240006162-1A_227NJFLT4_L7_dD-RG-BQSR_FINAL.bam',
    'B1_S02': 'LB20231211_2_EKDL240006163-1A_227NJFLT4_L7_dD-RG-BQSR_FINAL.bam',
    'B1_S03': 'LB20231211_3_EKDL240006164-1A_227NJFLT4_L7_dD-RG-BQSR_FINAL.bam',
    'B1_S04': 'LB20231211_4_EKDL240006165-1A_227NJFLT4_L7_dD-RG-BQSR_FINAL.bam',
    'B1_S05': 'LB20231211_6_EKDL240006167-1A_227NJFLT4_L7_dD-RG-BQSR_FINAL.bam',
    # Batch 2
    'B2_S01': 'LB_20251209_PRE_1_EKDL250035517-1A_23_merge_dD_SI.bam',
    'B2_S02': 'LB_20251209_PRE_2_EKDL250035518-1A_237N3GLT4_L2_dD_SI.bam',
    'B2_S03': 'LB_20251209_PRE_3_EKDL250035519-1A_237N3GLT4_L2_dD_SI.bam',
    'B2_S04': 'LB_20251209_PRE_4_EKDL250035520-1A_237N3GLT4_L2_dD_SI.bam',
    'B2_S05': 'LB_20251209_PRE_5_EKDL250035521-1A_237N3GLT4_L1_dD_SI.bam',
    'B2_S06': 'LB_20251209_PRE_6_EKDL250035522-1A_23_merge_dD_SI.bam',
    'B2_S07': 'LB_20251209_PRE_7_EKDL250035523-1A_23_merge_dD_SI.bam',
    'B2_S08': 'LB_20251209_PRE_8_EKDL250035524-1A_23_merge_dD_SI.bam',
    'B2_S09': 'LB_20251210_PRE_1_EKDL250035525-1A_23_merge_dD_SI.bam',
    'B2_S10': 'LB_20251210_PRE_2_EKDL250035526-1A_237N3GLT4_L8_dD_SI.bam',
    'B2_S11': 'LB_20251210_PRE_3_EKDL250035527-1A_237N3GLT4_L8_dD_SI.bam',
    'B2_S12': 'LB_20251210_PRE_4_EKDL250035528-1A_237N3GLT4_L8_dD_SI.bam',
    'B2_S13': 'LB_20251210_PRE_5_EKDL250035529-1A_23_merge_dD_SI.bam',
    'B2_S14': 'LB_20251210_PRE_6_EKDL250035530-1A_237N3GLT4_L5_dD_SI.bam',
    'B2_S15': 'LB_20251210_PRE_7_EKDL250035531-1A_237N3GLT4_L5_dD_SI.bam',
    'B2_S16': 'LB_20251210_PRE_8_EKDL250035532-1A_237N3GLT4_L5_dD_SI.bam',
    # Batch 3
    'B3_S01': 'LB_20241212_1_SI.bam',
    'B3_S02': 'LB_20241212_2_SI.bam',
    'B3_S03': 'LB_20241212_3_SI.bam',
    'B3_S04': 'LB_20241212_4_SI.bam',
    'B3_S05': 'LB_20241212_6_SI.bam',
    'B3_S06': 'LB_241114_1_SI.bam',
    'B3_S07': 'LB_241114_2_SI.bam',
    'B3_S08': 'LB_241114_3_SI.bam',
    'B3_S09': 'LB_241114_4_SI.bam',
    'B3_S10': 'LB_241114_5_SI.bam',
    'B3_S11': 'LB_241114_6_SI.bam',
    'B3_S12': 'LB_241114_7_SI.bam',
    'B3_S13': 'LB_241114_8_SI.bam',
}

# ── Read existing CSV ────────────────────────────────────────
with open(INPUT, 'r') as f:
    reader = csv.DictReader(f)
    rows = list(reader)
    original_fieldnames = reader.fieldnames

print(f"Read {len(rows)} rows from {INPUT}")
print(f"Existing columns: {original_fieldnames}")

# ── Backup original ──────────────────────────────────────────
import shutil
shutil.copy(INPUT, BACKUP)
print(f"Backup saved to {BACKUP}")

# ── Add original_filename column ─────────────────────────────
new_fieldnames = ['sample', 'original_filename', 'tumor_fraction', 'ploidy', 'MAD']

updated_rows = []
missing = []
for row in rows:
    sample = row['sample']
    orig = FILENAME_MAP.get(sample, 'UNKNOWN')
    if orig == 'UNKNOWN':
        missing.append(sample)
    row['original_filename'] = orig
    updated_rows.append(row)

if missing:
    print(f"WARNING: No mapping found for: {missing}")
else:
    print("All samples mapped successfully")

# ── Write updated CSV ────────────────────────────────────────
with open(OUTPUT, 'w', newline='') as f:
    writer = csv.DictWriter(f, fieldnames=new_fieldnames)
    writer.writeheader()
    writer.writerows(updated_rows)

print(f"\nUpdated CSV saved to {OUTPUT}")
print(f"Columns: {new_fieldnames}")
print(f"\nFinal content:")
with open(OUTPUT, 'r') as f:
    print(f.read())
