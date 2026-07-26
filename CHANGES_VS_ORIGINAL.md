# What differs between this review repo and the original

Compares `ptcl_ctdna_thesis_review` (branch `PTCL`) against `ptcl_ctdna_thesis` at
`4e59a75`, the state before the review began.

The original repository is **untouched**. Its `main` is still at `4e59a75` and no branch was
left behind on it.

## At a glance

| | |
|---|---|
| Commits | 63 |
| Files changed | 86 (11 new, 0 deleted) |
| Lines | +2,547 / −171 |
| Comment lines added to scripts | ~1,240 |
| Scripts touched | 74 of 77 |
| Bug fixes | 21 |
| Behaviour-preserving refactors | 2 |
| New reported values | 2 |

Commit types: 35 `docs`, 21 `fix`, 2 `refactor`, 2 `feat`, 2 `chore`, 1 `perf`.

## New files

| File | Purpose |
|---|---|
| `README.md` (rewritten) | Was a stub with an empty structure section. Now carries cohort derivation, execution order, environment, limitations. |
| `RERUN_REQUIRED.md` | Which outputs predate a fix, which were verified unaffected, and open data questions. |
| `CLEANUP_PROPOSAL.md` | Structural proposals, with applied items marked. |
| `CHANGES_VS_ORIGINAL.md` | This file. |
| `CLAUDE.md` | Was untracked in the original; now versioned. |
| `scripts/archive/` + README | Superseded concordance v1 and its wrapper, quarantined. |
| `results/ichorCNA/params/diagnostic_runs/` | 7 params files recovered from ALICE, closing the B3_S12 evidence gap. |

## The 21 fixes, by consequence

### Changed a committed number

| Fix | Effect |
|---|---|
| `fix(ichorcna): Call arm direction by fraction` | 7 of 780 arm rows flip GAIN→LOSS across 6 samples. AS totals unchanged — `scored_25` identical for every row. |
| `fix(rascal): Compare implied TF, not cellularity` | `tf_difference` wrong for most of the batch. Resolves a contradiction where two committed CSVs disagreed on the same fits (B3_S06: 0.337 vs 0.318). No agreement label changes. |
| `fix(fragmentomics): Recompute z after exclusion` | **Largest downstream effect.** Every per-bin p-value changes. See "Results" below. |
| `fix(fragmentomics): Correct omnibus p-values` | Kruskal-Wallis omnibus was never BH-corrected despite the Methods saying so. Primary analysis survives (0.037–0.038); the B2_S11-excluded sensitivity analysis does not (0.055). |
| `fix(fragmentomics): Add BH to correlations` | 8 Spearman tests were uncorrected while step 3 corrected. All three significant results survive (0.004–0.010). |
| `fix(fragmentomics): Use agreement statistics` | Pearson *r* replaced by Bland-Altman on the dedup figures. `delta_sl_ratio` is positive in 13/13 samples (sign test p=0.0002) while *r* read 0.9998. |
| `fix(ichorcna): Stop reporting zero p-value` | One cell: `proportional_bias_pvalue` was literally `0` from underflow. |
| `fix(rascal): Label discordant fits accurately` | One cell: B3_S09 read "IMPROVED" while 0.199 from ichorCNA against a 0.10 band. |

### Prevented a silent failure

| Fix | What it prevented |
|---|---|
| `fix(gistic): Fail on missing seg file` | A missing seg file silently shrank the cohort, and GISTIC q-values are computed against sample count. |
| `fix(fragmentomics): Detect failed chunks` | The submit driver treated FAILED/TIMEOUT as success and carried on. |
| `fix(fragmentomics): Warn on truncated histogram` | No completeness check on fragment histograms. |
| `fix(qdnaseq): Report seg files that fail to parse` ×2 | Malformed seg files dropped from the master table but still counted as completed. |
| `fix(rascal): Record grid-minimum fallback` | A failed fit was indistinguishable from a real one — `n_solutions=1`, `ambiguous=FALSE`, `status=SUCCESS`. |
| `fix(qdnaseq): Create sample dir before plotting` | Diagnostic plots vanished for exactly the samples that had failed. |
| `fix(fragmentomics): Unhardcode step 5f bin` | Hardcoded chr4 lookup would report the wrong bin after re-running. |

### Correctness and robustness

`fix(fragmentomics): Order chromosomes correctly` — bin index built from an unfactored
character column sorted chr1, chr10, chr11 … chr2. Fixed at source so step 5c no longer
needs to overwrite step 5's figures.

`fix(qdnaseq): Stop clipping raw read count plot` · `fix(qdnaseq): Repair installer on clean
account` (covered 3 of 13 required packages) · `perf(fragmentomics): Speed up GC counting`
(~3×10⁹ awk iterations under a 1-hour walltime) · `fix(fragmentomics): Write full exclusion
record` · `fix(fragmentomics): Repair step 5d print` (a bug introduced by the preceding fix
and caught by re-running).

### Refactors — behaviour preserved

`refactor(ichorcna): Archive concordance v1` and `refactor(fragmentomics): Make 5c
diagnostic-only`. Both removed cases where a superseded script silently overwrote current
results with an older analysis.

## Results status

**Complete.** All eight ALICE jobs completed, the outputs were synced and committed
(`ab72653`), and this repository now carries the regenerated results. Verified:

- `genomewide_perbin_significance.csv` shows chr4 and chr19 both at `p_adj = 0.0499`
- `tf_correlation_results.csv` carries `p_adjusted`
- `RASCAL_batch3_summary.csv` carries `rascal_expected_tf` and `fallback_used`

Not regenerated, and correctly so: `fragmentomics_metrics_summary.csv` and
`batch3_dedup_vs_original_comparison.csv` are step 2 outputs and step 2 was not re-run;
`group_comparison_summary_stats.csv` holds medians and quartiles, which the correction does
not touch. `fig_batch3_longfragfraction_dedup_sensitivity.png` also did not regenerate while
both siblings did — confirming it has no source script.

### What the re-run showed

| | before | after |
|---|---|---|
| chr4:5–10Mb | p_adj = 0.0280 | p_adj = 0.0499 |
| chr19:15–20Mb | not significant | p_adj = 0.0499 |
| chr19 in top 20 bins | 3 | 4 |

Both land at 0.0499 because Benjamini-Hochberg takes a running minimum: chr4's own value is
0.0579, pulled down by chr19. **Remove chr19 and chr4 is not significant.** The chr19
enrichment also strengthened (hypergeometric p 0.0046 → 0.00026), so the z-score bug was not
its cause — which leaves the GC correction, since chr19 is the most GC-rich chromosome and
this pipeline corrects the ratio rather than the counts.

## Deliberately not changed

Flagged in place rather than fixed, because changing them would be new analysis:

- **GC correction applied to the ratio, not the counts** — the decisive open question above.
- **`expected_logR` exact only at ploidy 2** — drives a hand-reviewed triage step; 20 of 34
  samples sit >0.1 from ploidy 2. Correcting it would retrospectively alter the basis of
  curation already performed.
- **RASCAL `min_ploidy = 1.5` binds** — 5 of 13 samples fit at exactly 1.50.
- **RASCAL ±0.15 constraint applied to cellularity, not tumour fraction.**
- **`step5f` measures terminal-bin position, not telomeric content** — 39 flagged bins are
  17 p-arm + 22 q-arm termini.
- **`long_fragment_fraction` never measured contamination** — the `-f 2` proper-pair filter
  caps observable length at 242–294bp per library (B2_S11: 864bp), removing the
  di-nucleosomal population. The six zero values are arithmetically correct.

## Open questions for supervision

1. **GC correction** — decides whether step 5e has any finding at all.
2. **B2_S15** — shows the same ploidy non-identifiability as B3_S12 (1.854 vs 2.598 by
   starting point) but was retained and curated into the most-altered genome in the cohort.
3. **`long_fragment_fraction`** — drop it, or re-extract without `-f 2`.
4. **FGA vs GISTIC** — ρ=0.450, p=0.0467 raw, 0.126 after BH.
5. **TF vs FGA is negative** (−0.295) — noise-driven FGA at the low-tumour-fraction end, or
   real biology?
6. **Coverage range** — documentation says 0.87–1.77x; `bam_qc_summary.tsv` gives 0.93 as
   the minimum.

## Not reviewed

`QC_results/` and the six root-level QC scripts were left out of scope. Two things were
observed before stopping: `sample_technical_metadata.tsv` contains only a failed heredoc
opener and holds no data, and 29 of 34 samples record `has_RG = NO` /
`BQSR_confirmed = NOT_CONFIRMED`.
