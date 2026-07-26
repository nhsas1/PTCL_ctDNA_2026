# Outputs requiring regeneration

Bug fixes applied during the code review changed analysis logic. The files in `results/`
were produced before those fixes and are listed here with what changed and what to re-run.

Raw data lives only on ALICE, so nothing in this list can be regenerated locally.

Status key: **stale** = must be regenerated; **verified unaffected** = checked and the fix
provably does not change these values.

---

## ichorCNA

### `calculate_FGA_AS.R` — arm direction at the 25% threshold

Fixed in `fix(ichorcna): Call arm direction by larger altered fraction`.

Arms where both the gained and lost fraction exceeded 25% were labelled GAIN regardless of
which was larger. Direction is now taken from the larger fraction.

| File | Status |
|---|---|
| `results/ichorCNA/FGA_AS_corrected/arm_level_scores_detail.csv` | **stale** — 7 of 780 rows change |
| `results/ichorCNA/FGA_AS_corrected/FGA_AS_summary.csv` | **stale** — `arms_gained_25` / `arms_lost_25` only |
| `results/ichorCNA/FGA_AS_corrected/AS_threshold_sensitivity_comparison.csv` | **stale** — carries `arms_gained_25` / `arms_lost_25`; its `AS_25`, `AS_difference` and `direction` columns are unaffected |
| `results/ichorCNA/FGA_AS_corrected/AS_sensitivity_comparison.png` | verified unaffected — plots `AS_50` vs `AS_25` and `direction` only |
| FGA columns, `AS_50`, `AS_25` totals, all correlations | verified unaffected |

The seven affected rows, all changing GAIN to LOSS:

| Sample | Arm | frac_gain | frac_loss |
|---|---|---|---|
| B1_S01 | 13q | 0.2688 | 0.7312 |
| B2_S05 | 3q  | 0.4106 | 0.5866 |
| B2_S08 | 17p | 0.3361 | 0.5809 |
| B3_S01 | 9q  | 0.2516 | 0.7443 |
| B3_S08 | 6p  | 0.3826 | 0.4783 |
| B3_S08 | 19q | 0.3702 | 0.4010 |
| B3_S09 | 9q  | 0.3459 | 0.3774 |

`scored_25` is identical for all 780 rows, so the aneuploidy score totals are unchanged.
Only the gained-versus-lost breakdown needs regenerating.

To re-run: `sbatch scripts/ichorCNA/run_FGA_AS.slurm` (or run `calculate_FGA_AS.R`
directly), then copy `FGA_AS_summary.csv` and `arm_level_scores_detail.csv` back into
`results/ichorCNA/FGA_AS_corrected/`.

---

## QDNAseq

No committed results are affected. The Stage 1 fixes touched failure-path reporting, the
per-sample plot directory, the raw-count plot axis, and the package installer. QDNAseq
diagnostic plots live on ALICE and are not tracked; regenerate them at leisure to pick up
the un-clipped y-axis.

---

## Data integrity issues — not code fixes, but blocking

These are not fixed by any commit. They need checking against the source data on ALICE.

**Six samples have `long_fragment_fraction` of exactly 0** in
`results/fragmentomics/tables/fragmentomics_metrics_summary.csv` (B1_S02, B2_S01, B2_S02,
B2_S03, B2_S05, B2_S07), over 7.5–11.3M fragments each. Zero fragments above 250bp is not
biologically achievable in plasma cfDNA, which always carries a di-nucleosomal population
near 330bp. The likely cause is fragment histogram files truncated near 250bp, which would
be invisible in `n_fragments` and `median_fragment_length`.

This gates the interpretation of B2_S11, whose flag rests on a long-fragment fraction 59x
the next highest — a ratio against a possibly artefactual baseline.

Diagnostic, on ALICE:

```
awk 'NR>1{if($1>m)m=$1}END{print FILENAME, m, NR}' fragment_histograms/*.tsv
```

If the six samples show `max(fragment_length)` near 250, the histograms are truncated and
every long-fragment result needs regenerating from the BAMs.

**`fig_batch3_longfragfraction_dedup_sensitivity.png` has no source script.** No script in
the repository writes that filename. It is a committed thesis figure with no provenance,
for the metric B2_S11 is flagged on.

---

## Reporting issues in committed output

**`results/ichorCNA/qdnaseq_concordance/concordance_statistics.csv`** reports
`proportional_bias_pvalue` as exactly `0`. No test returns a true zero; this is underflow
or rounding. Report as `< 2.2e-16`, or carry more precision, before it reaches a thesis
table.
