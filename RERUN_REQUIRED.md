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

## RASCAL

### `run_rascal_batch3.R` — cellularity compared against tumour fraction

Fixed in `fix(rascal): Compare implied tumour fraction, not cellularity`.

`tf_difference` subtracted ichorCNA tumour fraction from RASCAL cellularity directly. These
coincide only at ploidy 2; five samples fit at 1.50 and four above 3.3.

| File | Status |
|---|---|
| `results/RASCAL/RASCAL_batch3_summary.csv` | **stale** — `tf_difference` for all 13; gains `rascal_expected_tf` and `fallback_used` columns |
| `tf_agreement` labels | verified unaffected — all 13 were and remain DISCORDANT |
| `results/RASCAL/RASCAL_batch3_constrained_summary.csv` | already correct — it applied the conversion |

This fix **resolves a contradiction between two committed files**. They currently disagree
on the same fits:

| Sample | unconstrained CSV | constrained CSV (`unc_tf_diff`) | corrected |
|---|---|---|---|
| B3_S06 | 0.337 | 0.318 | 0.319 |
| B3_S08 | 0.317 | 0.299 | 0.299 |
| B3_S09 | 0.379 | 0.384 | 0.384 |

After the fix the unconstrained file reproduces the constrained one exactly.

`fallback_used` is a new column with no historical equivalent, so which existing rows came
from a genuine fit and which from the grid-minimum fallback cannot be recovered without
re-running.

### `run_rascal_constrained_batch3.R` — interpretation label

Fixed in `fix(rascal): Label improved-but-discordant fits accurately`.

| File | Status |
|---|---|
| `results/RASCAL/RASCAL_batch3_constrained_summary.csv` | **stale** — `interpretation` for B3_S09 only |
| all numeric columns | verified unaffected |

B3_S09 reads `IMPROVED — constrained better than unconstrained`; it should indicate that it
remains discordant, being 0.199 from ichorCNA against an acceptance band of 0.10. Editable
by hand.

### Not fixed — would change the analysis, not correct it

**The ploidy floor binds.** `min_ploidy = 1.5` in both RASCAL scripts. Five of 13 samples
return a fitted ploidy of exactly 1.50 (B3_S06, B3_S07, B3_S08, B3_S10, B3_S12), meaning
the optimum lies at or beyond the search boundary. RASCAL's own default lower bound is
1.25. These five ploidy values should not be quoted as estimates. Widening the floor
changes every fit.

**The ±0.15 constraint is applied to the wrong variable.** In the constrained script the
band is derived from tumour fraction but used as bounds on cellularity. B3_S09 searched
cellularity in [0.21, 0.51] and ended 0.199 from ichorCNA under a constraint labelled 0.15.
Correcting it means back-solving the bounds through the conversion at each candidate
ploidy, which changes the searched space and so is a new analysis.

**The constrained run left ploidy unconstrained.** It was motivated by the ploidy boundary
artefact but restricts only cellularity; the fitted ploidies move off the floor as a side
effect rather than by design.

---

## QDNAseq

No committed results are affected. The Stage 1 fixes touched failure-path reporting, the
per-sample plot directory, the raw-count plot axis, and the package installer. QDNAseq
diagnostic plots live on ALICE and are not tracked; regenerate them at leisure to pick up
the un-clipped y-axis.

---

## Evidence gaps — nothing to regenerate, but arguments that cannot currently be checked

**The B3_S12 exclusion cannot be verified from this repository.** B3_S12 is excluded from
the CNA cohort for ploidy non-identifiability. Its main run fits ploidy 2.831 with a 0.436
subclone fraction, and two diagnostic re-runs exist to test whether a near-diploid solution
explains the data comparably (`run_ichorCNA_B3_S12_diploid.slurm`,
`run_ichorCNA_B3_S12_noSubclone.slurm`).

Those runs completed, but **their params and seg outputs are not committed**, and the
captured logs echo only the input parameters. Nothing in the repository records what ploidy
or tumour fraction they produced.

Note also that `--ploidy "c(2)"` sets the *starting* values for ichorCNA's ploidy search,
not a hard constraint, so a run labelled "diploid" still estimates its own ploidy. The
argument needs the fitted output, not the requested parameter.

To close this: copy `B3_S12_diploid/` and `B3_S12_noSubclone/` params files from ALICE into
`results/ichorCNA/params/` under distinguishing names. The same applies to the B2_S08 and
B2_S15 diagnostic runs, which underpin those two curation decisions.

**The expected-logR screen is mis-calibrated for high-ploidy samples.** The formula in
`ichorCNA_anomaly_scan.R` is exact only at ploidy 2 and becomes progressively permissive as
ploidy departs from it. Twenty of 34 samples sit more than 0.1 from ploidy 2, including
B2_S15 at 2.598. Left unchanged deliberately, since it drives a hand-reviewed triage step
and correcting it would retrospectively alter the basis of curation already performed — it
should be revised alongside a re-review, not on its own. See the comment at the function.

**The gene mapping does not record its Ensembl release.** `CNA_gene_mapping.R` queries
biomaRt without pinning a version, so re-running can yield different coordinates from
identical input. The build is correct for hg38; only the version is unrecorded.

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

### `plot_qdnaseq_ichorcna_concordance_v2.R` — p-value underflow

Fixed in `fix(ichorcna): Stop reporting a p-value of exactly zero`.

| File | Status |
|---|---|
| `results/ichorCNA/qdnaseq_concordance/concordance_statistics.csv` | **stale** — one cell |
| everything else in that file | verified unaffected |

`proportional_bias_pvalue` currently reads `0`; it should read `2.22e-16` and be quoted in
the thesis as `p < 2.2e-16`. The slope, bias, limits of agreement, correlations and counts
are all unchanged, so this can be corrected by hand rather than re-running if preferred.
