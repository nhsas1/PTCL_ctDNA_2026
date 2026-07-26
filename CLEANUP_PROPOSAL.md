# Repository cleanup proposal

**Nothing in this document has been executed.** Renames touch paths in many scripts, so
every item here needs sign-off before anything moves. Items are ordered by value-to-risk:
the early ones are additive and safe, the later ones touch working code.

---

## 1. README structure — proposed

The repository has a `README.md` but no orientation for a reader arriving cold. A marker
assessing implementability needs to see the pipeline order and the environment without
reading 90 scripts. Proposed sections:

**Project** — one paragraph: sWGS of plasma cfDNA in PTCL, 34 samples, thesis title.

**Cohort** — the three batches (n=5, 16, 13), coverage range 0.87–1.77x, hg38/BWA-MEM, and
the Batch 3 deduplication asymmetry stated once, prominently, since it explains flag
differences throughout.

**The two analysis arms** — CNA profiling and fragmentomics, with the key distinction that
the CNA arm runs on the curated n=20 subset while fragmentomics uses all 34.

**Cohort derivation** — the exact arithmetic, because this is the number most likely to be
questioned:

```
34 total
-13 below the 3% tumour fraction detection floor
-1  excluded for ploidy non-identifiability (B3_S12)
=20 analysable for CNA
```

Worth noting explicitly that the "11 below floor + 2 borderline" phrasing used elsewhere
describes the same 13 samples, since the two framings currently look like a discrepancy.

**Pipeline execution order** — a numbered list naming the actual scripts, since the
dependency order is not obvious from the filenames. Particularly: step 5c must not be run
after step 5, and step 5e must be run after step 5d.

**Environment** — module versions (samtools/1.18, Picard 3.0.0, R/4.3.1, GISTIC2 via
Apptainer) and the note that `install_packages.R` bootstraps the R library but `rascal` must
be installed separately.

**Stated limitations** — mappability correction absent, Batch 3 duplicates retained in the
CNA arm, the two MAPQ thresholds (37 in QDNAseq, 20 in fragmentomics), coverage-only ploidy
estimation, and the reproducibility gaps listed in section 5 below.

**Results map** — which directory holds what, and a pointer to `RERUN_REQUIRED.md`.

---

## 2. `.gitignore` audit — three findings, all minor

The large-file rules are doing their job. During this review roughly 60 `.ORIGINAL` backup
files were created and **none** appeared as untracked, confirming the backup pattern is
correctly covered.

**a. Two rules are duplicated.** `*.bam` and `*.bai` each appear twice (lines 3–4 and
35–36). Harmless, but the second block is labelled "ichorCNA large outputs" and doesn't need
them.

**b. Three `.ORIGINAL` files are tracked despite the `*.ORIGINAL` rule.**

```
results/ichorCNA/corrected_seg_files/originals_backup/B2_S08.seg.txt.ORIGINAL
results/ichorCNA/corrected_seg_files/originals_backup/B2_S15.seg.txt.ORIGINAL
results/ichorCNA/corrected_seg_files/originals_backup/B3_S09.seg.txt.ORIGINAL
```

They were committed before the rule existed, and gitignore does not apply retroactively to
tracked files.

**Recommendation: keep them tracked, and add a negation rule.** These are not stray
backups — they are the pre-curation state of the three manually corrected samples, which is
exactly the provenance evidence a marker would want when assessing whether curation was
justified. They total a few KB. Suggested:

```gitignore
*.ORIGINAL
!results/ichorCNA/corrected_seg_files/originals_backup/*.ORIGINAL
```

This makes the intent explicit rather than leaving them as an accident of commit order.

**c. Nothing else is mis-tracked.** Every tracked file was checked against the ignore rules;
these three are the only cases.

---

## 3. Committed SLURM logs — 24 files, 272 KB

`scripts/ichorCNA/logs/` holds 24 `.out`/`.err` files from the anomaly scan, the B2_S08 /
B2_S15 / B3_S12 diagnostic runs, and the concordance jobs.

**Recommendation: keep them, but they are not currently doing the job they appear to do.**
As established in Stage 2, these logs echo *input parameters only* — they do not record what
ploidy or tumour fraction the diagnostic runs produced. So they document that a run happened
without documenting its outcome, which is precisely the evidence the B3_S12 exclusion needs.

The higher-value action is not deleting them but **adding the missing params files** (see
section 5). If space ever matters, `slurm_*.log` is already ignored, so extending that
pattern would be trivial — but 272 KB is not the problem here.

---

## 4. Naming and structure — the substantive proposals

These touch paths. Each is listed with what breaks.

### 4a. `export_rascal_input*` is misnamed (recommended)

RASCAL does not read these files. Both RASCAL scripts read the QDNAseq `.seg`/`.igv` files
directly. The only consumers are the two concordance scripts.

| Current | Proposed |
|---|---|
| `scripts/RASCAL/export_rascal_input.R` | `scripts/QDNAseq/export_binlevel_for_concordance.R` |
| `scripts/QDNAseq/export_rascal_input_batch3.R` | `scripts/QDNAseq/export_binlevel_for_concordance_batch3.R` |
| output dir `RASCAL/input/` | `concordance/qdnaseq_bins/` |
| SLURM job name `RASCAL_export` | `concordance_export` |

**Touches:** 2 R scripts, 2 SLURM wrappers, and the input paths in both
`plot_qdnaseq_ichorcna_concordance*.R`. Six files. Already documented in-place, so this is
cosmetic-but-clarifying rather than urgent.

### 4b. Batch 1–2 and Batch 3 counterparts live in different directories

`export_rascal_input.R` sits in `scripts/RASCAL/` while its Batch 3 twin sits in
`scripts/QDNAseq/`. Moving the Batch 1–2 one to `scripts/QDNAseq/` puts the pair together.
Folds naturally into 4a.

### 4c. Two project roots are in use

258 path references use `/scratch/alice/n/nhsas1/PTCL`, 43 use
`/scratch/alice/n/nhsas1/ptcl_ctdna_thesis`. The split is not by pipeline — it cuts across
ichorCNA and fragmentomics.

**Recommendation: settle on one and define it once per script.** The deeper fix is that
every script hardcodes absolute paths; a single `PROJECT_ROOT` variable at the top of each
would make the whole repo portable and make this class of drift impossible. That is a larger
change and should be a deliberate decision, not a side effect of cleanup.

### 4d. Five different SLURM log directories

```
11  ptcl_ctdna_thesis/scripts/ichorCNA/logs
13  PTCL/fragmentomics/logs
 5  PTCL/ichorCNA/logs
 6  PTCL/scripts/logs
 2  PTCL/scripts/RASCAL/logs
```

All must pre-exist or SLURM rejects the job at submission. **Recommendation: one
`logs/<stage>/` tree.** Touches every `.slurm` file but only the two `--output`/`--error`
lines in each, so it is mechanical and low-risk.

### 4e. Deployed script paths are flat for QDNAseq, nested elsewhere

Six wrappers call `scripts/<x>.R` while twenty call `scripts/<stage>/<x>.R`. The repository
is nested throughout, so the QDNAseq deployment does not match the repository layout.
**Recommendation: make QDNAseq nested to match.** Touches 6 wrapper lines.

### 4f. `step5c` overwrites `step5`'s output filenames

Re-running step 5 silently restores the known-buggy lexicographic chromosome ordering.
`step5d` handles this correctly with a `_final` suffix. **Recommendation: give step 5c's
output a `_chrfix` suffix**, or fold the fix into step 5 and retire 5c. Now documented in
both scripts, so the trap is at least visible.

Related: there is no `step5b` in the repository, so the numbering 5, 5c, 5d, 5e, 5f refers to
a state a reader cannot reconstruct.

### 4g. Two concordance scripts write identical filenames

`plot_qdnaseq_ichorcna_concordance.R` (v1) overwrites v2's committed figures with an older
version lacking the Bland-Altman panel. Marked superseded in Stage 2. **Recommendation:
either delete v1 or move it to an `archive/` subdirectory** — a superseded script that
silently clobbers current results is a live hazard, not provenance.

### 4h. Metadata files live in two places

`sample_metadata.csv` sits at the project root, `sample_metadata_batch3.csv` under
`scripts/`. Same kind of file, two conventions. Low priority.

---

## 5. Reproducibility gaps — worth a README entry even if unfixable

These are not cleanup items but belong in the same conversation.

**`hg38_bins_15kb_annotated.rds` has no generating script.** Six scripts depend on it. This
is the largest reproducibility hole in the CNA arm — the bin definitions, GC annotation and
blacklist flags underpinning every QDNAseq result cannot be regenerated from this
repository.

**`CNA_gene_mapping.R` does not pin an Ensembl release.** The build is correct for hg38, but
re-running can yield different coordinates. Record the release used, or pin it.

**Diagnostic run outputs are not committed.** The B2_S08, B2_S15 and B3_S12 alternative-run
params and seg files exist on ALICE but not here, so three curation decisions and one
exclusion cannot be verified from the repository. **This is the highest-value gap to close**
and needs only a file copy.

**`fig_batch3_longfragfraction_dedup_sensitivity.png` has no source script.** A committed
thesis figure with no provenance, for the metric B2_S11 is flagged on.

---

## 6. Suggested order

1. Commit `CLAUDE.md` (done alongside this proposal — it was untracked).
2. Add the README (additive, no risk).
3. `.gitignore` negation rule for the curation backups, drop the two duplicate lines.
4. Copy the missing diagnostic params files off ALICE.
5. Retire or archive concordance v1, and fix the step5c filename collision — these two are
   live hazards rather than tidiness.
6. Everything in section 4c–4e — path and log consolidation — as one deliberate pass,
   ideally with `PROJECT_ROOT` variables, after the thesis numbers are settled.

Items 1–5 are safe to do now. Item 6 changes paths in every script and is better done once,
after re-running, than piecemeal before.
