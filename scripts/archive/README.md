# Archived scripts

Superseded scripts kept for provenance. **Nothing here should be run.**

Each was replaced by a newer version that writes the same output filenames, so running an
archived script would silently overwrite current committed results with an older analysis.

| Script | Superseded by | Why |
|---|---|---|
| `plot_qdnaseq_ichorcna_concordance.R` | `scripts/ichorCNA/plot_qdnaseq_ichorcna_concordance_v2.R` | v2 adds the Bland-Altman agreement analysis and `concordance_statistics.csv`. v1 writes the same figure and per-sample CSV filenames, so running it replaces the committed figures with a two-panel version reporting correlation only. |
| `run_concordance.slurm` | `scripts/ichorCNA/run_concordance_v2.slurm` | Wrapper for the above. |

The deployed paths inside these files still point at their original locations and have
deliberately not been updated, since running them is not intended.
