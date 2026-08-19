# NucFlag results summarization

Documented copies of the scripts used to summarize this pipeline's own misassembly-
calling output across the Mtb151 (151 samples) and TBP22CI (22 samples) cohorts.

- **`nucflag_common.py`** -- shared config/loader module: discovers samples, loads
  per-sample misassembly BEDs, H37Rv-liftover BEDs, and the two curated region-set BEDs
  (NucDivHotspots windows, Paralogous Regions), plus an interval-overlap helper.
- **`analyze_nucflag_results.py`** -- per-cohort stats (total misassemblies,
  samples with >=1, overlap with NucDivHotspots windows, overlap with Paralogous
  Regions). Writes `{cohort}_per_sample_nucflag_summary.tsv` (copied into `results/`
  here) and a `.txt` narrative summary (not copied here).
- **`merge_all_misassemblies.py`** -- one call-level TSV per cohort with every
  misassembly call (lifted or not, via a `Lifted` boolean column), for deeper/notebook-
  driven analysis beyond the headline counts. Writes `{cohort}_all_misassemblies_merged.tsv`
  (copied into `results/` here).

## Copied result tables (`results/`)

- `Mtb151_per_sample_nucflag_summary.tsv`, `TBP22CI_per_sample_nucflag_summary.tsv`
- `Mtb151_all_misassemblies_merged.tsv`, `TBP22CI_all_misassemblies_merged.tsv`
