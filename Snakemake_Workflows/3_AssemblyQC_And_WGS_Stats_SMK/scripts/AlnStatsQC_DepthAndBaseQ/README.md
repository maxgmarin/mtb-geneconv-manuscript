# AlnStatsQC: depth and base-quality comparison (Paralogous Regions vs. genome-wide)

Documented copies of the two scripts used to compare long-read alignment depth and base
quality inside H37Rv Paralogous Regions against the rest of the genome, for both the
Mtb151 and TBP22CI cohorts. Original source and full results/docs:
`/n/data1/hms/dbmi/farhat/mm774/Projects/Mtb-GeneConv/AlnStatsQC/`.

- **`compute_genome_pr_nonpr_depth_summary.py`** -- genome-wide / Paralogous-Region-only /
  non-PR-only pooled mean+median depth per sample, from `mosdepth`'s whole-genome
  per-base output (`{sample}.per-base.bed.gz`, produced by this pipeline's earlier
  coverage-stats exploration -- not part of the core `AsmAndWGS.QC.smk` rules copied
  here) clipped to each region set via `bedtools intersect`.
- **`compute_baseq_pr_vs_genome_summary.py`** -- genome-wide mean base quality (from
  `samtools coverage`'s output on the LR-to-H37Rv alignment) vs. the mean/median of
  per-PR-region mean base quality. Flags (but doesn't exclude) samples whose reads carry
  a placeholder all-`!` (Q0) quality string -- a real property of some samples' input
  reads, not a pipeline bug.
- **`alnstats_common.py`** -- shared config/loader module both scripts import from.

