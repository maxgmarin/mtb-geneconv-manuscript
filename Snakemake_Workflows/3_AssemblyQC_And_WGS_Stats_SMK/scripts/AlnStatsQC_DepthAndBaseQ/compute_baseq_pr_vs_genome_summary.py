"""
Genome-wide vs. Paralogous-Region mean base quality, per sample, across both cohorts
(Mtb151 + TBP22CI) -- uses only already-generated, already-validated-as-real samtools
coverage outputs (BAM-based, not the CRAM-based InPRsOnly/NotInPRs files, which have a
broken meanbaseq=255 sentinel -- see docs/AlnStatsQC_Overview.md).

Method, per sample:
  1. Genome-wide meanbaseq: single value from {sample}.LR.AlnToH37Rv.coverage.genome.tsv
  2. Per-PR-region meanbaseq: ~200 values from
     {sample}.LR.AlnToH37Rv.coverage.paralogous_regions.tsv, aggregated by plain mean and
     median across regions.
  3. Compare PR aggregate to genome-wide value directly.

Known data-quality caveat (found while building this script, confirmed via raw BAM
inspection): a subset of Mtb151 samples have Genome_MeanBaseQ == 0 -- NOT a pipeline bug,
but because every base in that sample's BAM carries a placeholder quality string of all
'!' (Phred+33 = Q0), i.e. the input long reads simply lack real per-base quality scores.
These samples are flagged via HasZeroBaseQ and excluded from the printed cohort-level
summary stats (still included in the per-sample TSV, clearly flagged) so they don't
silently drag down the aggregate comparison.

Run with the nucflag_env Python (has pandas):
  /n/data1/hms/dbmi/farhat/mm774/Projects/Mtb-GeneConv/NucFlag_Testing/envs/nucflag_env/bin/python \
    compute_baseq_pr_vs_genome_summary.py
"""

from pathlib import Path

import pandas as pd

from alnstats_common import DATASETS, discover_samples_in, load_genome_coverage_for, load_pr_coverage_for

OUT_DIR = Path(__file__).parent.parent / "results"
OUT_DIR.mkdir(exist_ok=True)


def main():
    rows = []
    for cohort, results_dir in DATASETS.items():
        for sample in discover_samples_in(results_dir):
            genome_row = load_genome_coverage_for(results_dir, sample)
            pr_df = load_pr_coverage_for(results_dir, sample)

            genome_baseq = float(genome_row["meanbaseq"])
            pr_baseq_mean = float(pr_df["meanbaseq"].mean())
            pr_baseq_median = float(pr_df["meanbaseq"].median())

            rows.append(dict(
                SampleID=sample,
                Cohort=cohort,
                Genome_MeanBaseQ=genome_baseq,
                PR_MeanBaseQ_MeanOfRegions=pr_baseq_mean,
                PR_MeanBaseQ_MedianOfRegions=pr_baseq_median,
                PR_minus_Genome_MeanBaseQ=pr_baseq_mean - genome_baseq,
                PR_minus_Genome_MedianBaseQ=pr_baseq_median - genome_baseq,
                N_PR_Regions=len(pr_df),
                HasZeroBaseQ=(genome_baseq == 0),
            ))

    per_sample = pd.DataFrame(rows)
    out_path = OUT_DIR / "AllSamples_BaseQ_PR_vs_Genome_Summary.tsv"
    per_sample.to_csv(out_path, sep="\t", index=False)

    def fmt(series):
        return f"mean={series.mean():.3f}, median={series.median():.3f}, min={series.min():.3f}, max={series.max():.3f}"

    lines = []
    lines.append("Genome-wide vs. per-PR-region mean base quality -- both cohorts")
    lines.append("=" * 70)
    n_flagged_total = int(per_sample["HasZeroBaseQ"].sum())
    if n_flagged_total:
        lines.append(
            f"\nNOTE: {n_flagged_total} sample(s) have Genome_MeanBaseQ == 0 -- confirmed "
            f"via raw BAM inspection to be a placeholder all-'!' (Q0) quality string in the "
            f"input reads, not a pipeline bug. These are flagged (HasZeroBaseQ=True) in the "
            f"TSV but excluded from the cohort summary stats below."
        )
        flagged = per_sample[per_sample["HasZeroBaseQ"]]
        for cohort in DATASETS:
            ids = flagged[flagged["Cohort"] == cohort]["SampleID"].tolist()
            if ids:
                lines.append(f"  {cohort}: {', '.join(ids)}")
    for cohort in DATASETS:
        sub = per_sample[(per_sample["Cohort"] == cohort) & (~per_sample["HasZeroBaseQ"])]
        n = len(sub)
        lines.append(f"\n{cohort} (n={n}, excluding zero-baseq samples)")
        lines.append(f"  Genome meanbaseq:                {fmt(sub['Genome_MeanBaseQ'])}")
        lines.append(f"  PR meanbaseq (mean of regions):   {fmt(sub['PR_MeanBaseQ_MeanOfRegions'])}")
        lines.append(f"  PR meanbaseq (median of regions): {fmt(sub['PR_MeanBaseQ_MedianOfRegions'])}")
        lines.append(f"  PR-mean minus genome:   {fmt(sub['PR_minus_Genome_MeanBaseQ'])}")
        lines.append(f"  PR-median minus genome: {fmt(sub['PR_minus_Genome_MedianBaseQ'])}")
        n_pr_lower_mean = int((sub['PR_minus_Genome_MeanBaseQ'] < 0).sum())
        n_pr_lower_median = int((sub['PR_minus_Genome_MedianBaseQ'] < 0).sum())
        lines.append(f"  Samples where PR mean baseq < genome-wide baseq: {n_pr_lower_mean} / {n}")
        lines.append(f"  Samples where PR median baseq < genome-wide baseq: {n_pr_lower_median} / {n}")

    summary_text = "\n".join(lines)
    (OUT_DIR / "AllSamples_BaseQ_PR_vs_Genome_Summary.txt").write_text(summary_text + "\n")
    print(summary_text)
    print(f"\nWrote {out_path}")


if __name__ == "__main__":
    main()
