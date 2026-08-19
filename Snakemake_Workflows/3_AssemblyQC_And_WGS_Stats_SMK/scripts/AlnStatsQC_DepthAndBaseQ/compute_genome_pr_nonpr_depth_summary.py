"""
Genome-wide / Paralogous-Region-only / non-PR-only pooled mean+median depth, per sample,
across both cohorts (Mtb151 + TBP22CI) -- the new primary depth-summary method, using
`mosdepth`'s whole-genome per-base output plus `bedtools intersect` to clip it to the
exact PR / non-PR region-set boundaries. See docs/AlnStatsQC_Overview.md for the full
rationale and why this replaces two earlier, differently-confounded methods (an
arithmetic bp-weighted-subtraction approach, and a bedtools/CRAM-partitioned approach).

Method, per sample:
  1. Load {sample}.per-base.bed.gz (RLE: chrom, start, end, depth) -- ALL reads, no
     exclusion by whether a read also touches a region elsewhere.
  2. Genome-wide: pooled weighted mean + median directly over the whole file.
  3. PR-only: `bedtools intersect -a <per-base.bed.gz> -b <PR bed>` clips the RLE
     segments to exactly the PR regions' positions; pooled weighted mean + median over
     the clipped output.
  4. Non-PR-only: same, intersecting against the PR-complement bed instead.
  5. Weighted median: sort (depth, length) pairs by depth, walk cumulative length to the
     halfway point of total length.

Cross-check performed automatically: genome-wide pooled mean must match
{sample}.mosdepth.summary.txt's own reported mean (both computed by different code
paths -- mosdepth's C++ internals vs. this script's own aggregation over its output).

Run with the nucflag_env Python (has pandas) and bedtools on PATH:
  /n/data1/hms/dbmi/farhat/mm774/Projects/Mtb-GeneConv/NucFlag_Testing/envs/nucflag_env/bin/python \
    compute_genome_pr_nonpr_depth_summary.py
"""

import gzip
import subprocess
from pathlib import Path

import pandas as pd

from alnstats_common import (
    DATASETS, discover_samples_in, mosdepth_per_base_path, mosdepth_summary_path,
    load_pr_bed_for, load_nonpr_bed_for,
)

OUT_DIR = Path(__file__).parent.parent / "results"
OUT_DIR.mkdir(exist_ok=True)


def weighted_mean_median(rows):
    """rows: iterable of (depth, length). Returns (mean, median, total_length)."""
    total_len = 0
    total_bp = 0
    segs = []
    for d, L in rows:
        total_len += L
        total_bp += d * L
        segs.append((d, L))
    if total_len == 0:
        return float("nan"), float("nan"), 0
    mean = total_bp / total_len
    segs.sort(key=lambda x: x[0])
    half = total_len / 2
    cum = 0
    median = None
    for d, L in segs:
        cum += L
        if cum >= half:
            median = d
            break
    return mean, median, total_len


def genome_wide_stats(per_base_path):
    rows = []
    with gzip.open(per_base_path, "rt") as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            s, e, d = int(p[1]), int(p[2]), int(p[3])
            rows.append((d, e - s))
    return weighted_mean_median(rows)


def region_clipped_stats(per_base_path, region_bed_path):
    proc = subprocess.run(
        ["bedtools", "intersect", "-a", str(per_base_path), "-b", str(region_bed_path)],
        capture_output=True, text=True, check=True,
    )
    rows = []
    for line in proc.stdout.splitlines():
        if not line.strip():
            continue
        p = line.split("\t")
        s, e, d = int(p[1]), int(p[2]), int(p[3])
        rows.append((d, e - s))
    return weighted_mean_median(rows)


def read_mosdepth_summary_mean(summary_path):
    df = pd.read_csv(summary_path, sep="\t")
    total_row = df[df["chrom"] == "total"]
    return float(total_row["mean"].iloc[0])


def main():
    rows = []
    for cohort, results_dir in DATASETS.items():
        samples = discover_samples_in(results_dir)
        pr_bed = load_pr_bed_for(results_dir)
        nonpr_bed = load_nonpr_bed_for(results_dir)

        for sample in samples:
            per_base = mosdepth_per_base_path(results_dir, sample)

            genome_mean, genome_median, genome_len = genome_wide_stats(per_base)
            pr_mean, pr_median, pr_len = region_clipped_stats(per_base, pr_bed)
            nonpr_mean, nonpr_median, nonpr_len = region_clipped_stats(per_base, nonpr_bed)

            known_mean = read_mosdepth_summary_mean(mosdepth_summary_path(results_dir, sample))
            assert abs(genome_mean - known_mean) < 0.01, (
                f"{sample}: pooled genome mean {genome_mean} doesn't match "
                f"mosdepth.summary.txt's mean {known_mean}"
            )

            rows.append(dict(
                SampleID=sample,
                Cohort=cohort,
                Genome_MeanDepth=genome_mean,
                Genome_MedianDepth=genome_median,
                Genome_Length=genome_len,
                PR_MeanDepth=pr_mean,
                PR_MedianDepth=pr_median,
                PR_Length=pr_len,
                NonPR_MeanDepth=nonpr_mean,
                NonPR_MedianDepth=nonpr_median,
                NonPR_Length=nonpr_len,
            ))

    per_sample = pd.DataFrame(rows)
    out_path = OUT_DIR / "AllSamples_Genome_PR_NonPR_Depth_Summary.tsv"
    per_sample.to_csv(out_path, sep="\t", index=False)

    def fmt(series):
        return f"mean={series.mean():.3f}, median={series.median():.3f}, min={series.min():.3f}, max={series.max():.3f}"

    lines = []
    lines.append("Genome-wide / PR-only / non-PR-only pooled depth summary -- both cohorts")
    lines.append("=" * 70)
    for cohort in DATASETS:
        sub = per_sample[per_sample["Cohort"] == cohort]
        n = len(sub)
        lines.append(f"\n{cohort} (n={n})")
        lines.append(f"  Genome mean depth:  {fmt(sub['Genome_MeanDepth'])}")
        lines.append(f"  Genome median depth: {fmt(sub['Genome_MedianDepth'])}")
        lines.append(f"  PR mean depth:      {fmt(sub['PR_MeanDepth'])}")
        lines.append(f"  PR median depth:    {fmt(sub['PR_MedianDepth'])}")
        lines.append(f"  NonPR mean depth:   {fmt(sub['NonPR_MeanDepth'])}")
        lines.append(f"  NonPR median depth: {fmt(sub['NonPR_MedianDepth'])}")
        mean_ratio = sub["PR_MeanDepth"] / sub["NonPR_MeanDepth"]
        median_ratio = sub["PR_MedianDepth"] / sub["NonPR_MedianDepth"]
        lines.append(f"  PR/NonPR mean-depth ratio:   {fmt(mean_ratio)}")
        lines.append(f"  PR/NonPR median-depth ratio: {fmt(median_ratio)}")
        n_pr_lower_mean = int((mean_ratio < 1).sum())
        n_pr_lower_median = int((median_ratio < 1).sum())
        lines.append(f"  Samples where PR mean depth < NonPR mean depth: {n_pr_lower_mean} / {n}")
        lines.append(f"  Samples where PR median depth < NonPR median depth: {n_pr_lower_median} / {n}")

    summary_text = "\n".join(lines)
    (OUT_DIR / "AllSamples_Genome_PR_NonPR_Depth_Summary.txt").write_text(summary_text + "\n")
    print(summary_text)
    print(f"\nWrote {out_path}")


if __name__ == "__main__":
    main()
