"""
Manuscript-critical NucFlag misassembly analysis -- generalized to run over both
datasets (Mtb151, TBP22CI). Supersedes the earlier Mtb151-only version of this script
(same method/logic, just parameterized by cohort via nucflag_common.DATASETS).

Answers, per dataset, directly against live pipeline output paths:
  Q1. How many misassembly errors does NucFlag find across the entire dataset?
  Q2. How many samples have ANY assembly error flagged? Which ones?
  Q3. How many of those errors overlap a NucDiv hotspot window?
  Q4. How many overlap a Paralogous Region (PR) window?

Method (deliberately simple, per 2026-08-05 decision -- see
docs/NucFlag_Pipeline_Overview.md):
  - Q1/Q2: count non-"correct" rows in each sample's {sample}.misassemblies.bed.
  - Q3/Q4: of the misassemblies that successfully lifted over to H37Rv coordinates
    (misassemblies_liftover_to_H37Rv.bed), directly test each against the two
    reference region BEDs via a plain half-open-interval overlap. Misassemblies that
    did NOT lift over are excluded from Q3/Q4 entirely.

Run with the nucflag_env Python (has pandas):
  /n/data1/hms/dbmi/farhat/mm774/Projects/Mtb-GeneConv/NucFlag_Testing/envs/nucflag_env/bin/python \
    analyze_nucflag_results.py
"""

from pathlib import Path
import pandas as pd

from nucflag_common import (
    DATASETS, discover_samples, load_misassemblies, load_liftover,
    load_region_bed, overlapping_call_mask,
)

OUT_DIR = Path(__file__).parent.parent / "results"
OUT_DIR.mkdir(exist_ok=True)


def analyze_dataset(cohort, results_dir):
    samples = discover_samples(results_dir)
    hotspot_bed = load_region_bed("hotspot", results_dir)
    pr_bed = load_region_bed("pr", results_dir)

    rows = []
    for sample in samples:
        mis = load_misassemblies(sample, results_dir)
        flagged = mis[mis["name"] != "correct"]
        n_misassemblies = len(flagged)
        type_counts = flagged["name"].value_counts().to_dict()

        lift = load_liftover(sample, results_dir)
        n_lifted = len(lift)
        n_not_lifted = n_misassemblies - n_lifted

        hotspot_hit = overlapping_call_mask(lift, hotspot_bed)
        pr_hit = overlapping_call_mask(lift, pr_bed)
        n_overlap_hotspot = int(hotspot_hit.sum())
        n_overlap_pr = int(pr_hit.sum())

        rows.append(dict(
            SampleID=sample,
            Cohort=cohort,
            n_misassemblies=n_misassemblies,
            misassembly_type_counts=type_counts,
            n_misassemblies_lifted_to_h37rv=n_lifted,
            n_misassemblies_not_lifted=n_not_lifted,
            n_misassemblies_overlapping_hotspot=n_overlap_hotspot,
            n_misassemblies_overlapping_pr=n_overlap_pr,
        ))

    per_sample = pd.DataFrame(rows)
    per_sample.to_csv(OUT_DIR / f"{cohort}_per_sample_nucflag_summary.tsv", sep="\t", index=False)

    n_total = len(per_sample)
    n_with_any_error = int((per_sample["n_misassemblies"] > 0).sum())
    samples_with_error = per_sample.loc[per_sample["n_misassemblies"] > 0, "SampleID"].tolist()

    total_misassemblies = int(per_sample["n_misassemblies"].sum())
    total_lifted = int(per_sample["n_misassemblies_lifted_to_h37rv"].sum())
    total_not_lifted = int(per_sample["n_misassemblies_not_lifted"].sum())

    total_overlap_hotspot = int(per_sample["n_misassemblies_overlapping_hotspot"].sum())
    total_overlap_pr = int(per_sample["n_misassemblies_overlapping_pr"].sum())

    n_samples_hotspot_overlap = int((per_sample["n_misassemblies_overlapping_hotspot"] > 0).sum())
    n_samples_pr_overlap = int((per_sample["n_misassemblies_overlapping_pr"] > 0).sum())
    samples_hotspot_overlap = per_sample.loc[
        per_sample["n_misassemblies_overlapping_hotspot"] > 0, "SampleID"].tolist()
    samples_pr_overlap = per_sample.loc[
        per_sample["n_misassemblies_overlapping_pr"] > 0, "SampleID"].tolist()

    lines = []
    lines.append(f"{cohort} NucFlag misassembly analysis -- dataset summary")
    lines.append("=" * 60)
    lines.append(f"Total samples analyzed: {n_total}")
    lines.append("")
    lines.append(f"Q1. Total misassembly errors detected (all samples): {total_misassemblies}")
    lines.append("")
    lines.append(f"Q2. Samples with >=1 flagged assembly error: {n_with_any_error} / {n_total} "
                  f"({100 * n_with_any_error / n_total:.1f}%)")
    lines.append(f"    Sample list: {', '.join(samples_with_error)}")
    lines.append("")
    lines.append(f"Of {total_misassemblies} total misassemblies, {total_lifted} successfully "
                  f"lifted over to H37Rv coordinates ({total_not_lifted} did not and are "
                  f"excluded from Q3/Q4 below).")
    lines.append("")
    lines.append(f"Q3. Misassemblies (of the {total_lifted} lifted) overlapping a NucDivHotspots "
                  f"window: {total_overlap_hotspot}")
    lines.append(f"    Samples with >=1 hotspot-overlapping misassembly: "
                 f"{n_samples_hotspot_overlap} / {n_total} -- {', '.join(samples_hotspot_overlap)}")
    lines.append("")
    lines.append(f"Q4. Misassemblies (of the {total_lifted} lifted) overlapping a Paralogous "
                  f"Region (PR) window: {total_overlap_pr}")
    lines.append(f"    Samples with >=1 PR-overlapping misassembly: "
                 f"{n_samples_pr_overlap} / {n_total} -- {', '.join(samples_pr_overlap)}")
    lines.append("")
    lines.append("Caveat (deliberate simplification, first pass): misassemblies that did not "
                 "lift over to H37Rv coordinates are excluded entirely from Q3/Q4 (not counted "
                 "as either overlapping or non-overlapping).")

    summary_text = "\n".join(lines)
    (OUT_DIR / f"{cohort}_nucflag_dataset_summary.txt").write_text(summary_text + "\n")
    print(summary_text)
    print()


def main():
    for cohort, results_dir in DATASETS.items():
        analyze_dataset(cohort, results_dir)


if __name__ == "__main__":
    main()
