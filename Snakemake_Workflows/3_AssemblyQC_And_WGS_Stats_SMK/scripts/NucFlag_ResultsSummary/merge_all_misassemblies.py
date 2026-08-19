"""
Builds one big, call-level TSV per dataset (Mtb151, TBP22CI) with every misassembly
call NucFlag flagged across every sample -- for downstream deep-dive analysis (not
just the headline counts produced by analyze_nucflag_results.py).

One row per misassembly call. Columns:
  SampleID, Cohort          -- which sample/dataset this call is from
  chrom, chromStart, chromEnd, name, score, zscore, af
                             -- the call itself, in the sample's OWN assembly coordinates
                                (chrom here is just the sample/contig name NucFlag used,
                                not H37Rv)
  Lifted                    -- bool, did this call successfully lift over to H37Rv coords?
  H37Rv_chrom, H37Rv_start, H37Rv_end, ClampFlags
                             -- H37Rv coordinates if Lifted else NaN
  overlaps_hotspot, overlaps_pr
                             -- bool if Lifted (tested against the 37 NucDivHotspots /
                                ~200 Paralogous Region windows), else NA (not evaluable)

Includes ALL misassemblies, lifted or not -- unlike analyze_nucflag_results.py's Q3/Q4,
which only counts lifted ones. Filter on `Lifted` yourself if you want the
lifted-only subset that matches those headline numbers exactly.

Run with the nucflag_env Python (has pandas):
  /n/data1/hms/dbmi/farhat/mm774/Projects/Mtb-GeneConv/NucFlag_Testing/envs/nucflag_env/bin/python \
    merge_all_misassemblies.py
"""

import warnings
from pathlib import Path
import pandas as pd

# Benign: some samples' overlaps_hotspot/overlaps_pr columns are all-NA (no lifted
# calls), which triggers a pandas future-dtype-inference warning on concat. Verified
# this doesn't affect the actual values written out (cross-checked against
# analyze_nucflag_results.py's independently-computed totals and bedtools -- all match).
warnings.filterwarnings("ignore", category=FutureWarning)

from nucflag_common import (
    DATASETS, discover_samples, load_misassemblies, load_liftover,
    load_region_bed, intervals_overlap,
)

OUT_DIR = Path(__file__).parent.parent / "results"
OUT_DIR.mkdir(exist_ok=True)


OUTPUT_COLS = ["SampleID", "Cohort", "chrom", "chromStart", "chromEnd", "name", "score",
               "zscore", "af", "Lifted", "H37Rv_chrom", "H37Rv_start", "H37Rv_end",
               "ClampFlags", "overlaps_hotspot", "overlaps_pr"]


def calls_for_sample(sample, cohort, results_dir, hotspot_bed, pr_bed):
    mis = load_misassemblies(sample, results_dir)
    flagged = mis[mis["name"] != "correct"].copy()
    if flagged.empty:
        return pd.DataFrame(columns=OUTPUT_COLS)

    lift = load_liftover(sample, results_dir).rename(
        columns={"asm_start": "chromStart", "asm_end": "chromEnd", "clamp_flags": "ClampFlags"}
    )[["chromStart", "chromEnd", "name", "H37Rv_chrom", "H37Rv_start", "H37Rv_end", "ClampFlags"]]

    n_before = len(flagged)
    merged = flagged.merge(lift, on=["chromStart", "chromEnd", "name"], how="left")
    assert len(merged) == n_before, (
        f"{sample}: merge changed row count ({n_before} -> {len(merged)}) -- "
        "join key (chromStart, chromEnd, name) may not be unique for this sample"
    )

    merged["Lifted"] = merged["H37Rv_chrom"].notna()
    merged["overlaps_hotspot"] = pd.array([pd.NA] * len(merged), dtype="boolean")
    merged["overlaps_pr"] = pd.array([pd.NA] * len(merged), dtype="boolean")

    lifted_mask = merged["Lifted"]
    if lifted_mask.any():
        lifted_rows = merged.loc[lifted_mask]
        merged.loc[lifted_mask, "overlaps_hotspot"] = lifted_rows.apply(
            lambda r: bool(intervals_overlap(r.H37Rv_start, r.H37Rv_end,
                                              hotspot_bed.Start, hotspot_bed.End).any()),
            axis=1,
        ).astype("boolean")
        merged.loc[lifted_mask, "overlaps_pr"] = lifted_rows.apply(
            lambda r: bool(intervals_overlap(r.H37Rv_start, r.H37Rv_end,
                                              pr_bed.Start, pr_bed.End).any()),
            axis=1,
        ).astype("boolean")

    merged["SampleID"] = sample
    merged["Cohort"] = cohort
    return merged


def merge_dataset(cohort, results_dir):
    samples = discover_samples(results_dir)
    hotspot_bed = load_region_bed("hotspot", results_dir)
    pr_bed = load_region_bed("pr", results_dir)

    per_sample_calls = [calls_for_sample(s, cohort, results_dir, hotspot_bed, pr_bed) for s in samples]
    non_empty = [df for df in per_sample_calls if not df.empty]
    all_calls = pd.concat(non_empty, ignore_index=True) if non_empty else pd.DataFrame(columns=OUTPUT_COLS)

    all_calls = all_calls[OUTPUT_COLS]

    out_path = OUT_DIR / f"{cohort}_all_misassemblies_merged.tsv"
    all_calls.to_csv(out_path, sep="\t", index=False)
    print(f"{cohort}: {len(all_calls)} total misassembly calls across {len(samples)} samples "
          f"-> {out_path}")
    print(f"  Lifted: {int(all_calls['Lifted'].sum())} / {len(all_calls)}")
    print(f"  overlaps_hotspot=True: {int((all_calls['overlaps_hotspot'] == True).sum())}")
    print(f"  overlaps_pr=True: {int((all_calls['overlaps_pr'] == True).sum())}")
    return all_calls


def main():
    for cohort, results_dir in DATASETS.items():
        merge_dataset(cohort, results_dir)


if __name__ == "__main__":
    main()
