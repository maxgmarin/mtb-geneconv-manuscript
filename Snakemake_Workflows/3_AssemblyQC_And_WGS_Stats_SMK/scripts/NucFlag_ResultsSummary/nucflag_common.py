"""
Shared config/loaders for the Mtb151 + TBP22CI NucFlag analysis scripts in this
directory. See docs/NucFlag_Pipeline_Overview.md for exact file-format semantics.
"""

from pathlib import Path
import pandas as pd

BASE = Path("/n/data1/hms/dbmi/farhat/mm774/Projects/Mtb-GeneConv")

DATASETS = {
    "Mtb151": BASE / "NucFlag_Snakemake_TestOutputs_Mtb151",
    "TBP22CI": BASE / "NucFlag_Snakemake_TestOutputs_TBP22CI",
}

MISASSEMBLY_COLS = ["chrom", "chromStart", "chromEnd", "name", "score", "strand",
                     "thickStart", "thickEnd", "itemRgb", "zscore", "af"]

LIFTOVER_COLS = ["H37Rv_chrom", "H37Rv_start", "H37Rv_end", "sample", "asm_start",
                  "asm_end", "name", "clamp_flags"]


def discover_samples(results_dir):
    asm_dir = results_dir / "AsmAnalysis"
    samples = []
    for d in sorted(asm_dir.iterdir()):
        if (d / "nucflag_output" / f"{d.name}.misassemblies.bed").exists():
            samples.append(d.name)
    return samples


def load_misassemblies(sample, results_dir):
    path = results_dir / "AsmAnalysis" / sample / "nucflag_output" / f"{sample}.misassemblies.bed"
    return pd.read_csv(path, sep="\t", names=MISASSEMBLY_COLS, header=0)


def load_liftover(sample, results_dir):
    path = (results_dir / "AsmAnalysis" / sample / "nucflag_output"
            / "NucFlag_Misassemblies_LiftoverToH37Rv"
            / f"{sample}.misassemblies_liftover_to_H37Rv.bed")
    if not path.exists() or path.stat().st_size == 0:
        return pd.DataFrame(columns=LIFTOVER_COLS)
    return pd.read_csv(path, sep="\t", names=LIFTOVER_COLS, header=None)


def load_region_bed(kind, results_dir):
    if kind == "hotspot":
        path = results_dir / "_shared" / "NucDivHotspots.37.1kb.bed"
        cols = ["Chrom", "Start", "End", "Middle", "OverlapGenes", "NucDiv_kb",
                "NucDiv_ModZ", "OvrlapWi_HmMapAln_PR", "OvrlapWi_HmMapAln_LR",
                "OvrlapWi_LowPmap", "PLC_Mask_Ovrlap"]
    elif kind == "pr":
        path = results_dir / "_shared" / "H37Rv.ParalogousRegions.bed"
        cols = ["Chrom", "Start", "End", "HmRegionID", "Overlap_Genes", "Length", "PR_SetID"]
    else:
        raise ValueError(kind)
    return pd.read_csv(path, sep="\t", names=cols, header=None)


def intervals_overlap(s1, e1, s2, e2):
    return (s1 < e2) & (e1 > s2)


def overlapping_call_mask(liftover_df, region_df):
    """Boolean mask, one per row of liftover_df: does this lifted call's H37Rv interval
    overlap >=1 window in region_df?"""
    if liftover_df.empty:
        return pd.Series([], dtype=bool)
    return liftover_df.apply(
        lambda call: bool(intervals_overlap(call.H37Rv_start, call.H37Rv_end,
                                             region_df.Start, region_df.End).any()),
        axis=1,
    )
