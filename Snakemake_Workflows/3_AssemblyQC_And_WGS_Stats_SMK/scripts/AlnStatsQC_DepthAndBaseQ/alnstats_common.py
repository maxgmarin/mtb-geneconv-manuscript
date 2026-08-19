"""
Shared config/loaders for the AlnStatsQC scripts. Reads only -- never touches the
NucFlag pipeline or its outputs. See docs/AlnStatsQC_Overview.md for exact file-format
semantics and known gaps.
"""

from pathlib import Path
import pandas as pd

RESULTS_DIR = Path("/n/data1/hms/dbmi/farhat/mm774/Projects/Mtb-GeneConv/NucFlag_Snakemake_TestOutputs_Mtb151")
ASM_DIR = RESULTS_DIR / "AsmAnalysis"

# Both-cohort dataset registry, added 2026-08-08 for the mosdepth-based genome/PR/non-PR
# depth summary (see compute_genome_pr_nonpr_depth_summary.py) -- same pattern as
# MtbGC_NucFlag_AsmAnalysis/scripts/nucflag_common.py's DATASETS dict. The pre-existing
# functions below this point (load_genome_coverage, load_pr_coverage, etc.) remain
# Mtb151-only via the module-level RESULTS_DIR/ASM_DIR globals -- not touched/extended
# here, since nothing in this task needs them generalized.
_BASE = Path("/n/data1/hms/dbmi/farhat/mm774/Projects/Mtb-GeneConv")
DATASETS = {
    "Mtb151": _BASE / "NucFlag_Snakemake_TestOutputs_Mtb151",
    "TBP22CI": _BASE / "NucFlag_Snakemake_TestOutputs_TBP22CI",
}

GENOME_COLS = ["rname", "startpos", "endpos", "numreads", "covbases", "coverage",
               "meandepth", "meanbaseq", "meanmapq"]


def discover_samples_in(results_dir):
    """Cohort-generic version of discover_samples() -- any sample with a mosdepth
    per-base file, for the new genome/PR/non-PR depth summary."""
    samples = []
    asm_dir = results_dir / "AsmAnalysis"
    for d in sorted(asm_dir.iterdir()):
        if (d / "mosdepth_coverage_stats" / f"{d.name}.per-base.bed.gz").exists():
            samples.append(d.name)
    return samples


def mosdepth_per_base_path(results_dir, sample):
    return results_dir / "AsmAnalysis" / sample / "mosdepth_coverage_stats" / f"{sample}.per-base.bed.gz"


def mosdepth_summary_path(results_dir, sample):
    return results_dir / "AsmAnalysis" / sample / "mosdepth_coverage_stats" / f"{sample}.mosdepth.summary.txt"


def load_pr_bed_for(results_dir):
    return results_dir / "_shared" / "H37Rv.ParalogousRegions.bed"


def load_nonpr_bed_for(results_dir):
    return results_dir / "_shared" / "H37Rv.NonParalogousRegions.bed"


def load_genome_coverage_for(results_dir, sample):
    """Cohort-generic version of load_genome_coverage() -- BAM-based, real meanbaseq."""
    path = results_dir / "AsmAnalysis" / sample / "samtools_coverage_stats" / f"{sample}.LR.AlnToH37Rv.coverage.genome.tsv"
    df = pd.read_csv(path, sep="\t", names=GENOME_COLS, header=0)
    assert len(df) == 1, f"{sample}: expected exactly 1 genome-wide row, got {len(df)}"
    return df.iloc[0]


def load_pr_coverage_for(results_dir, sample):
    """Cohort-generic version of load_pr_coverage() -- BAM-based, real meanbaseq."""
    path = results_dir / "AsmAnalysis" / sample / "samtools_coverage_stats" / f"{sample}.LR.AlnToH37Rv.coverage.paralogous_regions.tsv"
    return pd.read_csv(path, sep="\t")


def discover_samples():
    samples = []
    for d in sorted(ASM_DIR.iterdir()):
        if (d / "samtools_coverage_stats" / f"{d.name}.LR.AlnToH37Rv.coverage.genome.tsv").exists():
            samples.append(d.name)
    return samples


def load_genome_coverage(sample):
    path = ASM_DIR / sample / "samtools_coverage_stats" / f"{sample}.LR.AlnToH37Rv.coverage.genome.tsv"
    df = pd.read_csv(path, sep="\t", names=GENOME_COLS, header=0)
    assert len(df) == 1, f"{sample}: expected exactly 1 genome-wide row, got {len(df)}"
    return df.iloc[0]


def load_pr_coverage(sample):
    path = ASM_DIR / sample / "samtools_coverage_stats" / f"{sample}.LR.AlnToH37Rv.coverage.paralogous_regions.tsv"
    return pd.read_csv(path, sep="\t")


def load_1kb_windows_coverage(sample):
    """Returns None if this sample is missing its 1kb-windows file (known gap -- just
    TB3113 as of 2026-08-07, see docs/AlnStatsQC_Overview.md)."""
    path = ASM_DIR / sample / "samtools_coverage_stats" / f"{sample}.LR.AlnToH37Rv.coverage.1kb_windows.tsv"
    if not path.exists() or path.stat().st_size == 0:
        return None
    return pd.read_csv(path, sep="\t")


def load_inprsonly_coverage(sample):
    """Real (not arithmetic-derived) coverage stats from samtools coverage run on just
    the reads overlapping any Paralogous Region (CoverageStats_PRs_Vs_NonPRs.smk)."""
    path = ASM_DIR / sample / "CoverageStats_PRs_Vs_NonPRs" / f"{sample}.LR.AlnToH37Rv.coverage.InPRsOnly.tsv"
    df = pd.read_csv(path, sep="\t", names=GENOME_COLS, header=0)
    assert len(df) == 1, f"{sample}: expected exactly 1 row, got {len(df)}"
    return df.iloc[0]


def load_notinprs_coverage(sample):
    """Real (not arithmetic-derived) coverage stats from samtools coverage run on just
    the reads that do NOT overlap any Paralogous Region -- the exact read-level
    complement of load_inprsonly_coverage (guaranteed by construction, see
    CoverageStats_PRs_Vs_NonPRs.smk's header comment)."""
    path = ASM_DIR / sample / "CoverageStats_PRs_Vs_NonPRs" / f"{sample}.LR.AlnToH37Rv.coverage.NotInPRs.tsv"
    df = pd.read_csv(path, sep="\t", names=GENOME_COLS, header=0)
    assert len(df) == 1, f"{sample}: expected exactly 1 row, got {len(df)}"
    return df.iloc[0]


def load_nucdiv_hotspot_bed():
    path = RESULTS_DIR / "_shared" / "NucDivHotspots.37.1kb.bed"
    cols = ["Chrom", "Start", "End", "Middle", "OverlapGenes", "NucDiv_kb", "NucDiv_ModZ",
            "OvrlapWi_HmMapAln_PR", "OvrlapWi_HmMapAln_LR", "OvrlapWi_LowPmap", "PLC_Mask_Ovrlap"]
    return pd.read_csv(path, sep="\t", names=cols, header=None)


def bp_weighted_mean(values, weights):
    """Exact for additive per-base quantities like depth (mean_i * length_i summed,
    divided by summed length): weighted mean of `values` by `weights`."""
    total_weight = weights.sum()
    if total_weight == 0:
        return float("nan")
    return (values * weights).sum() / total_weight
