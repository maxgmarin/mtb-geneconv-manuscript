# Joins nucflag_status_hotspots' per-window output (assembly coordinates) back against
# the full 37-window NucDivHotspots list (H37Rv coordinates + metadata), by way of
# NucDivHotspotsLiftover/{sample}.NucDivHotspots.liftover_to_Asm.bed.
#
# Anchored on the FULL window list (AsmAnalysis/_shared/NucDivHotspots.37.1kb.bed), NOT on
# the liftover file or the hotspot_status file -- paftools.js liftover silently drops rows
# that don't lift over (not just reorders them), so a plain join keyed off either of those
# would silently drop non-lifted windows from the output entirely, making "no misassembly
# found" and "couldn't even check this window" look identical. Instead every one of the 37
# windows always gets a row here: Asm_chrom/Asm_start/Asm_end/ClampFlags and every
# hotspot_status column (status, QV, per-type percentages/counts) are NaN when the window
# didn't lift over for this sample's assembly -- the `Lifted` column makes this explicit
# and filterable rather than something you have to infer from NaN presence.
#
# Safe as an exact-coordinate join once a window does have Asm coordinates:
# nucflag_status_hotspots' `-g region` grouping trims/merges calls to exactly the input
# region's [chromStart, chromEnd) bounds (confirmed from nucflag/call/status.py), so a
# lifted window's hotspot_status row always has the same coordinates as its liftover row.
#
# NOTE the output can have MORE than 37 rows for a sample: paftools.js liftover emits one
# row per alignment block an input region overlaps, so a window spanning a breakpoint
# between two blocks (observed on a PE_PGRS54 window, mada_2-46, Mtb151) produces TWO
# liftover rows for that one H37Rv window -- one "start_clamped", one "end_clamped",
# each landing in a different assembly region. Both are kept (not deduplicated/merged):
# they can carry genuinely different hotspot_status results, and collapsing them would
# hide a real repeat-region alignment fragmentation signal.

import logging

import pandas as pd

logging.basicConfig(
    filename=snakemake.log[0],
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(message)s",
)

FULL_HOTSPOT_COLS = [
    "H37Rv_chrom", "H37Rv_start", "H37Rv_end",
    "Middle", "OverlapGenes", "NucDiv_kb", "NucDiv_ModZ",
    "OvrlapWi_HmMapAln_PR", "OvrlapWi_HmMapAln_LR", "OvrlapWi_LowPmap",
    "PLC_Mask_Ovrlap",
]
LIFTOVER_COLS = [
    "Asm_chrom", "Asm_start", "Asm_end",
    "H37Rv_chrom", "H37Rv_start", "H37Rv_end",
    "Middle", "OverlapGenes", "NucDiv_kb", "NucDiv_ModZ",
    "OvrlapWi_HmMapAln_PR", "OvrlapWi_HmMapAln_LR", "OvrlapWi_LowPmap",
    "PLC_Mask_Ovrlap", "ClampFlags",
]
H37RV_KEY = ["H37Rv_chrom", "H37Rv_start", "H37Rv_end"]

full_hotspots = pd.read_csv(
    snakemake.input.full_hotspots, sep="\t", names=FULL_HOTSPOT_COLS, header=None
)
liftover = pd.read_csv(snakemake.input.liftover, sep="\t", names=LIFTOVER_COLS, header=None)

base = full_hotspots.merge(
    liftover[H37RV_KEY + ["Asm_chrom", "Asm_start", "Asm_end", "ClampFlags"]],
    on=H37RV_KEY,
    how="left",
)
base["Asm_start"] = base["Asm_start"].astype("Int64")
base["Asm_end"] = base["Asm_end"].astype("Int64")
base["Lifted"] = base["Asm_chrom"].notna()

n_total = len(base)
n_not_lifted = int((~base["Lifted"]).sum())
if n_not_lifted:
    logging.warning(
        f"{n_not_lifted} of {n_total} NucDivHotspots windows did not lift over for this "
        f"sample's assembly -- Asm coordinates and all hotspot_status columns are NaN for "
        f"those rows (see the 'Lifted' column)."
    )


def annotate(status_path, output_path):
    status = pd.read_csv(status_path, sep="\t").rename(
        columns={"#chrom": "Asm_chrom", "chromStart": "Asm_start", "chromEnd": "Asm_end"}
    )
    status["Asm_start"] = status["Asm_start"].astype("Int64")
    status["Asm_end"] = status["Asm_end"].astype("Int64")

    merged = base.merge(status, on=["Asm_chrom", "Asm_start", "Asm_end"], how="left")

    n_lifted_no_status = int((merged["Lifted"] & merged.drop(columns=base.columns).isna().all(axis=1)).sum())
    if n_lifted_no_status:
        logging.warning(
            f"{n_lifted_no_status} of {n_total} windows lifted over but had no matching "
            f"row in {status_path} (unexpected -- nucflag status should cover every "
            f"lifted window)."
        )

    merged.to_csv(output_path, sep="\t", index=False)
    logging.info(
        f"Wrote {len(merged)} rows to {output_path} "
        f"({n_total - n_not_lifted} lifted, {n_not_lifted} not lifted)."
    )


annotate(snakemake.input.length, snakemake.output.length_annotated)
annotate(snakemake.input.count_bed, snakemake.output.count_annotated)
