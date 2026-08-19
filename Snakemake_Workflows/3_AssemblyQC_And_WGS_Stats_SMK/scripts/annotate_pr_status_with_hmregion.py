# Joins nucflag_status_pr_regions' per-region output (assembly coordinates only, no
# identity info) back against the full Paralogous Regions list (H37Rv coordinates +
# metadata), by way of the liftover-to-Asm bed -- same pattern as
# annotate_hotspot_status_with_nucdiv.py, anchored on the FULL region list so a PR that
# fails to lift over still gets a row (Asm coordinates and status columns NaN for that
# row, Lifted=False) instead of silently vanishing.
#
# Joined on HmRegionID (a real unique ID already present in the source data) rather
# than a (chrom, start, end) coordinate triple like the NucDivHotspots version --
# NucDivHotspots had no dedicated ID column so its annotate script had to key on
# coordinates out of necessity; PRs don't have that constraint.
#
# Validated ad hoc in NucFlag_Testing/ToolTesting/ParalogousRegionsLiftoverTest/ before
# writing this rule -- see that directory's README.md for a caveat observed there:
# paftools.js liftover can degenerately collapse two different H37Rv source regions
# onto the identical 1bp assembly coordinate when the true target isn't confidently
# alignable (with no clamp-flag set on either). The merge below handles the resulting
# fan-out correctly regardless -- both regions keep their own row, each correctly
# annotated with the (identical) status for that shared coordinate.

import logging

import pandas as pd

logging.basicConfig(
    filename=snakemake.log[0],
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(message)s",
)

FULL_PR_COLS = ["H37Rv_chrom", "H37Rv_start", "H37Rv_end", "HmRegionID", "Overlap_Genes", "Length", "PR_SetID"]
LIFTOVER_COLS = [
    "Asm_chrom", "Asm_start", "Asm_end",
    "H37Rv_chrom", "H37Rv_start", "H37Rv_end",
    "HmRegionID", "Overlap_Genes", "Length", "PR_SetID", "ClampFlags",
]

full_prs = pd.read_csv(snakemake.input.full_prs, sep="\t", names=FULL_PR_COLS, header=None)
liftover = pd.read_csv(snakemake.input.liftover, sep="\t", names=LIFTOVER_COLS, header=None)

base = full_prs.merge(
    liftover[["HmRegionID", "Asm_chrom", "Asm_start", "Asm_end", "ClampFlags"]],
    on="HmRegionID",
    how="left",
)
base["Asm_start"] = base["Asm_start"].astype("Int64")
base["Asm_end"] = base["Asm_end"].astype("Int64")
base["Lifted"] = base["Asm_chrom"].notna()

n_total = len(base)
n_not_lifted = int((~base["Lifted"]).sum())
if n_not_lifted:
    logging.warning(
        f"{n_not_lifted} of {n_total} Paralogous Regions did not lift over for this "
        f"sample's assembly -- Asm coordinates and all status columns are NaN for "
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
            f"{n_lifted_no_status} of {n_total} regions lifted over but had no matching "
            f"row in {status_path} (unexpected -- nucflag status should cover every "
            f"lifted region)."
        )

    merged.to_csv(output_path, sep="\t", index=False)
    logging.info(
        f"Wrote {len(merged)} rows to {output_path} "
        f"({n_total - n_not_lifted} lifted, {n_not_lifted} not lifted)."
    )


annotate(snakemake.input.length, snakemake.output.length_annotated)
annotate(snakemake.input.count_bed, snakemake.output.count_annotated)
