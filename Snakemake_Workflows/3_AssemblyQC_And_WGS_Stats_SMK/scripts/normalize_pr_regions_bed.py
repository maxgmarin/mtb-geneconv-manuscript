# Converts the user's Paralogous Regions TSV (H37Rv coordinates; header + Chr/Start/End
# among many other columns, plus a duplicate lowercase chrom/start/end at the end) into
# a header-free BED-with-metadata file: chrom, start, end, then metadata columns.
#
# Unlike NucDivHotspots' source TSV (already Chrom/Start/End as the first 3 columns --
# make_nucdiv_hotspots_bed just does `awk 'NR>1'`), the PR TSV's Chr/Start/End are
# columns 2-4, not 1-3, so this needs a real column reorder for paftools.js liftover to
# parse the output as a BED at all. Kept metadata columns: HmRegionID (unique ID),
# Overlap_Genes, Length, PR_SetID -- matches the extra_cols already used for
# samtools_coverage_lr_paralogous_regions (AlnToH37Rv.AsmAndReads.smk).
#
# Validated ad hoc in NucFlag_Testing/ToolTesting/ParalogousRegionsLiftoverTest/ before
# writing this rule.

import pandas as pd

df = pd.read_csv(snakemake.input.tsv, sep="\t")
out = df[["Chr", "Start", "End", "HmRegionID", "Overlap_Genes", "Length", "PR_SetID"]]
out.to_csv(snakemake.output.bed, sep="\t", header=False, index=False)
