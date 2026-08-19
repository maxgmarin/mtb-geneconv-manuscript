# Reattaches original metadata columns to a paftools.js liftover output.
#
# paftools.js liftover's output BED does NOT preserve the input file's own column-4 value --
# it overwrites column 4 with an auto-generated "{orig_chrom}_{orig_start}_{orig_end}" label
# encoding the INPUT coordinates. This script uses that label to join back to the original
# (pre-liftover) file's full row, so the final output carries both the lifted and original
# coordinates plus whatever metadata columns the original file had (misassembly type, or
# NucDivHotspots gene/nucdiv-stat columns, etc. -- this script is agnostic to what they are).
#
# Robust to chrom/sample names that themselves contain underscores (e.g. "M0016395_7",
# "mada_2-31"): parses the encoded label with rsplit("_", 2) from the right, since start/end
# are always the last two purely-numeric underscore-delimited tokens.
#
# paftools.js also appends "_t5"/"_t3" suffixes (verified against its actual source,
# lh3/minimap2/misc/paftools.js's paf_liftover function) when a region's start (5') or end
# (3') coordinate falls outside any confidently-alignable block and gets CLAMPED to that
# block's boundary instead of being exactly projected -- i.e. a partial/approximate liftover,
# not (as originally assumed here) an "ambiguous multiple target hits" indicator. Both can
# appear together ("_t5_t3") if neither boundary lifts cleanly. Observed on a PE_PGRS repeat
# region in NucDivHotspots testing: a 1000bp H37Rv window lifted to only 446bp on the
# assembly, tagged "_t5" -- the lifted start is just the nearest alignment block's edge, not
# a real projected coordinate (consistent with PE_PGRS's known length/repeat variability
# between Mtb strains). Both flags are stripped before parsing chrom/start/end and recorded
# in their own output column so partially-lifted regions stay visible rather than being
# silently treated as exact.
#
# Also robust to rows that fail to lift over entirely (paftools.js silently drops them, not
# just reorders) -- logs how many of the original input rows made it through.

import logging
import re

logging.basicConfig(
    filename=snakemake.log[0],
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(message)s",
)

BOUNDARY_CLAMP_SUFFIX = re.compile(r"(_t5|_t3)+$")

# Load original (pre-liftover) rows, keyed by (chrom, start, end) -> full row's extra columns.
orig_by_coord = {}
with open(snakemake.input.orig_bed) as fh:
    for line in fh:
        fields = line.rstrip("\n").split("\t")
        chrom, start, end, *extra = fields
        orig_by_coord[(chrom, start, end)] = extra
n_total = len(orig_by_coord)

n_lifted = 0
n_clamped = 0
with (
    open(snakemake.input.liftover_raw) as fh_in,
    open(snakemake.output.bed, "w") as fh_out,
):
    for line in fh_in:
        fields = line.rstrip("\n").split("\t")
        lifted_chrom, lifted_start, lifted_end, encoded_name = fields[0:4]

        clamp_flags = []
        m = BOUNDARY_CLAMP_SUFFIX.search(encoded_name)
        if m:
            suffix = m.group(0)
            if "_t5" in suffix:
                clamp_flags.append("start_clamped")
            if "_t3" in suffix:
                clamp_flags.append("end_clamped")
            encoded_name = encoded_name[: m.start()]
            n_clamped += 1

        orig_chrom, orig_start, orig_end = encoded_name.rsplit("_", 2)
        extra = orig_by_coord.get((orig_chrom, orig_start, orig_end))
        if extra is None:
            logging.warning(
                f"Liftover output referenced unrecognized original coordinates "
                f"{orig_chrom}:{orig_start}-{orig_end} (encoded name {encoded_name!r}); skipping."
            )
            continue
        row = [
            lifted_chrom,
            lifted_start,
            lifted_end,
            orig_chrom,
            orig_start,
            orig_end,
            *extra,
            ";".join(clamp_flags),
        ]
        fh_out.write("\t".join(row) + "\n")
        n_lifted += 1

logging.info(
    f"{n_lifted} of {n_total} input regions successfully lifted over "
    f"({n_clamped} had a start and/or end coordinate clamped to the nearest alignment block "
    f"boundary rather than exactly projected -- see the boundary-clamp-flags column)."
)
