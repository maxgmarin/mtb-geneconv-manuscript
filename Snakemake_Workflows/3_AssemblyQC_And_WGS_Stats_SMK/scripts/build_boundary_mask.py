# Snakemake `script:` target for rule build_boundary_mask.
# SKETCH ONLY -- not yet run.
#
# Reads the contig length from the assembly's .fai and writes a 2-line BED masking
# the first/last `params.mask_bp` bases -- the circular-genome linearization
# artifact. Assumes a single-contig (complete) assembly, matching our current
# Mtb hybrid assemblies; would need to loop over contigs for multi-contig input.

with open(snakemake.input.fai) as fh:
    fields = fh.readline().strip().split("\t")
    contig, contig_len = fields[0], int(fields[1])

mask_bp = snakemake.params.mask_bp

with open(snakemake.output.bed, "w") as out:
    out.write(f"{contig}\t0\t{mask_bp}\n")
    out.write(f"{contig}\t{contig_len - mask_bp}\t{contig_len}\n")
