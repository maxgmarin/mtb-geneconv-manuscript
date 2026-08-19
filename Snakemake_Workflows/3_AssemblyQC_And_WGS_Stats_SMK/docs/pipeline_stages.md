# Pipeline stages

`AsmAndWGS.QC.smk` runs 5 stages per sample. Each stage lists the rules involved, the
key tools/parameters, and its outputs.

## Stage 1: Align long reads to the sample's own assembly

| Rules | Key tools & parameters | Outputs |
|---|---|---|
| `stage_assembly`, `faidx_assembly` | `samtools faidx` | `assembly/{sample}.fasta(.fai)` |
| `minimap2_align_sort` | `minimap2 -ax <preset> --MD --eqx` (preset by `LR_Technology`: `lr:hqae` HiFi, `map-ont` ONT/Subreads) + `samtools sort` | `alignment/{sample}.LR.AlnToOwnAssembly.bam(.bai)` |
| `filter_and_index_primary_only` | `samtools view -F 2308` (drop unmapped/secondary/supplementary) | `alignment/{sample}.LR.AlnToOwnAssembly.primary_only.bam(.bai)` |
| `qc_full_bam`, `qc_filtered_bam` | `samtools flagstat` / `samtools coverage` | `*.flagstat.txt`, `*.coverage.txt` (full and primary-only) |

`--MD --eqx` is required for NucFlag's per-base mismatch/identity pileup tracks to
populate; it does not affect misassembly calls themselves.

## Stage 2: Run NucFlag

| Rules | Key tools & parameters | Outputs |
|---|---|---|
| `build_boundary_mask` | fixed 3bp mask at each contig end (circular-genome linearization artifact) | `nucflag_output/ignore_boundary.bed` |
| `nucflag_call` | `nucflag call --overlap_calls --ignore_regions <boundary mask>` with a technology-specific config (`config/nucflag_final_hifi.toml` or `nucflag_final_ont_r9.toml`) | `{sample}.misassemblies.bed`, `{sample}.status.bed`, `{sample}_plots/` |
| `nucflag_qv` | `nucflag qv` | `{sample}.qv.bed` |

## Stage 3: Assembly<->H37Rv alignment

| Rules | Key tools & parameters | Outputs |
|---|---|---|
| `mm2_asm_to_h37rv` | `minimap2 -ax asm20 --MD --cs` (target=H37Rv, query=Assembly) + mapq==60 filter | `MM2_AsmToH37Rv_asm20/{sample}.mm2.AsmToH37Rv.asm20.{bam,paf}` |
| `mm2_h37rv_to_asm_liftover` | `minimap2 -cx asm20 --cs` (target=Assembly, query=H37Rv -- the OPPOSITE orientation) | `MM2_H37RvToAsm_asm20_ForLiftOver/{sample}.mm2.H37RvToAsm.asm20.paf` |

Two alignments are needed because `paftools.js liftover <PAF> <BED>` always lifts from
the PAF's query to its target, and coordinates need to flow in both directions
(misassembly calls: assembly -> H37Rv; region-of-interest BEDs: H37Rv -> assembly).

## Stage 4: Liftover (both directions)

| Rules | Direction | Outputs |
|---|---|---|
| `misassemblies_only_bed`, `liftover_misassemblies_raw`, `liftover_misassemblies_to_h37rv` | Misassembly calls: assembly -> H37Rv | `NucFlag_Misassemblies_LiftoverToH37Rv/{sample}.misassemblies_liftover_to_H37Rv.bed` |
| `make_nucdiv_hotspots_bed`, `liftover_nucdiv_hotspots_raw`, `liftover_nucdiv_hotspots_to_asm` | NucDivHotspots windows: H37Rv -> assembly | `37NucDivHotspots_LiftoverToAsm/{sample}.NucDivHotspots.liftover_to_Asm.bed` |
| `normalize_pr_regions_bed`, `liftover_pr_regions_raw`, `liftover_pr_regions_to_asm` | Paralogous Regions: H37Rv -> assembly | `ParalogousRegions_LiftoverToAsm/{sample}.ParalogousRegions.liftover_to_Asm.bed` |

Each direction follows the same 3-step pattern: normalize the source table to a
header-free BED, run `paftools.js liftover` (which does not preserve extra metadata
columns), then rejoin the lifted coordinates with the original row's metadata
(`scripts/merge_liftover_with_metadata.py`) -- anchored on the *full* region list rather
than the liftover output, since `paftools.js liftover` silently drops rows that fail to
lift (see the `Lifted` column in later outputs).

## Stage 5: Overlap misassemblies against both region sets, then annotate

| Rules | Key tools & parameters | Outputs |
|---|---|---|
| `nucflag_status_hotspots`, `nucflag_status_pr_regions` | `nucflag status -g region -m {length,count}` | `*_status.{length,count}.bed` |
| `annotate_hotspot_status` (`scripts/annotate_hotspot_status_with_nucdiv.py`) | rejoin H37Rv coords + NucDiv metadata, anchored on all 37 windows | `*_status.{length,count}.annotated.bed` |
| `annotate_pr_status` (`scripts/annotate_pr_status_with_hmregion.py`) | rejoin H37Rv coords + PR metadata (keyed on `HmRegionID`), anchored on the full PR list | `*_status.{length,count}.annotated.bed` |

`length` mode reports the percentage of each region covered by each call type plus a QV
verdict; `count` mode reports the literal number of call segments intersecting each
region -- both are cheap to generate together and answer slightly different questions.
