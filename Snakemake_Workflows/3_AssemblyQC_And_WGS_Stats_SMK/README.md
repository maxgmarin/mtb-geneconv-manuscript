# 3_AssemblyQC_And_WGS_Stats_SMK

## Overview

Assembly QC pipeline for the hybrid long-read + short-read genome assemblies used
throughout this manuscript. 

The key steps that are run for each sample are: <br> 
- align its long-read WGS data to its own
assembly
- run [NucFlag](https://github.com/logsdon-lab/NucFlag) to detect
misassemblies
- lift over and annotate the results against H37Rv coordinates and two curated
region sets (NucDivHotspots windows, Paralogous Regions). This is done so misassemblies can be
checked for overlap with regions of interest identified in the manuscript.

The snakemake pipeline is defined by the **`AsmAndWGS.QC.smk`** file. 

## Inputs

| Input | Description |
|---|---|
| Sample sheet (`samples_tsv`) | TSV with columns `SampleID`, `LR_FQ_PATH`, `LR_Technology` (`PacBio_HiFi`, `ONT_r941`, or `PacBio_Subreads`), `AssemblyPath`. See `config/samples_example.tsv`. |
| H37Rv reference FASTA (`h37rv_fasta`) | RefSeq ASM19595v2, contig `NC_000962.3`. |
| Paralogous Regions table (`h37rv_paralogous_regions_tsv`) | `../../Data/H37Rv.NonUniqueSeqRegions.V1/H37Rv.Minimap2.HomologyMap.ParalogousRegions.tsv.gz`. |
| NucDivHotspots table | `resources/NucDivHotspots.37.1kb.tsv` (bundled with this pipeline, not per-run config). |
| NucFlag tool configs | `config/nucflag_final_hifi.toml`, `config/nucflag_final_ont_r9.toml` (one per long-read technology). |

## Pipeline Stages

See `docs/pipeline_stages.md` for the full breakdown of each stage:

1. Align long reads to the sample's own assembly (minimap2, primary-only filter, QC)
2. Run NucFlag (contig-boundary mask, `nucflag call`, `nucflag qv`)
3. Assembly<->H37Rv alignment (both orientations, feeding both liftover directions)
4. Liftover: misassembly calls -> H37Rv; NucDivHotspots windows and Paralogous Regions -> each assembly
5. Overlap misassembly calls against both region sets, then annotate with H37Rv coordinates + metadata

## Outputs

Per sample, under `{results_dir}/AsmAnalysis/{sample}/`:

| Output | Path |
|---|---|
| Misassembly calls + QV | `nucflag_output/{sample}.misassemblies.bed`, `{sample}.qv.bed` |
| Misassemblies lifted to H37Rv | `nucflag_output/NucFlag_Misassemblies_LiftoverToH37Rv/{sample}.misassemblies_liftover_to_H37Rv.bed` |
| NucDivHotspots overlap (annotated) | `nucflag_output/NucFlagStatus_Eval_37NucDivHotspots/{sample}.hotspot_status.{length,count}.annotated.bed` |
| Paralogous Regions overlap (annotated) | `nucflag_output/NucFlagStatus_Eval_ParalogousRegions/{sample}.pr_status.{length,count}.annotated.bed` |
| Primary-only LR-to-assembly alignment + QC | `alignment/{sample}.LR.AlnToOwnAssembly.primary_only.bam(.bai)`, `*.flagstat.txt`, `*.coverage.txt` |

## Key Tools

| Tool | Version | Purpose |
|---|---|---|
| minimap2 | 2.28 | Long-read-to-assembly and assembly-to-H37Rv alignment |
| samtools | 1.20 | BAM sort/index/filter/QC |
| nucflag | 1.0.0 | Misassembly calling, QV, region-overlap status |
| snakemake | 9.23.1 | Workflow orchestration |

Exact versions are pinned in `envs/environment.yml`.

## Environment setup

```bash
conda env create -f envs/environment.yml
conda activate nucflag_env
```

`envs/environment.yml` is an exact copy of the environment used to develop and run this
pipeline (also used by the downstream scripts in `scripts/AlnStatsQC_DepthAndBaseQ/` and
`scripts/NucFlag_ResultsSummary/`, which is why it also includes `bedtools`/`mosdepth`
even though the core pipeline rules here don't need them directly).

## Usage

### Dry run (check pipeline)

```bash
cd 3_AssemblyQC_And_WGS_Stats_SMK
snakemake -s AsmAndWGS.QC.smk -n
```

### Local execution

```bash
snakemake -s AsmAndWGS.QC.smk --cores 4 \
  --config samples_tsv=config/samples_example.tsv results_dir=/path/to/results
```

### HMS O2 Slurm cluster submission

```bash
snakemake -s AsmAndWGS.QC.smk --profile profiles/o2_slurm_generic --jobs 50 \
  --config samples_tsv=config/samples_example.tsv results_dir=/path/to/results
```

See `docs/run_commands.md` for more detail, plus the real invocation used to generate
the results summarized below.

## Utility Python scripts used by the snakemake pipeline
The five `scripts/*.py` files at the top level of `scripts/` (`build_boundary_mask.py`,
`merge_liftover_with_metadata.py`, `annotate_hotspot_status_with_nucdiv.py`,
`annotate_pr_status_with_hmregion.py`, `normalize_pr_regions_bed.py`) are `script:`
targets used directly by `AsmAndWGS.QC.smk` itself (part of the core pipeline, not a
downstream analysis).

## Downstream analysis scripts (`scripts/`)

There are two sub-directories that describe the downstream analysis scripts that consume
this pipeline's outputs.

- **`scripts/AlnStatsQC_DepthAndBaseQ/`** -- long-read alignment depth and base-quality
  comparison (Paralogous Regions vs. genome-wide).
- **`scripts/NucFlag_ResultsSummary/`** -- misassembly summary statistics across both
  cohorts, plus the already-generated result tables (`results/`) those scripts produced.


## Input Manifests

See `config/samples_example.tsv` for the sample-sheet format. Real manifests used for
this manuscript's cohorts are tracked elsewhere in this repo -- see
`../../Data/` and `../1_GenomeAssembly_SMK/Input_WGS_FQ_PATHs/`.
