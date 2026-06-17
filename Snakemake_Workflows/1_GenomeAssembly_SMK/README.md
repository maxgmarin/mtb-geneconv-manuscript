# Genome Assembly Snakemake Pipelines

## Overview

These pipelines were used to generate complete *Mtb* genome assemblies for two datasets used in this study. The **Mtb151CI** dataset (151 complete genomes) is the primary assembly dataset, matching exactly the assemblies from [Marin et al. (2025)](https://academic.oup.com/bioinformatics/article/41/5/btaf219/8127202) — *Pitfalls of bacterial pan-genome analysis approaches: a case study of Mycobacterium tuberculosis and two less clonal bacterial species* ([github.com/farhat-lab/mtb-pg-benchmarking-2024paper](https://github.com/farhat-lab/mtb-pg-benchmarking-2024paper)). The **TBP-22-LR** dataset (22 genomes) was generated specifically for this study: PacBio HiFi long-read sequencing was performed on selected isolates to validate candidate gene conversion events initially identified from short-read WGS data alone.

All genomes were assembled using a hybrid approach combining long-read and short-read WGS data for each isolate. The general strategy was: (1) *de novo* assembly and long-read polishing with [Flye](https://github.com/fenderglass/Flye) (v2.6 or v2.9), followed by (2) short-read polishing with [Pilon](https://github.com/broadinstitute/pilon) (v1.23) using Illumina reads aligned to the draft assembly. For Oxford Nanopore assemblies, the higher base-call error rate of ONT v9.4.1 chemistry required additional polishing steps: a long-read polishing pass with [Medaka](https://github.com/nanoporetech/medaka) was performed before Pilon, and [PolyPolish](https://github.com/rrwick/Polypolish) was applied as a further short-read polishing step. All final assemblies were annotated with [Bakta](https://github.com/oschwengers/bakta) (v4.8).

Three separate pipelines were used to accommodate differences in long-read platform and chemistry:

| Pipeline | Long-read data | Key tools |
|---|---|---|
| `1.Mtb.Generate.HybridAsm.PacBioRSII.smk` | PacBio Subreads (RS II & Sequel II) | Flye v2.6, Pilon v1.23 |
| `2.Mtb.Generate.HybridAsm.PacBioHiFi.smk` | PacBio CCS/HiFi (Sequel II) | Flye v2.9, Pilon v1.23 |
| `3.Mtb.Generate.HybridAsm.ONT9.4.1.smk` | Oxford Nanopore v9.4.1 | Flye v2.6, Medaka v1.5.0, Pilon v1.23, PolyPolish v0.5.0 |

All pipelines were orchestrated using the [Snakemake](https://snakemake.readthedocs.io/) workflow system and run on the HMS O2 cluster (Slurm job scheduler). Each Snakemake file contains the exact commands and parameters used for each step.

## Pipelines

1. `1.Mtb.Generate.HybridAsm.PacBioRSII.smk` — Hybrid assembly of *Mtb* isolates sequenced with PacBio (Subreads, RS II & Sequel II) and Illumina WGS
2. `2.Mtb.Generate.HybridAsm.PacBioHiFi.smk` — Hybrid assembly of *Mtb* isolates sequenced with PacBio (CCS/HiFi, Sequel II) and Illumina WGS
3. `3.Mtb.Generate.HybridAsm.ONT9.4.1.smk` — Hybrid assembly of *Mtb* isolates sequenced with Oxford Nanopore (chemistry 9.4.1) and Illumina WGS

## Inputs

Each pipeline takes three key inputs:

| Input | Description |
|---|---|
| `output_dir` | Target directory where all pipeline results will be written |
| `--configfile` | JSON config file (`Mtb-WGA.config.json`) specifying tool parameters and settings |
| `inputSampleData_TSV` | Sample sheet TSV specifying per-isolate input data (see below) |

### Sample sheet format

The sample sheet is a TSV file with one row per isolate containing:
- `SampleID` — unique isolate identifier
- Path to paired-end short-read FASTQs (Illumina WGS)
- Path to long-read raw reads in FASTQ format (PacBio or ONT)

### Input manifests used in this study

The following TSVs in `Input_WGS_FQ_PATHs/` are the exact sample sheets each pipeline was run on. Each file specifies the SampleID and input FASTQ paths for all isolates in that dataset.

- [Mtb151.SetAandB.WiPacBio_Subreads.48Samples.LRandSR.InputWGS.PATHs.tsv](Input_WGS_FQ_PATHs/Mtb151.SetAandB.WiPacBio_Subreads.48Samples.LRandSR.InputWGS.PATHs.tsv)
- [Mtb151.SetC.WiPacBio_HiFi_CCS.8Samples.LRandSR.InputWGS.PATHs.tsv](Input_WGS_FQ_PATHs/Mtb151.SetC.WiPacBio_HiFi_CCS.8Samples.LRandSR.InputWGS.PATHs.tsv)
- [Mtb151.SetD.Hall2022.WiONT94_Reads.78Samples.LRandSR.InputWGS.PATHs.tsv](Input_WGS_FQ_PATHs/Mtb151.SetD.Hall2022.WiONT94_Reads.78Samples.LRandSR.InputWGS.PATHs.tsv)
- [Mtb151.SetE.Peker2021.WiONT94_Reads.18Samples.LRandSR.InputWGS.PATHs.tsv](Input_WGS_FQ_PATHs/Mtb151.SetE.Peker2021.WiONT94_Reads.18Samples.LRandSR.InputWGS.PATHs.tsv)
- [TBP-22-LR.TGEN_Reseq_Isolates.WiPacBio_HiFi_CCS.22Samples.InputWGS_PATHs.tsv](Input_WGS_FQ_PATHs/TBP-22-LR.TGEN_Reseq_Isolates.WiPacBio_HiFi_CCS.22Samples.InputWGS_PATHs.tsv)

## Usage

```bash
# Define key paths
targetOutput_Dir="/path/to/output/directory"
inputConfigFile="/path/to/Mtb-WGA.config.json"
inputSampleData_TSV="/path/to/sample_sheet.tsv"

mkdir -p ${targetOutput_Dir}

# Run pipeline (example using PacBio HiFi)
snakemake -s 2.Mtb.Generate.HybridAsm.PacBioHiFi.smk \
    --config output_dir=${targetOutput_Dir} \
             inputSampleData_TSV=${inputSampleData_TSV} \
    --configfile ${inputConfigFile} \
    -np
```

> Remove the `-np` (dry-run) flag to execute the pipeline.

## Slurm Cluster Submission (HMS O2)

All pipelines were run on the HMS O2 cluster using Slurm. The `SlurmClusterConfigs/` subdirectory contains per-pipeline JSON files that define the resource requests (partition, CPUs, memory, walltime) for each Snakemake rule:

| Config file | Pipeline |
|---|---|
| `clusterConfig.HybridAsm.PacBioRSII.json` | PacBio Subreads (RS II) hybrid assembly |
| `clusterConfig.HybridAsm.PacBioCCS.json` | PacBio CCS/HiFi hybrid assembly |
| `clusterConfig.HybridAsm.ONTAssembly.json` | Oxford Nanopore hybrid assembly |

These are passed to Snakemake via `--cluster-config`, which maps rule-level resource keys (e.g. `{cluster.mem}`) into the `sbatch` submission command. Example for the PacBio CCS pipeline:

```bash
targetOutput_Dir="/path/to/output/directory"
inputConfigFile="/path/to/Mtb-WGA.config.json"
inputSampleData_TSV="/path/to/sample_sheet.tsv"

mkdir -p ${targetOutput_Dir}/O2logs/cluster/

snakemake -s 2.Mtb.Generate.HybridAsm.PacBioHiFi.smk \
    --config output_dir=${targetOutput_Dir} \
             inputSampleData_TSV=${inputSampleData_TSV} \
    --configfile ${inputConfigFile} \
    -p --use-conda -j 250 \
    --cluster-config SlurmClusterConfigs/clusterConfig.HybridAsm.PacBioCCS.json \
    --cluster "sbatch -p {cluster.p} -n {cluster.n} -t {cluster.t} --mem {cluster.mem} \
               -o ${targetOutput_Dir}/O2logs/cluster/{cluster.o} \
               -e ${targetOutput_Dir}/O2logs/cluster/{cluster.e}" \
    --latency-wait 35 -k
```

Key flags:
- `-j 250` — allow up to 250 jobs to be submitted concurrently
- `--use-conda` — activate per-rule Conda environments as specified in the workflow
- `--latency-wait 35` — wait up to 35 seconds for output files to appear on the shared filesystem after a job completes
- `-k` — keep going with independent rules if one job fails

## References

The `references/` subdirectory contains reference files used by the assembly pipelines:
- Trimmomatic adapter lists for Illumina WGS trimming
- *M. tuberculosis* H37Rv *dnaA* gene/sequence (used for genome reorientation)
