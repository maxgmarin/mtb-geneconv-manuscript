# `mtb-geneconv-manuscript`

Supporting code, data, and metadata for the manuscript "Gene conversion is a key driver of diversity hotspots in *Mycobacterium tuberculosis* antigens and virulence-associated loci"

**Maximillian G. Marin, Natalia Quinones-Olvera, Hu Jin, Michael A. Harris, Brendan M. Jeffrey, Alex Rosenthal, Kenan C. Murphy, Christopher Sassetti, Heng Li, Maha R. Farhat**

[![DOI:10.64898/2026.02.26.708061](https://img.shields.io/badge/bioRxiv-10.64898%2F2026.02.26.708061-b31b1b)](https://doi.org/10.64898/2026.02.26.708061)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18946124.svg)](https://doi.org/10.5281/zenodo.18946124)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

---

## Overview

Gene conversion — the non-reciprocal transfer of sequence information between homologous loci — has long been suspected to contribute to diversity in *Mycobacterium tuberculosis* (Mtb), yet its genome-wide extent and functional consequences have remained unclear. This study uses 151 complete hybrid Mtb genome assemblies to systematically characterize gene conversion events across the Mtb genome, revealing that gene conversion is a major driver of diversity hotspots in antigenic and virulence-associated loci including PPE18 and other PE/PPE family genes.

This repository contains all code, bioinformatics workflows, and analysis notebooks used to generate the results in the manuscript. 

---

## Table of Contents

- [Repository Structure](#repository-structure)
- [Data Availability](#data-availability)
- [Installation](#installation)
- [Snakemake Workflows](#snakemake-workflows)
- [Analysis Notebooks](#analysis-notebooks)
- [`gcutils` Python Package](#gcutils-python-package)
- [Citation](#citation)

---

## Repository Structure

| Directory | Description |
|---|---|
| [`0.GenomeData/`](0.GenomeData/) | Dataset descriptions, per-sample metadata, SRA accessions, and links to genome assemblies deposited on Zenodo |
| [`Snakemake_Workflows/`](Snakemake_Workflows/) | Snakemake pipelines for genome assembly and core whole-genome analysis (variant calling, phylogenetics, recombination detection) |
| [`Analysis/`](Analysis/) | Jupyter notebooks for all downstream analyses (nucleotide diversity, gene conversion events, epitope analysis, validation) |
| [`Data/`](Data/) | Input metadata tables and processed results files used by the analysis notebooks |
| [`References/`](References/) | H37Rv reference genome annotations, masking files, and other reference resources |
| [`gcutils/`](gcutils/) | Custom Python utility package installed as a dependency for the analysis notebooks |

---

## Data Availability

Complete genome assemblies for the **Mtb151CI** and **TBP-22-LR** datasets are deposited on Zenodo:

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18946124.svg)](https://doi.org/10.5281/zenodo.18946124)

Raw sequencing data accessions (SRA/ENA run accessions and BioProject IDs) for all three datasets are provided in **Supplemental Data Table 1 (SD1)**, documented in [`0.GenomeData/`](0.GenomeData/).

---


## Snakemake Workflows

Two Snakemake workflows were used to process raw sequencing data and generate the inputs for downstream analysis:

| Workflow | Description |
|---|---|
| [`Snakemake_Workflows/1_GenomeAssembly_SMK/`](Snakemake_Workflows/1_GenomeAssembly_SMK/) | Hybrid genome assembly from long-read + short-read WGS data (Flye + Pilon ± Medaka/PolyPolish) |
| [`Snakemake_Workflows/2_CoreAnalysis_SMK/`](Snakemake_Workflows/2_CoreAnalysis_SMK/) | Whole-genome analysis: reference alignment, variant calling, phylogenetics, nucleotide diversity, and recombination detection (Gubbins) |

See each workflow's `README.md` for inputs, usage, and exact run commands.

---

## Analysis Notebooks

All downstream analyses are implemented as Jupyter notebooks in `Analysis/`, organized by topic:

| Directory | Description |
|---|---|
| [`Analysis/1_Reference_Preprocessing/`](Analysis/1_Reference_Preprocessing/) | H37Rv reference annotation processing and homology/mappability masking |
| [`Analysis/2_Mtb_Prior_AsmQCandComparison/`](Analysis/2_Mtb_Prior_AsmQCandComparison/) | Assembly QC and WGS summary statistics |
| [`Analysis/3_Mtb_NucDivAnalysis_V5/`](Analysis/3_Mtb_NucDivAnalysis_V5/) | Nucleotide diversity analysis across the Mtb genome |
| [`Analysis/4_Mtb_RecombEvent_Explore_V6/`](Analysis/4_Mtb_RecombEvent_Explore_V6/) | Gene conversion event detection, mapping to paralogous regions, and visualization |
| [`Analysis/4_Mtb_ParalogRegions_MutPatternComparison_V4/`](Analysis/4_Mtb_ParalogRegions_MutPatternComparison_V4/) | Mutation pattern comparison between paralogous and non-paralogous regions |
| [`Analysis/4_Mtb_DNARepairGenes_LoF_Analysis_V1/`](Analysis/4_Mtb_DNARepairGenes_LoF_Analysis_V1/) | Loss-of-function analysis of DNA repair genes |
| [`Analysis/5_Mtb_EpitopeAnalysis_V5/`](Analysis/5_Mtb_EpitopeAnalysis_V5/) | Antigenic variation and HLA epitope analysis including PPE18 |
| [`Analysis/7_TGEN_GCVerf_Part1_SRWGS_Processing/`](Analysis/7_TGEN_GCVerf_Part1_SRWGS_Processing/) | Short-read WGS processing for the TGEN-937-SR validation dataset |
| [`Analysis/7_TGEN_GCVerf_Part2_OrgDataForSelectedIsolates/`](Analysis/7_TGEN_GCVerf_Part2_OrgDataForSelectedIsolates/) | Data organization for the 22 selected validation isolates |
| [`Analysis/7_TGEN_GCVerf_Part3_AnalyzeValidationData/`](Analysis/7_TGEN_GCVerf_Part3_AnalyzeValidationData/) | Analysis of the TBP-22-LR PacBio HiFi validation dataset |

---

## `gcutils` Python Package

The [`gcutils/`](gcutils/) directory contains a small Python package of utility functions developed for this project. It is used throughout the analysis notebooks for tasks such as processing Gubbins output, homology map analysis, and genome visualization.

It is not published on PyPI and must be installed directly from this repository:

```bash
pip install .
```

If you are modifying the package code, install in editable mode so changes are reflected immediately without reinstalling:

```bash
pip install -e .
```

---

## Installation

All dependencies for the downstream Python-based analysis notebooks can be installed via [Conda](https://docs.conda.io/en/latest/):

```bash
# 1) Clone the repository
git clone https://github.com/maxgmarin/mtb-geneconv-manuscript
cd mtb-geneconv-manuscript/

# 2) Create and activate the conda environment
conda create --name bfds_v1 -c conda-forge -c bioconda \
    python=3.10 pandas numpy matplotlib seaborn scipy tqdm pip \
    bioframe screed mmh3 pycirclize ete3 biopython plotly \
    ipykernel samtools minimap2 bcftools vcftools snakemake=7.2

conda activate bfds_v1

# 3) Register the kernel for Jupyter
python -m ipykernel install --user

# 4) Install the gcutils utility package in editable mode
pip install -e .
```

> The Snakemake workflows in `Snakemake_Workflows/` use their own per-rule Conda environments defined in each workflow's `CondaEnvs/` or `envs/` subdirectory.

---


## Citation

Marin MG, Quinones-Olvera N, Jin H, Harris MA, Jeffrey BM, Rosenthal A, Murphy KC, Sassetti C, Li H, Farhat MR. Gene conversion is a key driver of diversity hotspots in *Mycobacterium tuberculosis* antigens and virulence-associated loci. *bioRxiv* 2026. https://doi.org/10.64898/2026.02.26.708061
