# 0. Genome Metadata

This directory contains a summary of genomic data and their associated metadata used in this project.


---

## Datasets

Three genomic datasets were used in this study.

**Mtb151CI** is the primary dataset of 151 complete, circularized *M. tuberculosis* genomes spanning lineages 1–8, assembled from raw sequencing data drawn from multiple published studies. All assemblies are hybrid long-read/short-read assemblies (Flye + Pilon or Flye + Pilon + PolyPolish, depending on the sequencing platform) and were used for the main gene conversion analyses in the manuscript.

**TGEN-937-SR** is a collection of 937 *M. tuberculosis* Illumina short-read WGS samples (TGEN isolates) that were analyzed for partial signatures of gene conversion detectable with short-read data. From this screen, 22 isolates were selected for PacBio HiFi resequencing to validate a subset of candidate gene conversion events.

**TBP-22-LR** is the resequencing validation dataset comprising those 22 TGEN isolates sequenced with PacBio HiFi long-read WGS (Sequel II). This dataset was used to confirm putative gene conversion events initially identified from short-read data.

---

## Metadata Table — SD1 (Supplemental Data 1)

Per-sample metadata for all three datasets, along with raw sequencing data accessions (SRA/ENA run accessions and BioProject IDs), are provided in **Supplemental Data Table 1 (SD1)**:

For convience `SD1.xlsx` is also provided in this directory;
- **SD1.xlsx** — Excel workbook with one sheet per dataset and per sequencing data subset

SD1 includes the following sheets:

| Sheet | Contents |
|---|---|
| SD1 - Overview | Dataset descriptions and notes |
| Mtb151 dataset Overview (n=151) | Assembly metadata, lineage calls, and assembly pipeline used for each genome |
| TGEN-937-SR Dataset (n=937) | Short-read WGS accessions, lineage calls, and coverage metrics |
| TBP-22-LR Dataset (n=22) | Long-read and short-read accessions for resequencing validation isolates |
| Mtb151CI - SetA – PB RSII Subreads | Raw PacBio RSII read accessions for the Mtb151CI subset from Marin et al. (2022) (BioProject PRJNA719670) |
| Mtb151CI - SetB – PB Sequel2 Subreads | Raw PacBio Sequel II subread accessions for the Mtb151CI subset from TB Portals (PRJNA421446) |
| Mtb151CI - SetC – PB Sequel2 CCS (HiFi) | Raw PacBio HiFi CCS read accessions for the Mtb151CI subset from the TRUST cohort — Marin et al. (2025) |
| Mtb151CI - SetD – ONT (Hall 2022) | Raw Oxford Nanopore read accessions for the Mtb151CI subset from Hall et al. (2022) |
| Mtb151CI - SetE – ONT (Peker 2022) | Raw Oxford Nanopore read accessions for the Mtb151CI subset from Peker et al. (2022) |

---

## Genome Assembly Sequences

Complete genome assemblies for both datasets (Mtb151CI and TBP-22-LR) are deposited on Zenodo and are available in two equivalent formats:

- **`.agc`** — Highly compressed [AGC](https://github.com/refresh-bio/agc) archive for efficient storage and retrieval
- **`.tar.gz`** — Standard archive containing the same genomes as individual FASTA files

| File | Dataset | Format |
|---|---|---|
| `Marin2026.Mtb.CompleteHybridAsms.Mtb151CI.WiH37Rv.agc` | Mtb151CI | AGC archive |
| `Marin2026.Mtb.CompleteHybridAsms.Mtb151CI.tar.gz` | Mtb151CI | FASTA (tar.gz) |
| `Marin2026.Mtb.CompleteHybridAsms.TBP22CI.WiH37Rv.agc` | TBP-22-LR | AGC archive |
| `Marin2026.Mtb.CompleteHybridAsms.TBP22CI.tar.gz` | TBP-22-LR | FASTA (tar.gz) |

> **Zenodo record:** [https://doi.org/10.64898/2026.02.26.708061](https://doi.org/10.64898/2026.02.26.708061) *(update with Zenodo DOI)*

