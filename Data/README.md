# Data

This directory contains input metadata tables and processed results files used by the analysis notebooks in `Analysis/`. It is organized into subdirectories by dataset or analysis type.

## Contents

| Directory | Description |
|---|---|
| `220813_MtbEpitopes/` | IEDB T-cell epitope data mapped to H37Rv, including curated epitope sets from Lindestam 2016 and Panda 2024 with per-gene epitope mapping statistics |
| `221017_TBPortals_LRandSR_InputDataTSVs/` | Input WGS path manifests and metadata for the TB Portals long-read + short-read dataset (31 isolates, QC-pass) used for initial assembly |
| `231002.InputAsmTSVs.TBP22.29I.Complete/` | Assembly QC stats and FASTA path manifests for the 29-isolate TBP-22 complete hybrid assembly set |
| `231121_HybridMtbAsm_QCPass_Meta_Set3/` | Assembly summary and QC tables for the Mtb151CI dataset subsets (Hall 2022, Peker 2021, TRUST cohort, TB Portals Sequel II, combined 48-isolate SR set) |
| `231121.InputAsmTSVs.MtbSetV3.151CI/` | Input sample sheet and assembly summary for the final 151-isolate Mtb151CI dataset used as the primary analysis input to `2_CoreAnalysis_SMK` |
| `241024.MutSpectraAnalysis.MUSICAL/` | Mutational spectra inputs and MUSICAL output for SNPs in paralogous vs. non-paralogous regions of Mtb151CI |
| `241030.Mtb151CI.AllVariants.Anno.V1/` | All-variant TSV for Mtb151CI with CDS annotations and inter-SNP distances (compressed; used by recombination and nucleotide diversity notebooks) |
| `250910.PPE18.netMHCpanII.PredictionsForRecombinantSeqs/` | netMHCpanII HLA class II binding predictions for recombinant PPE18 amino acid sequences inferred from gene conversion events (27-allele and 7-allele Euro panels) |
| `H37Rv.NonUniqueSeqRegions.V1/` | Derived masking files for H37Rv: paralogous regions, local repeats, low-complexity regions (longdust), and low-mappability regions — the union mask used for SNV filtering |
| `TBP22.22CI.GCEVerfIsolates.Metadata/` | Metadata and input path manifests for the 22 TBP-22-LR validation isolates, including per-isolate gene conversion event mappings |
| `Tgen1K_WGS_RunMetadata/` | Short-read WGS run metadata and sample information for the TGEN-937-SR dataset |
