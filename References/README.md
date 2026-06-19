# References

This directory contains H37Rv reference genome annotations, masking files, and other small reference datasets used throughout the analysis. All files are derived from publicly available sources (NCBI, WHO, or published literature) and are versioned here for reproducibility.

## Contents

| Path | Description |
|---|---|
| `190927_H37rv_GeneAnnotationsAndLists/` | Curated gene lists for PE/PPE family genes (Ates 2020 annotations), ESX secretion system genes, and standard exclusion lists (Coscolla 2015) used for filtering |
| `190927_H37rv_ListOf_ESXgenes.tsv` | List of ESX (ESX-1 through ESX-5) secretion system gene IDs in H37Rv |
| `201027_H37rv_AnnotatedGenes_And_IntergenicRegions/` | H37Rv gene-level annotation tables with gene centers, lengths, and locus tags; used for genome-wide annotation of variants and diversity estimates |
| `H37Rv_GenomeWindows/` | H37Rv genome segmented into 1,000 bp non-overlapping windows (BED and annotated TSV formats) |
| `H37Rv_Longdust_LowComplexityRegions/` | Low-complexity region (LCR) calls from [longdust](https://github.com/lh3/longdust) under default and relaxed stringency settings (BED format) |
| `H37Rv_MappabilityAnalysis_K50E4/` | Short-read mappability map for H37Rv (k=50, e=4): per-position pileup map, bigWig track, and a BED of regions with mappability < 1 |
| `Mtb_H37Rv_MaskingSchemes/` | Combined masking BEDs merging paralogous/low-complexity regions (from homology mapping) with Coscolla 2015 gene exclusions |
| `WHO_MtbAMR_Catalog/` | WHO *M. tuberculosis* antimicrobial resistance variant catalog (cleaned CSV format) used for variant annotation |
