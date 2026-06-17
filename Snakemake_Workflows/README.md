# Snakemake Workflows

All major computational steps in this study were orchestrated using the [Snakemake](https://snakemake.readthedocs.io/) workflow management system and executed on the HMS O2 cluster (Slurm). Each subdirectory contains a self-contained workflow covering a distinct stage of the analysis, along with its configuration files, Conda environment definitions, and cluster submission configs.

| Directory | Description |
|---|---|
| [`1_GenomeAssembly_SMK/`](1_GenomeAssembly_SMK/) | Hybrid genome assembly from long-read + short-read WGS data (Flye + Pilon ± Medaka/PolyPolish) |
| [`2_CoreAnalysis_SMK/`](2_CoreAnalysis_SMK/) | Whole-genome analysis: reference alignment, variant calling, phylogenetics, nucleotide diversity, and recombination detection |
