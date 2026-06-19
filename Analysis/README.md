# Analysis Notebooks

All downstream data processing and analysis in this study was implemented as Python-based [Jupyter](https://jupyter.org/) notebooks. Notebooks are organized into subdirectories by analysis topic, numbered roughly in the order they were run. Each notebook contains the code, parameters, and embedded outputs for a specific analysis step.

## Subdirectories

| Directory | Description |
|---|---|
| [`1_Reference_Preprocessing/`](1_Reference_Preprocessing/) | Processing of H37Rv reference annotations; generation of homology maps and low-complexity/low-mappability masking files used throughout the analysis |
| [`2_Mtb_Prior_AsmQCandComparison/`](2_Mtb_Prior_AsmQCandComparison/) | Assembly QC evaluation and summary statistics across long-read and short-read WGS datasets |
| [`3_Mtb_NucDivAnalysis_V5/`](3_Mtb_NucDivAnalysis_V5/) | Nucleotide diversity (π) analysis across the Mtb genome; processing of windowed diversity estimates and identification of high-diversity hotspot regions |
| [`4_Mtb_RecombEvent_Explore_V6/`](4_Mtb_RecombEvent_Explore_V6/) | Core gene conversion analysis: processing Gubbins recombination predictions, mapping events to paralogous regions, paralog network analysis, and visualization of gene conversion events |
| [`4_Mtb_ParalogRegions_MutPatternComparison_V4/`](4_Mtb_ParalogRegions_MutPatternComparison_V4/) | Comparison of mutation patterns between paralogous and non-paralogous genomic regions |
| [`4_Mtb_DNARepairGenes_LoF_Analysis_V1/`](4_Mtb_DNARepairGenes_LoF_Analysis_V1/) | Analysis of loss-of-function variants in DNA repair genes across the Mtb151CI dataset |
| [`5_Mtb_EpitopeAnalysis_V5/`](5_Mtb_EpitopeAnalysis_V5/) | Antigenic variation analysis: mapping gene conversion events to T-cell epitopes, HLA binding predictions for recombinant PPE18 sequences, and mutation enrichment analysis |
| [`7_TGEN_GCVerf_Part1_SRWGS_Processing/`](7_TGEN_GCVerf_Part1_SRWGS_Processing/) | Short-read WGS processing and Gubbins analysis for the TGEN-937-SR dataset (937 isolates screened for gene conversion signatures) |
| [`7_TGEN_GCVerf_Part2_OrgDataForSelectedIsolates/`](7_TGEN_GCVerf_Part2_OrgDataForSelectedIsolates/) | Organization and preparation of input data for the 22 isolates selected for PacBio HiFi resequencing validation |
| [`7_TGEN_GCVerf_Part3_AnalyzeValidationData/`](7_TGEN_GCVerf_Part3_AnalyzeValidationData/) | Analysis of the TBP-22-LR validation dataset: evaluating gene conversion events detected in complete assemblies against short-read predictions |

---

## Installation

A shared Conda environment was used to run all notebooks. To reproduce the analysis environment:

```bash
# Create the conda environment
conda create --name bfds_v1 -c conda-forge -c bioconda \
    python=3.10 pandas numpy matplotlib seaborn scipy tqdm pip \
    bioframe screed mmh3 pycirclize ete3 biopython plotly \
    ipykernel samtools minimap2 bcftools vcftools snakemake=7.2

conda activate bfds_v1

# Register the kernel for Jupyter
python -m ipykernel install --user

# Install the gcutils utility package (from repo root)
cd ..
pip install -e .
```

The `gcutils` package (in the `gcutils/` directory at the repo root) contains shared utility functions used across multiple notebooks. It must be installed before running the notebooks.
