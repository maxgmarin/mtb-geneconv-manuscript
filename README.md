# Project Code Repository: Analysis & characterization of gene conversion within the Mtb genome

This repository contains all code and bioinformatics pipelines for the ongoing analysis of complete Mtb genome assemblies <br>


## Contents
- [Installation](#Installation)
- [Snakemake Pipelines](#Snakemake-pipelines)
- [Data Analysis](#Data-Analysis)
- [Results](#Results)
- [License](#License)

## Installation
All dependencies needed to reproduce the downsteram (python based) analysis can be installed via [Conda](https://docs.conda.io/en/latest/) .
```
# 1) Clone repository
git clone https://github.com/maxgmarin/mtb-GE-analysis

# 2) Create a conda environment named 'CoreEnv_PG_V1'
cd mtb-GE-analysis/

#conda env create --file CondaEnvs/bfds_v1.yml -n bfds_v1

# 3) Activate environment (for data analysis)

### This ENV will have jupyter, python, common libraries & Snakemake (v7.2)

conda create --name bfds_v1 -c conda-forge -c bioconda python=3.10 pandas numpy matplotlib seaborn scipy tqdm pip bioframe screed mmh3 pycirclize ete3 biopython plotly ipykernel samtools minimap2 bcftools vcftools snakemake=7.2

conda activate bfds_v1 

python -m ipykernel install --user


# 4) Pip install (in editable mode) the `ge_analysis_utils` package
# NOTE: THis package contains utility scripts and functions used for analysis

pip install -e .
```





