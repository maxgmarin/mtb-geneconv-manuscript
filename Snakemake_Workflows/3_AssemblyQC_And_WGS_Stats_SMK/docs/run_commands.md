# Run commands

## Environment

```bash
source /home/mm774/miniforge/etc/profile.d/conda.sh
conda activate /n/data1/hms/dbmi/farhat/mm774/Projects/Mtb-GeneConv/NucFlag_Testing/envs/nucflag_env
# or: conda env create -f envs/environment.yml -n asmqc_env && conda activate asmqc_env
```

## Dry run

```bash
cd 3_AssemblyQC_And_WGS_Stats_SMK
snakemake -s AsmAndWGS.QC.smk -n
```

## Local execution

```bash
snakemake -s AsmAndWGS.QC.smk --cores 4 \
  --config samples_tsv=config/samples_example.tsv results_dir=/path/to/results
```

## HMS O2 SLURM cluster submission

```bash
snakemake -s AsmAndWGS.QC.smk --profile profiles/o2_slurm_generic --jobs 50 \
  --config samples_tsv=config/samples_example.tsv results_dir=/path/to/results
```

The profile uses `snakemake-executor-plugin-cluster-generic` (see `envs/environment.yml`)
to submit each rule as its own `sbatch` job, using each rule's `params.partition` and
`resources.mem_mb`/`runtime_slurm`/`threads`.

## Real example invocation (Mtb151 cohort, original working pipeline)

The command used to generate the Mtb151 results summarized in
`scripts/NucFlag_ResultsSummary/results/`, run from the original working directory
above:

```bash
snakemake -s workflow/Snakefile --profile workflow/profiles/o2_slurm_generic --jobs 50 \
  --config samples_tsv=<Mtb151 151-sample TSV> results_dir=<Mtb151 results dir>
```
