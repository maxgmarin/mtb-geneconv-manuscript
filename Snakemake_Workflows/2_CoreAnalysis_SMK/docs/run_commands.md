# Exact Run Commands — Mtb.WGA.Core.V3.smk

This document contains the exact commands used to run the core WGA pipeline on each dataset. Input sample sheet TSVs are tracked in this repository — see the [Input Manifests](../README.md#input-manifests) section of the main README.

---

## Dataset 1: Mtb151CI (151 hybrid assemblies)

```bash
Farhat_Lab_Dir="/n/data1/hms/dbmi/farhat"
SMK_RepoDir="${Farhat_Lab_Dir}/mm774/Snakemake_Pipelines/Mtb-WGA-SMK"

targetOutput_Dir="${Farhat_Lab_Dir}/mm774/Projects/Mtb-WGA-SMK-Output/231121_MtbSetV3_151CI"
mkdir -p ${targetOutput_Dir}

inputConfigFile="${SMK_RepoDir}/Mtb-WGA.config.json"
input_ClusterConfig="${SMK_RepoDir}/Mtb-WGA.clusterConfig.json"

input_SampleInfo_TSV="${Farhat_Lab_Dir}/mm774/Projects/Mtb-VCI-MGM/Data/231121.InputAsmTSVs.MtbSetV3.151CI/231121.MtbSetV3.151CI.HybridAndSRAsm.FAPATHs.V1.tsv"

cd ${SMK_RepoDir}
```

### Dry run
```bash
snakemake -s Mtb.WGA.Core.V3.smk \
    --config output_dir=${targetOutput_Dir} inputSampleData_TSV=${input_SampleInfo_TSV} \
    --configfile ${inputConfigFile} \
    -np
```

### Local execution
```bash
snakemake -s Mtb.WGA.Core.V3.smk \
    --config output_dir=${targetOutput_Dir} inputSampleData_TSV=${input_SampleInfo_TSV} \
    --configfile ${inputConfigFile} \
    -p --cores 4 --rerun-incomplete
```

### HMS O2 Slurm cluster submission
```bash
mkdir -p ${targetOutput_Dir}/O2logs/cluster/

time snakemake -s Mtb.WGA.Core.V3.smk \
    --config output_dir=${targetOutput_Dir} inputSampleData_TSV=${input_SampleInfo_TSV} \
    --configfile ${inputConfigFile} \
    -p --use-conda -j 500 \
    --cluster-config ${input_ClusterConfig} \
    --cluster "sbatch -p {cluster.p} -n {cluster.n} -t {cluster.t} --mem {cluster.mem} \
               -o ${targetOutput_Dir}/O2logs/cluster/{cluster.o} \
               -e ${targetOutput_Dir}/O2logs/cluster/{cluster.e}" \
    --latency-wait 45 -k
```

---

## Dataset 2: TBP-22-LR (22 PacBio HiFi assemblies)

```bash
Farhat_Lab_Dir="/n/data1/hms/dbmi/farhat"
SMK_RepoDir="${Farhat_Lab_Dir}/mm774/Snakemake_Pipelines/Mtb-WGA-SMK"

targetOutput_Dir="${Farhat_Lab_Dir}/mm774/Projects/Mtb-WGA-SMK-Output/250730_TBP22_22CI_V2"
mkdir -p ${targetOutput_Dir}

inputConfigFile="${SMK_RepoDir}/Mtb-WGA.config.json"
input_ClusterConfig="${SMK_RepoDir}/Mtb-WGA.clusterConfig.json"

input_SampleInfo_TSV="${Farhat_Lab_Dir}/mm774/Projects/Mtb-GeneConv/mtb-GE-analysis/Data/TBP22.22CI.GCEVerfIsolates.Metadata/250801.TBP22.22CI.GCEVerfIsolates.Asm_LR_SR.InputPATHs.tsv"

cd ${SMK_RepoDir}
```

### Dry run
```bash
snakemake -s Mtb.WGA.Core.V3.smk \
    --config output_dir=${targetOutput_Dir} inputSampleData_TSV=${input_SampleInfo_TSV} \
    --configfile ${inputConfigFile} \
    -np
```

### Local execution
```bash
snakemake -s Mtb.WGA.Core.V3.smk \
    --config output_dir=${targetOutput_Dir} inputSampleData_TSV=${input_SampleInfo_TSV} \
    --configfile ${inputConfigFile} \
    -p --cores 4 --rerun-incomplete
```

### HMS O2 Slurm cluster submission
```bash
mkdir -p ${targetOutput_Dir}/O2logs/cluster/

time snakemake -s Mtb.WGA.Core.V3.smk \
    --config output_dir=${targetOutput_Dir} inputSampleData_TSV=${input_SampleInfo_TSV} \
    --configfile ${inputConfigFile} \
    -p --use-conda -j 500 \
    --cluster-config ${input_ClusterConfig} \
    --cluster "sbatch -p {cluster.p} -n {cluster.n} -t {cluster.t} --mem {cluster.mem} \
               -o ${targetOutput_Dir}/O2logs/cluster/{cluster.o} \
               -e ${targetOutput_Dir}/O2logs/cluster/{cluster.e}" \
    --latency-wait 45 -k
```
