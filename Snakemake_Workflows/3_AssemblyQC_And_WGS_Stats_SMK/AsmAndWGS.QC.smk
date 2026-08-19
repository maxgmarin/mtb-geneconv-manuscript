# Assembly QC and WGS-stats pipeline.

# See README.md for the pipeline overview and how to run it.

import os

import pandas as pd

configfile: "config/config.yaml"

samples_df = pd.read_csv(config["samples_tsv"], sep="\t", dtype=str)
samples_df = samples_df.set_index("SampleID", drop=False)

# Define alignment parameters and NucFlag config parameters depending on long-read data type
TECH_PARAMS = {
    "PacBio_HiFi": {"mm2_preset": "lr:hqae", "nucflag_config": "nucflag_final_hifi.toml"},
    "ONT_r941": {"mm2_preset": "map-ont", "nucflag_config": "nucflag_final_ont_r9.toml"},
    "PacBio_Subreads": {"mm2_preset": "map-ont", "nucflag_config": "nucflag_final_ont_r9.toml"},
}


def get_fastq(wildcards):
    return samples_df.loc[wildcards.sample, "LR_FQ_PATH"]


def get_assembly(wildcards):
    return samples_df.loc[wildcards.sample, "AssemblyPath"]


def get_technology(wildcards):
    return samples_df.loc[wildcards.sample, "LR_Technology"]


def get_mm2_preset(wildcards):
    return TECH_PARAMS[get_technology(wildcards)]["mm2_preset"]


def get_nucflag_config(wildcards):
    tech = get_technology(wildcards)
    return os.path.join(config["nucflag_config_dir"], TECH_PARAMS[tech]["nucflag_config"])


def format_slurm_time(minutes):
    """Format a minute count as SLURM's D-HH:MM:SS --time string."""
    minutes = int(minutes)
    days, rem_minutes = divmod(minutes, 24 * 60)
    hours, mins = divmod(rem_minutes, 60)
    return f"{days}-{hours:02d}:{mins:02d}:00"


rule all:
    input:
        expand(
            "{results_dir}/AsmAnalysis/{sample}/nucflag_output/{sample}.qv.bed",
            results_dir=config["results_dir"], sample=samples_df["SampleID"],
        ),
        expand(
            "{results_dir}/AsmAnalysis/{sample}/nucflag_output/NucFlag_Misassemblies_LiftoverToH37Rv/{sample}.misassemblies_liftover_to_H37Rv.bed",
            results_dir=config["results_dir"], sample=samples_df["SampleID"],
        ),
        expand(
            "{results_dir}/AsmAnalysis/{sample}/nucflag_output/NucFlagStatus_Eval_37NucDivHotspots/{sample}.hotspot_status.length.annotated.bed",
            results_dir=config["results_dir"], sample=samples_df["SampleID"],
        ),
        expand(
            "{results_dir}/AsmAnalysis/{sample}/nucflag_output/NucFlagStatus_Eval_37NucDivHotspots/{sample}.hotspot_status.count.annotated.bed",
            results_dir=config["results_dir"], sample=samples_df["SampleID"],
        ),
        expand(
            "{results_dir}/AsmAnalysis/{sample}/nucflag_output/NucFlagStatus_Eval_ParalogousRegions/{sample}.pr_status.length.annotated.bed",
            results_dir=config["results_dir"], sample=samples_df["SampleID"],
        ),
        expand(
            "{results_dir}/AsmAnalysis/{sample}/nucflag_output/NucFlagStatus_Eval_ParalogousRegions/{sample}.pr_status.count.annotated.bed",
            results_dir=config["results_dir"], sample=samples_df["SampleID"],
        ),
        # Primary-only LR-to-own-assembly alignment, useful on its own (e.g. loading into
        # IGV) independent of the NucFlag calls, plus flagstat/coverage QC on both the
        # full and primary-only BAMs -- without requesting these explicitly, qc_full_bam/
        # qc_filtered_bam would never run since nothing else depends on their outputs.
        expand(
            "{results_dir}/AsmAnalysis/{sample}/alignment/{sample}.LR.AlnToOwnAssembly.primary_only.bam",
            results_dir=config["results_dir"], sample=samples_df["SampleID"],
        ),
        expand(
            "{results_dir}/AsmAnalysis/{sample}/alignment/{sample}.LR.AlnToOwnAssembly.flagstat.txt",
            results_dir=config["results_dir"], sample=samples_df["SampleID"],
        ),
        expand(
            "{results_dir}/AsmAnalysis/{sample}/alignment/{sample}.LR.AlnToOwnAssembly.coverage.txt",
            results_dir=config["results_dir"], sample=samples_df["SampleID"],
        ),
        expand(
            "{results_dir}/AsmAnalysis/{sample}/alignment/{sample}.LR.AlnToOwnAssembly.primary_only.flagstat.txt",
            results_dir=config["results_dir"], sample=samples_df["SampleID"],
        ),
        expand(
            "{results_dir}/AsmAnalysis/{sample}/alignment/{sample}.LR.AlnToOwnAssembly.primary_only.coverage.txt",
            results_dir=config["results_dir"], sample=samples_df["SampleID"],
        ),


### Stage 1: align long reads to each sample's own assembly ####################

rule stage_assembly:
    input:
        asm=get_assembly,
    output:
        fasta="{results_dir}/AsmAnalysis/{sample}/assembly/{sample}.fasta",
    params:
        partition=config["partition"],
        slurm_log_dir=config["slurm_log_dir"],
    resources:
        mem_mb=500, runtime=5, runtime_slurm="0-00:05:00",
    threads: 1
    shell:
        "ln -sf {input.asm} {output.fasta}"


rule faidx_assembly:
    input:
        fasta="{results_dir}/AsmAnalysis/{sample}/assembly/{sample}.fasta",
    output:
        fai="{results_dir}/AsmAnalysis/{sample}/assembly/{sample}.fasta.fai",
    params:
        partition=config["partition"],
        slurm_log_dir=config["slurm_log_dir"],
    resources:
        mem_mb=1000, runtime=5, runtime_slurm="0-00:05:00",
    threads: 1
    shell:
        "samtools faidx {input.fasta}"


def minimap2_runtime_minutes(wildcards, input):
    return max(90, int(input.size_mb / 250))


def minimap2_runtime_slurm(wildcards, input):
    return format_slurm_time(minimap2_runtime_minutes(wildcards, input))


rule minimap2_align_sort:
    input:
        asm="{results_dir}/AsmAnalysis/{sample}/assembly/{sample}.fasta",
        fai="{results_dir}/AsmAnalysis/{sample}/assembly/{sample}.fasta.fai",
        fq=get_fastq,
    output:
        bam="{results_dir}/AsmAnalysis/{sample}/alignment/{sample}.LR.AlnToOwnAssembly.bam",
        bai="{results_dir}/AsmAnalysis/{sample}/alignment/{sample}.LR.AlnToOwnAssembly.bam.bai",
        sam=temp("{results_dir}/AsmAnalysis/{sample}/alignment/{sample}.LR.AlnToOwnAssembly.sam"),
    params:
        preset=get_mm2_preset,
        partition=config["partition"],
        slurm_log_dir=config["slurm_log_dir"],
    resources:
        # Flat, not size-scaled: minimap2 processes reads in fixed-size batches, so peak
        # RSS plateaus rather than scaling with total input size (observed 3.8-5.5GB
        # across the largest FASTQs in this project's own batches; 10000 gives ~2x margin).
        mem_mb=10000,
        runtime=minimap2_runtime_minutes,
        runtime_slurm=minimap2_runtime_slurm,
    threads: 8
    log:
        "{results_dir}/AsmAnalysis/{sample}/alignment/{sample}.LR.AlnToOwnAssembly.minimap2.log",
    shell:
        # --MD --eqx required for NucFlag's per-base mismatch/identity pileup tracks.
        """
        minimap2 -t {threads} -ax {params.preset} --MD --eqx {input.asm} {input.fq} > {output.sam} 2> {log}
        samtools sort -@ {threads} -O bam -o {output.bam} {output.sam}
        samtools index {output.bam}
        """


def bam_qc_mem_mb(wildcards, input):
    return max(1000, int(0.15 * input.size_mb))


def bam_qc_runtime_minutes(wildcards, input):
    return max(30, int(input.size_mb / 300))


def bam_qc_runtime_slurm(wildcards, input):
    return format_slurm_time(bam_qc_runtime_minutes(wildcards, input))


rule qc_full_bam:
    input:
        bam="{results_dir}/AsmAnalysis/{sample}/alignment/{sample}.LR.AlnToOwnAssembly.bam",
        bai="{results_dir}/AsmAnalysis/{sample}/alignment/{sample}.LR.AlnToOwnAssembly.bam.bai",
    output:
        flagstat="{results_dir}/AsmAnalysis/{sample}/alignment/{sample}.LR.AlnToOwnAssembly.flagstat.txt",
        coverage="{results_dir}/AsmAnalysis/{sample}/alignment/{sample}.LR.AlnToOwnAssembly.coverage.txt",
    params:
        partition=config["partition"],
        slurm_log_dir=config["slurm_log_dir"],
    resources:
        mem_mb=bam_qc_mem_mb, runtime=bam_qc_runtime_minutes, runtime_slurm=bam_qc_runtime_slurm,
    threads: 1
    shell:
        "samtools flagstat {input.bam} > {output.flagstat} && "
        "samtools coverage {input.bam} > {output.coverage}"


rule filter_and_index_primary_only:
    # -F 2308 drops unmapped/secondary/supplementary alignments.
    input:
        bam="{results_dir}/AsmAnalysis/{sample}/alignment/{sample}.LR.AlnToOwnAssembly.bam",
        bai="{results_dir}/AsmAnalysis/{sample}/alignment/{sample}.LR.AlnToOwnAssembly.bam.bai",
    output:
        bam="{results_dir}/AsmAnalysis/{sample}/alignment/{sample}.LR.AlnToOwnAssembly.primary_only.bam",
        bai="{results_dir}/AsmAnalysis/{sample}/alignment/{sample}.LR.AlnToOwnAssembly.primary_only.bam.bai",
    params:
        partition=config["partition"],
        slurm_log_dir=config["slurm_log_dir"],
    resources:
        mem_mb=bam_qc_mem_mb, runtime=bam_qc_runtime_minutes, runtime_slurm=bam_qc_runtime_slurm,
    threads: 1
    shell:
        """
        samtools view -b -F 2308 -o {output.bam} {input.bam}
        samtools index {output.bam}
        """


rule qc_filtered_bam:
    input:
        bam="{results_dir}/AsmAnalysis/{sample}/alignment/{sample}.LR.AlnToOwnAssembly.primary_only.bam",
        bai="{results_dir}/AsmAnalysis/{sample}/alignment/{sample}.LR.AlnToOwnAssembly.primary_only.bam.bai",
    output:
        flagstat="{results_dir}/AsmAnalysis/{sample}/alignment/{sample}.LR.AlnToOwnAssembly.primary_only.flagstat.txt",
        coverage="{results_dir}/AsmAnalysis/{sample}/alignment/{sample}.LR.AlnToOwnAssembly.primary_only.coverage.txt",
    params:
        partition=config["partition"],
        slurm_log_dir=config["slurm_log_dir"],
    resources:
        mem_mb=bam_qc_mem_mb, runtime=bam_qc_runtime_minutes, runtime_slurm=bam_qc_runtime_slurm,
    threads: 1
    shell:
        "samtools flagstat {input.bam} > {output.flagstat} && "
        "samtools coverage {input.bam} > {output.coverage}"


### Stage 2: run NucFlag on the primary-only alignment ##########################

rule build_boundary_mask:
    # Masks the first/last N bp of the contig -- a circular-genome linearization
    # artifact, not something `[general] ignore_boundaries` in the NucFlag config
    # itself handles despite the similar name. mask_bp is a fixed pipeline
    # implementation detail, not exposed via config.yaml.
    input:
        fai="{results_dir}/AsmAnalysis/{sample}/assembly/{sample}.fasta.fai",
    output:
        bed="{results_dir}/AsmAnalysis/{sample}/nucflag_output/ignore_boundary.bed",
    params:
        mask_bp=3,
        partition=config["partition"],
        slurm_log_dir=config["slurm_log_dir"],
    resources:
        mem_mb=500, runtime=5, runtime_slurm="0-00:05:00",
    script:
        "scripts/build_boundary_mask.py"


def nucflag_call_mem_mb(wildcards, input):
    # Coefficient tuned against the deepest-coverage sample observed in this project
    # (PacBio_Subreads with meandepth ~3400x on a 4.4Mb genome) -- NucFlag's per-position
    # pileup memory scales disproportionately to raw BAM size on very-high-depth samples.
    return max(4000, int(4.0 * input.size_mb))


def nucflag_call_runtime_minutes(wildcards, input):
    return max(30, int(input.size_mb / 500))


def nucflag_call_runtime_slurm(wildcards, input):
    return format_slurm_time(nucflag_call_runtime_minutes(wildcards, input))


rule nucflag_call:
    # Always plots with --overlap_calls (calls drawn on the coverage plot rather than a
    # separate track); confirmed ad hoc that this only changes plot rendering, never the
    # misassembly classification itself.
    input:
        bam="{results_dir}/AsmAnalysis/{sample}/alignment/{sample}.LR.AlnToOwnAssembly.primary_only.bam",
        bai="{results_dir}/AsmAnalysis/{sample}/alignment/{sample}.LR.AlnToOwnAssembly.primary_only.bam.bai",
        asm="{results_dir}/AsmAnalysis/{sample}/assembly/{sample}.fasta",
        fai="{results_dir}/AsmAnalysis/{sample}/assembly/{sample}.fasta.fai",
        ignore_bed="{results_dir}/AsmAnalysis/{sample}/nucflag_output/ignore_boundary.bed",
    output:
        bed="{results_dir}/AsmAnalysis/{sample}/nucflag_output/{sample}.misassemblies.bed",
        status="{results_dir}/AsmAnalysis/{sample}/nucflag_output/{sample}.status.bed",
        plot_dir=directory("{results_dir}/AsmAnalysis/{sample}/nucflag_output/{sample}_plots"),
    params:
        config=get_nucflag_config,
        partition=config["partition"],
        slurm_log_dir=config["slurm_log_dir"],
    resources:
        mem_mb=nucflag_call_mem_mb,
        runtime=nucflag_call_runtime_minutes,
        runtime_slurm=nucflag_call_runtime_slurm,
    threads: 1
    log:
        "{results_dir}/AsmAnalysis/{sample}/nucflag_output/{sample}.nucflag_call.log",
    shell:
        "nucflag call -i {input.bam} -f {input.asm} -c {params.config} "
        "--ignore_regions {input.ignore_bed} "
        "-o {output.bed} -s {output.status} -d {output.plot_dir} "
        "--overlap_calls "
        "-t {threads} -p {threads} > {log} 2>&1"


rule nucflag_qv:
    input:
        bed="{results_dir}/AsmAnalysis/{sample}/nucflag_output/{sample}.misassemblies.bed",
    output:
        qv="{results_dir}/AsmAnalysis/{sample}/nucflag_output/{sample}.qv.bed",
    params:
        partition=config["partition"],
        slurm_log_dir=config["slurm_log_dir"],
    resources:
        mem_mb=500, runtime=5, runtime_slurm="0-00:05:00",
    threads: 1
    shell:
        "nucflag qv -i {input.bed} -o {output.qv}"


### Stage 3: assembly<->H37Rv alignment (feeds both liftover directions) ########

rule mm2_asm_to_h37rv:
    # target=H37Rv, query=Assembly -- the orientation needed to lift misassembly calls
    # (assembly coordinates) to H37Rv coordinates.
    input:
        asm="{results_dir}/AsmAnalysis/{sample}/assembly/{sample}.fasta",
        h37rv=config["h37rv_fasta"],
    output:
        bam="{results_dir}/AsmAnalysis/{sample}/MM2_AsmToH37Rv_asm20/{sample}.mm2.AsmToH37Rv.asm20.bam",
        bai="{results_dir}/AsmAnalysis/{sample}/MM2_AsmToH37Rv_asm20/{sample}.mm2.AsmToH37Rv.asm20.bam.bai",
        paf="{results_dir}/AsmAnalysis/{sample}/MM2_AsmToH37Rv_asm20/{sample}.mm2.AsmToH37Rv.asm20.paf",
        sam=temp("{results_dir}/AsmAnalysis/{sample}/MM2_AsmToH37Rv_asm20/{sample}.mm2.AsmToH37Rv.asm20.sam"),
    params:
        partition=config["partition"],
        slurm_log_dir=config["slurm_log_dir"],
    resources:
        mem_mb=2000, runtime=10, runtime_slurm="0-00:10:00",
    threads: 1
    log:
        "{results_dir}/AsmAnalysis/{sample}/MM2_AsmToH37Rv_asm20/{sample}.mm2.AsmToH37Rv.log",
    shell:
        """
        minimap2 -ax asm20 --MD --cs {input.h37rv} {input.asm} 2> {log} | awk '$1 ~ /^@/ || ($5 == 60)' > {output.sam}
        samtools view -bS {output.sam} | samtools sort -o {output.bam} -
        samtools index {output.bam}
        minimap2 -cx asm20 --cs {input.h37rv} {input.asm} 2>> {log} | awk '$1 ~ /^R/ || ($12 == 60)' > {output.paf}
        """


rule mm2_h37rv_to_asm_liftover:
    # target=Assembly, query=H37Rv -- the OPPOSITE orientation, needed to lift
    # H37Rv-coordinate regions of interest (NucDivHotspots, Paralogous Regions) onto
    # each sample's own assembly.
    input:
        asm="{results_dir}/AsmAnalysis/{sample}/assembly/{sample}.fasta",
        h37rv=config["h37rv_fasta"],
    output:
        paf="{results_dir}/AsmAnalysis/{sample}/MM2_H37RvToAsm_asm20_ForLiftOver/{sample}.mm2.H37RvToAsm.asm20.paf",
    params:
        partition=config["partition"],
        slurm_log_dir=config["slurm_log_dir"],
    resources:
        mem_mb=2000, runtime=10, runtime_slurm="0-00:10:00",
    threads: 1
    log:
        "{results_dir}/AsmAnalysis/{sample}/MM2_H37RvToAsm_asm20_ForLiftOver/{sample}.mm2.H37RvToAsm.log",
    shell:
        "minimap2 -cx asm20 --cs {input.asm} {input.h37rv} 2> {log} | "
        "awk '$1 ~ /^R/ || ($12 == 60)' > {output.paf}"


### Stage 4a: lift misassembly calls (assembly coords) -> H37Rv #################

rule misassemblies_only_bed:
    input:
        bed="{results_dir}/AsmAnalysis/{sample}/nucflag_output/{sample}.misassemblies.bed",
    output:
        bed="{results_dir}/AsmAnalysis/{sample}/nucflag_output/NucFlag_Misassemblies_LiftoverToH37Rv/{sample}.misassemblies_only.bed",
    params:
        partition=config["partition"],
        slurm_log_dir=config["slurm_log_dir"],
    resources:
        mem_mb=500, runtime=5, runtime_slurm="0-00:05:00",
    threads: 1
    shell:
        """awk -F'\\t' 'BEGIN{{OFS="\\t"}} NR>1 && $4!="correct" {{print $1,$2,$3,$4}}' {input.bed} > {output.bed}"""


rule liftover_misassemblies_raw:
    input:
        paf="{results_dir}/AsmAnalysis/{sample}/MM2_AsmToH37Rv_asm20/{sample}.mm2.AsmToH37Rv.asm20.paf",
        bed="{results_dir}/AsmAnalysis/{sample}/nucflag_output/NucFlag_Misassemblies_LiftoverToH37Rv/{sample}.misassemblies_only.bed",
    output:
        bed=temp("{results_dir}/AsmAnalysis/{sample}/nucflag_output/NucFlag_Misassemblies_LiftoverToH37Rv/{sample}.misassemblies_liftover_to_H37Rv.raw.bed"),
    params:
        partition=config["partition"],
        slurm_log_dir=config["slurm_log_dir"],
    resources:
        mem_mb=500, runtime=5, runtime_slurm="0-00:05:00",
    threads: 1
    shell:
        "paftools.js liftover {input.paf} {input.bed} > {output.bed}"


rule liftover_misassemblies_to_h37rv:
    input:
        liftover_raw="{results_dir}/AsmAnalysis/{sample}/nucflag_output/NucFlag_Misassemblies_LiftoverToH37Rv/{sample}.misassemblies_liftover_to_H37Rv.raw.bed",
        orig_bed="{results_dir}/AsmAnalysis/{sample}/nucflag_output/NucFlag_Misassemblies_LiftoverToH37Rv/{sample}.misassemblies_only.bed",
    output:
        bed="{results_dir}/AsmAnalysis/{sample}/nucflag_output/NucFlag_Misassemblies_LiftoverToH37Rv/{sample}.misassemblies_liftover_to_H37Rv.bed",
    params:
        partition=config["partition"],
        slurm_log_dir=config["slurm_log_dir"],
    resources:
        mem_mb=500, runtime=5, runtime_slurm="0-00:05:00",
    threads: 1
    log:
        "{results_dir}/AsmAnalysis/{sample}/nucflag_output/NucFlag_Misassemblies_LiftoverToH37Rv/{sample}.merge_liftover.log",
    script:
        "scripts/merge_liftover_with_metadata.py"


### Stage 4b: lift NucDivHotspots windows (H37Rv coords) -> each assembly #######

rule make_nucdiv_hotspots_bed:
    # One-time (not per-sample) conversion of the NucDivHotspots TSV into a header-free
    # BED-with-metadata file.
    input:
        tsv="resources/NucDivHotspots.37.1kb.tsv",
    output:
        bed="{results_dir}/_shared/NucDivHotspots.37.1kb.bed",
    params:
        partition=config["partition"],
        slurm_log_dir=config["slurm_log_dir"],
    resources:
        mem_mb=500, runtime=5, runtime_slurm="0-00:05:00",
    threads: 1
    shell:
        """awk -F'\\t' 'BEGIN{{OFS="\\t"}} NR>1' {input.tsv} > {output.bed}"""


rule liftover_nucdiv_hotspots_raw:
    input:
        paf="{results_dir}/AsmAnalysis/{sample}/MM2_H37RvToAsm_asm20_ForLiftOver/{sample}.mm2.H37RvToAsm.asm20.paf",
        bed="{results_dir}/_shared/NucDivHotspots.37.1kb.bed",
    output:
        bed=temp("{results_dir}/AsmAnalysis/{sample}/nucflag_output/37NucDivHotspots_LiftoverToAsm/{sample}.NucDivHotspots.liftover_to_Asm.raw.bed"),
    params:
        partition=config["partition"],
        slurm_log_dir=config["slurm_log_dir"],
    resources:
        mem_mb=500, runtime=5, runtime_slurm="0-00:05:00",
    threads: 1
    shell:
        "paftools.js liftover {input.paf} {input.bed} > {output.bed}"


rule liftover_nucdiv_hotspots_to_asm:
    input:
        liftover_raw="{results_dir}/AsmAnalysis/{sample}/nucflag_output/37NucDivHotspots_LiftoverToAsm/{sample}.NucDivHotspots.liftover_to_Asm.raw.bed",
        orig_bed="{results_dir}/_shared/NucDivHotspots.37.1kb.bed",
    output:
        bed="{results_dir}/AsmAnalysis/{sample}/nucflag_output/37NucDivHotspots_LiftoverToAsm/{sample}.NucDivHotspots.liftover_to_Asm.bed",
    params:
        partition=config["partition"],
        slurm_log_dir=config["slurm_log_dir"],
    resources:
        mem_mb=500, runtime=5, runtime_slurm="0-00:05:00",
    threads: 1
    log:
        "{results_dir}/AsmAnalysis/{sample}/nucflag_output/37NucDivHotspots_LiftoverToAsm/{sample}.merge_liftover.log",
    script:
        "scripts/merge_liftover_with_metadata.py"


### Stage 4c: lift Paralogous Regions (H37Rv coords) -> each assembly ###########

rule normalize_pr_regions_bed:
    # One-time (not per-sample) conversion of the Paralogous Regions TSV into a
    # header-free BED-with-metadata file (real column reorder, since the source TSV's
    # Chr/Start/End are columns 2-4, not 1-3).
    input:
        tsv=config["h37rv_paralogous_regions_tsv"],
    output:
        bed="{results_dir}/_shared/H37Rv.ParalogousRegions.bed",
    params:
        partition=config["partition"],
        slurm_log_dir=config["slurm_log_dir"],
    resources:
        mem_mb=500, runtime=5, runtime_slurm="0-00:05:00",
    threads: 1
    script:
        "scripts/normalize_pr_regions_bed.py"


rule liftover_pr_regions_raw:
    input:
        paf="{results_dir}/AsmAnalysis/{sample}/MM2_H37RvToAsm_asm20_ForLiftOver/{sample}.mm2.H37RvToAsm.asm20.paf",
        bed="{results_dir}/_shared/H37Rv.ParalogousRegions.bed",
    output:
        bed=temp("{results_dir}/AsmAnalysis/{sample}/nucflag_output/ParalogousRegions_LiftoverToAsm/{sample}.ParalogousRegions.liftover_to_Asm.raw.bed"),
    params:
        partition=config["partition"],
        slurm_log_dir=config["slurm_log_dir"],
    resources:
        mem_mb=500, runtime=5, runtime_slurm="0-00:05:00",
    threads: 1
    shell:
        "paftools.js liftover {input.paf} {input.bed} > {output.bed}"


rule liftover_pr_regions_to_asm:
    input:
        liftover_raw="{results_dir}/AsmAnalysis/{sample}/nucflag_output/ParalogousRegions_LiftoverToAsm/{sample}.ParalogousRegions.liftover_to_Asm.raw.bed",
        orig_bed="{results_dir}/_shared/H37Rv.ParalogousRegions.bed",
    output:
        bed="{results_dir}/AsmAnalysis/{sample}/nucflag_output/ParalogousRegions_LiftoverToAsm/{sample}.ParalogousRegions.liftover_to_Asm.bed",
    params:
        partition=config["partition"],
        slurm_log_dir=config["slurm_log_dir"],
    resources:
        mem_mb=500, runtime=5, runtime_slurm="0-00:05:00",
    threads: 1
    log:
        "{results_dir}/AsmAnalysis/{sample}/nucflag_output/ParalogousRegions_LiftoverToAsm/{sample}.merge_liftover.log",
    script:
        "scripts/merge_liftover_with_metadata.py"


### Stage 5: overlap misassembly calls against both region sets, then annotate ##

rule nucflag_status_hotspots:
    # Reuses NucFlag's own `nucflag status` subcommand to intersect a sample's full
    # misassembly calls against the NucDivHotspots windows (already lifted to this
    # sample's assembly coordinates), rather than reimplementing interval overlap logic.
    # Two output modes: "length" (% of each window covered by each call type + a QV
    # verdict) and "count" (literal count of call segments intersecting each window).
    input:
        bed="{results_dir}/AsmAnalysis/{sample}/nucflag_output/{sample}.misassemblies.bed",
        regions="{results_dir}/AsmAnalysis/{sample}/nucflag_output/37NucDivHotspots_LiftoverToAsm/{sample}.NucDivHotspots.liftover_to_Asm.bed",
    output:
        length="{results_dir}/AsmAnalysis/{sample}/nucflag_output/NucFlagStatus_Eval_37NucDivHotspots/{sample}.hotspot_status.length.bed",
        count_bed="{results_dir}/AsmAnalysis/{sample}/nucflag_output/NucFlagStatus_Eval_37NucDivHotspots/{sample}.hotspot_status.count.bed",
    params:
        partition=config["partition"],
        slurm_log_dir=config["slurm_log_dir"],
    resources:
        mem_mb=500, runtime=5, runtime_slurm="0-00:05:00",
    threads: 1
    log:
        "{results_dir}/AsmAnalysis/{sample}/nucflag_output/NucFlagStatus_Eval_37NucDivHotspots/{sample}.nucflag_status.log",
    shell:
        """
        nucflag status -i {input.bed} -b {input.regions} -g region -m length -o {output.length} > {log} 2>&1
        nucflag status -i {input.bed} -b {input.regions} -g region -m count -o {output.count_bed} >> {log} 2>&1
        """


rule annotate_hotspot_status:
    # Rejoins the original H37Rv coordinates + NucDivHotspots metadata onto each status
    # row, anchored on the FULL 37-window list (not just windows that happened to lift
    # over) -- paftools.js liftover silently drops rows that fail to lift, so anchoring
    # on the liftover/status file instead would silently omit those windows. See the
    # `Lifted` boolean column in the output.
    input:
        length="{results_dir}/AsmAnalysis/{sample}/nucflag_output/NucFlagStatus_Eval_37NucDivHotspots/{sample}.hotspot_status.length.bed",
        count_bed="{results_dir}/AsmAnalysis/{sample}/nucflag_output/NucFlagStatus_Eval_37NucDivHotspots/{sample}.hotspot_status.count.bed",
        liftover="{results_dir}/AsmAnalysis/{sample}/nucflag_output/37NucDivHotspots_LiftoverToAsm/{sample}.NucDivHotspots.liftover_to_Asm.bed",
        full_hotspots="{results_dir}/_shared/NucDivHotspots.37.1kb.bed",
    output:
        length_annotated="{results_dir}/AsmAnalysis/{sample}/nucflag_output/NucFlagStatus_Eval_37NucDivHotspots/{sample}.hotspot_status.length.annotated.bed",
        count_annotated="{results_dir}/AsmAnalysis/{sample}/nucflag_output/NucFlagStatus_Eval_37NucDivHotspots/{sample}.hotspot_status.count.annotated.bed",
    params:
        partition=config["partition"],
        slurm_log_dir=config["slurm_log_dir"],
    resources:
        mem_mb=500, runtime=5, runtime_slurm="0-00:05:00",
    threads: 1
    log:
        "{results_dir}/AsmAnalysis/{sample}/nucflag_output/NucFlagStatus_Eval_37NucDivHotspots/{sample}.annotate_hotspot_status.log",
    script:
        "scripts/annotate_hotspot_status_with_nucdiv.py"


rule nucflag_status_pr_regions:
    # Same approach as nucflag_status_hotspots, applied to the Paralogous Regions
    # instead of NucDivHotspots.
    input:
        bed="{results_dir}/AsmAnalysis/{sample}/nucflag_output/{sample}.misassemblies.bed",
        regions="{results_dir}/AsmAnalysis/{sample}/nucflag_output/ParalogousRegions_LiftoverToAsm/{sample}.ParalogousRegions.liftover_to_Asm.bed",
    output:
        length="{results_dir}/AsmAnalysis/{sample}/nucflag_output/NucFlagStatus_Eval_ParalogousRegions/{sample}.pr_status.length.bed",
        count_bed="{results_dir}/AsmAnalysis/{sample}/nucflag_output/NucFlagStatus_Eval_ParalogousRegions/{sample}.pr_status.count.bed",
    params:
        partition=config["partition"],
        slurm_log_dir=config["slurm_log_dir"],
    resources:
        mem_mb=500, runtime=5, runtime_slurm="0-00:05:00",
    threads: 1
    log:
        "{results_dir}/AsmAnalysis/{sample}/nucflag_output/NucFlagStatus_Eval_ParalogousRegions/{sample}.nucflag_status.log",
    shell:
        """
        nucflag status -i {input.bed} -b {input.regions} -g region -m length -o {output.length} > {log} 2>&1
        nucflag status -i {input.bed} -b {input.regions} -g region -m count -o {output.count_bed} >> {log} 2>&1
        """


rule annotate_pr_status:
    # Rejoins the original H37Rv coordinates + PR metadata (HmRegionID, Overlap_Genes,
    # Length, PR_SetID) onto each status row, keyed on HmRegionID (a real unique ID),
    # anchored on the FULL PR list for the same reason as annotate_hotspot_status.
    input:
        length="{results_dir}/AsmAnalysis/{sample}/nucflag_output/NucFlagStatus_Eval_ParalogousRegions/{sample}.pr_status.length.bed",
        count_bed="{results_dir}/AsmAnalysis/{sample}/nucflag_output/NucFlagStatus_Eval_ParalogousRegions/{sample}.pr_status.count.bed",
        liftover="{results_dir}/AsmAnalysis/{sample}/nucflag_output/ParalogousRegions_LiftoverToAsm/{sample}.ParalogousRegions.liftover_to_Asm.bed",
        full_prs="{results_dir}/_shared/H37Rv.ParalogousRegions.bed",
    output:
        length_annotated="{results_dir}/AsmAnalysis/{sample}/nucflag_output/NucFlagStatus_Eval_ParalogousRegions/{sample}.pr_status.length.annotated.bed",
        count_annotated="{results_dir}/AsmAnalysis/{sample}/nucflag_output/NucFlagStatus_Eval_ParalogousRegions/{sample}.pr_status.count.annotated.bed",
    params:
        partition=config["partition"],
        slurm_log_dir=config["slurm_log_dir"],
    resources:
        mem_mb=500, runtime=5, runtime_slurm="0-00:05:00",
    threads: 1
    log:
        "{results_dir}/AsmAnalysis/{sample}/nucflag_output/NucFlagStatus_Eval_ParalogousRegions/{sample}.annotate_pr_status.log",
    script:
        "scripts/annotate_pr_status_with_hmregion.py"
