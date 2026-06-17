# Pipeline Stages — Mtb.WGA.Core.V3.smk

This document describes the 8 logical processing stages of the core WGA pipeline. All exact commands and parameters are in [`../Mtb.WGA.Core.V3.smk`](../Mtb.WGA.Core.V3.smk).

---

## Stage 1: Per-sample alignment & variant calling

Each complete assembly FASTA is aligned to the H37Rv reference genome using minimap2 in assembly-to-reference mode (`asm10`). Variants are called via two complementary approaches: paftools (from PAF alignments) and bcftools mpileup (from BAM alignments).

**Rules:** `MM2_AsmToH37rv`, `filter_MM2_paftools_VCF_GC3_PP_AlignTo_H37rv_RemoveIndelsGreaterThan15bp`, `bcftools_mpileup_MM2_AsmToH37rv`, `getAll_SNPpositions_mpileup_MM2_AsmToH37rv`

**Key tools & parameters:**
- `minimap2 -ax asm10` — alignment optimized for highly similar sequences (~1 divergence)
- `paftools.js call` — calls variants from PAF format
- `bcftools mpileup` + `bcftools call -c` — classical pileup-based variant calling
- Indels >15 bp are removed from the paftools VCF

**Outputs:** Per-sample BAM, PAF, paftools VCF (SNPs + indels ≤15bp; SNPs only), mpileup VCF, SNP position TSV

---

## Stage 2: Lineage calling

*Mtb* lineage is assigned for each isolate using the per-sample paftools VCF.

**Rules:** `FastLinCaller_AsmToH37Rv_MM2`

**Key tools:** `fast-lineage-caller`

**Outputs:** Per-sample lineage call TSV

---

## Stage 3: Reference homology mapping

The H37Rv reference genome is self-aligned using minimap2 with short k-mer parameters (k=19, w=19) to identify repetitive low-complexity (RLC) and paralogous low-complexity (PLC) regions. A reverse alignment (H37Rv → assembly, asm20 mode) is also produced for downstream liftover use.

**Rules:** `HomologyMapping_Ref_k19w19Param`, `MM2_H37RvToAsm_asm20_ForLiftOver`

**Key tools & parameters:**
- `minimap2 -k19 -w19` — short k-mer self-alignment to find near-identical homologous regions
- `minimap2 -cx asm20` — reverse alignment for liftover
- bedtools merge — consolidates overlapping homologous intervals into a BED mask

**Outputs:** PAF, SAM/BAM, variant calls between homologous copies, merged RLC/PLC BED mask files

---

## Stage 4: Multi-sample SNV merging & filtering

Per-sample mpileup VCFs are reconciled to a common set of positions (union of all SNP positions across the cohort), then merged into a single multi-sample VCF. The joint callset is filtered under four conditions to enable sensitivity analyses:

| Filter | Description |
|---|---|
| Unfiltered | All SNVs |
| 10% ambiguity | Remove sites with >10% missing data (`F_MISSING > 0.10`) |
| RLC removed | Exclude sites overlapping repetitive/low-complexity regions |
| PLC removed | Exclude sites overlapping paralogous low-complexity regions |
| RLC + low-mappability | Combined mask (10% ambiguity + RLC + low-mappability) |

**Rules:** `combineAll_SNPpositions_mpileup_MM2_AsmToH37rv`, `Filter_MM2_mpileup_VarCalling_To_OnlySNPpositionsInUnionOfAllSNPs`, `merge_BCFs_Renamed_PATHs_To_TXT`, `merge_All_BCFs_Renamed_To_VCF`, `filter_SNVs_10AmbThresh_MergeSNVs_mpileup`, `filter_SNVs_RLC_Regions_MergeSNVs_mpileup`, `filter_SNVs_PLC_Regions_MergeSNVs_mpileup`, `filter_SNVs_RLCandLowPMap_Regions_MergeSNVs_mpileup`

**Key tools:** `bcftools merge`, `bcftools view`, `bedtools intersect -v`, `bgzip`, `tabix`

**Outputs:** Multi-sample merged VCF (+ bgzipped + tabix indexed) under each filter condition; SNP position lists

---

## Stage 5: FASTA alignment generation

Filtered multi-sample VCFs are converted to FASTA alignments for phylogenetic inference. Additionally, per-sample consensus sequences (SNVs inserted into H37Rv) are generated and concatenated into a full-genome MSA used by Gubbins.

**Rules:** `convert_MergedVCF_To_FASTA_ALN` (and variants for each filter condition), `insert_SNVs_IntoRefGenome_PerSample_10AmbThresh`, `create_FullGenome_MSA_SNVs_10AmbThresh`

**Key tools & parameters:**
- `vcf2phylip.py -i -f -m -100` — converts VCF to FASTA; `-m -100` retains all sites
- `bcftools consensus` — inserts SNVs per sample into H37Rv FASTA
- `bioawk` — renames FASTA headers to sample IDs

**Outputs:** Multi-sample FASTA alignments (SNV-only and full-genome MSA) per filter condition

---

## Stage 6: Phylogenetic inference

Phylogenetic trees are inferred under multiple SNV filter conditions using both FastTree (for rapid exploratory trees) and IQ-TREE (for final bootstrap-supported trees).

**Rules:** `fasttree_GTR_from_MergedSNPs`, `fasttree_GTR_from_MergedSNVs_10AmbFilt`, `IQtree_GTR_from_MergedSNVs_10AmbFilt`, `IQtree_GTR_from_MergedSNVs_PLCFilt`, `IQtree_GTR_from_MergedSNVs_10AmbFilt_RLCandLowPmap_Filt`

**Key tools & parameters:**
- `FastTree -nt -gtr` — GTR model, nucleotide data
- `snp-sites` — extract invariant-site-free FASTA for ASC correction
- `iqtree -m GTR+ASC -nt 1 -bb 10000` — GTR with ascertainment bias correction, 10,000 ultrafast bootstrap replicates

**Outputs:** Newick tree files (.treefile) per filter condition

---

## Stage 7: Nucleotide diversity

Windowed nucleotide diversity (π) is calculated across the genome for the merged SNV callset.

**Rules:** `calculate_NucDiversity_MergeSNVs_mpileup`

**Key tools & parameters:**
- `vcftools --window-pi 1000` — sliding window π in 1,000 bp windows

**Outputs:** Per-window π TSV

---

## Stage 8: Recombination detection (Gubbins)

Putative recombination events are detected using Gubbins on the full-genome MSA (10% ambiguity filtered). The starting tree from Stage 6 (FastTree) is provided to Gubbins to guide iteration.

**Rules:** `run_Gubbins_v321_ExtensiveSearch_SmallWin_MS4_Wi_mpileup_SNVs_10AmbThresh`, `reformat_GubbinsEvents_v321_ExtensiveSearch_SmallWin_MS4_Wi_SNV_10AmbThresh`

**Key tools & parameters:**
- `run_gubbins.py --extensive-search --min-window-size 25 --max-window-size 1000 --min-snps 4` (Gubbins v3.2.1)
- `sed` + `gff2bed` — reformats Gubbins GFF output (chromosome renamed to NC_000962.3) to BED format

**Outputs:** Node-labeled phylogeny, per-branch recombination statistics CSV, recombination predictions GFF and BED
