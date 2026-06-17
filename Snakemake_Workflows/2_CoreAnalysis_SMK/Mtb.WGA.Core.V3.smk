
# Analysis pipeline of complete Mtb genome assemblies.

### This snakemake workflow defines the analysis workflow starting with genome metadata + a FASTA sequence of each genome assembly.




### Import Statements ###
import pandas as pd


### Define PATHs to files defined in thoe config file ###
refGenome_FA_PATH = config["RefGenome_FA_PATH"]
refGenome_GFF_PATH = config["RefGenome_GFF_PATH"]

H37rv_DnaA_FA_PATH = config["H37rv_DnaA_FA_PATH"]
H37rv_GBK_PATH = config["H37rv_GBK_PATH"]

H37Rv_AA_FA_PATH = "/n/data1/hms/dbmi/farhat/mm774/References/190619_Mycobrowser_H37rv_ReferenceFiles/Mycobacterium_tuberculosis_H37Rv_proteins_v3_TrimmedHeader.fasta"



# Define PATH of main OUTPUT directory
output_Dir = config["output_dir"]

i_log_outdir = output_Dir  + "/slurm_logs"


# Define analysis name to use a prefix for output files
#AnalysisName = config["analysis_name"]

# Define window sizes for nucleotide diversity calculation
NucDiv_WindowSizes_bp = ['1000']

# Read in data regarding input 
input_DataInfo_DF = pd.read_csv( config["inputSampleData_TSV"], sep='\t')


# Create a python list of Sample IDs
input_All_SampleIDs = list( input_DataInfo_DF["SampleID"].values )

SampleID_to_Asm_FA_Dict = dict(input_DataInfo_DF[['SampleID', 'HybridAsm_FA_PATH']].values)

print("List of input sampleIDs:", len(input_All_SampleIDs), input_All_SampleIDs)






############ Define Dicts w/ Short-read PE FASTQ Paths ############

# SampleID_To_FQ1_PathDict = input_DataInfo_DF.set_index('SampleID')['Illumina_PE_FQs_PATH'].str.split(";").str[0].to_dict()
# SampleID_To_FQ2_PathDict = input_DataInfo_DF.set_index('SampleID')['Illumina_PE_FQs_PATH'].str.split(";").str[1].to_dict()

# print("# of sampleIDs w/ SR-WGS FQ1 PATH DEFINED:", len(list(SampleID_To_FQ1_PathDict.keys())),)
# print("# of sampleIDs w/ SR-WGS FQ2 PATH DEFINED:", len(list(SampleID_To_FQ2_PathDict.keys())),)

#print(SampleID_To_FQ1_PathDict)
#print(SampleID_To_FQ2_PathDict)

print()

#################################################################




rule all:
    input:
        expand(output_Dir + "/AsmAnalysis/{sampleID}/VariantCallingVersusH37Rv/MM2_AsmToH37rv/{sampleID}.mm2.AsmToH37Rv.var.tsv", sampleID = input_All_SampleIDs),
        expand(output_Dir + "/AsmAnalysis/{sampleID}/VariantCallingVersusH37Rv/MM2_AsmToH37rv/{sampleID}.mm2.AsmToH37Rv.paf", sampleID = input_All_SampleIDs),


        expand(output_Dir + "/HomologyMapping/H37Rv_Ref_HomologyMapping_k19w19Param/H37Rv.HomologyMap.VizUpdated.bam", sampleID = input_All_SampleIDs),

        expand(output_Dir + "/AsmAnalysis/{sampleID}/FastANI/FastANI_AsmToH37Rv/{sampleID}.AsmToH37Rv.FastANI.txt", sampleID = input_All_SampleIDs),
        expand(output_Dir + "/AsmAnalysis/{sampleID}/LineageCalling/LinCall_Paftools_AsmToH37Rv/{sampleID}.AsmToH37Rv.lineage_call.tsv", sampleID = input_All_SampleIDs),


        output_Dir + "/Phylogenies/fasttree_mpileupSNVs_NoFilt/MM2.mpileup.call.Merged.SNVs.min-100.fasttree.newick",   
        output_Dir + "/Phylogenies/fasttree_mpileupSNVs_10AmbFilt/MM2.mpileup.call.Merged.SNVs.10AmbThresh.NoMask.min-100.fasttree.newick",

        expand(output_Dir + "/NucDiversity/NucDiv_SNVs_mpileup/MM2.mpileup.call.Merged.SNVs.NucDiv.{NucDiv_WindowSize}bp.windowed.pi", NucDiv_WindowSize = NucDiv_WindowSizes_bp),

        output_Dir + "/RecombDetection/Gubbins_v321_ExtSearch_SW_MS4_mpileup_SNVs_10AmbThresh/Gubbins.recombination_predictions.RenamedCHR.gff",


        expand(output_Dir + "/HomologyMapping/H37Rv_Ref_HomologyMapping_k19w19Param/H37Rv.HomologyMap.VizUpdated.bam", sampleID = input_All_SampleIDs),


        expand(output_Dir + "/AsmAnalysis/{sampleID}/VariantCallingVersusH37Rv/MM2_H37RvToAsm_asm20_ForLiftOver/{sampleID}.mm2.H37RvToAsm.asm20.paf", sampleID = input_All_SampleIDs),

        output_Dir + "/Asm_MergeSNPs_mpileup/AllSNPpositions.MM2.mpileup.SNPs.Union.AllSamples.tsv",
        
        output_Dir + "/Asm_MergeVar_mpileup/MM2.mpileup.call.Merged.SNVs.vcf",


        
        output_Dir + "/Asm_MergeSNPs_mpileup/AllSNPpositions.MM2.mpileup.SNPs.Union.AllSamples.tsv",
        

        output_Dir + "/Asm_MergeVar_mpileup/MM2.mpileup.call.Merged.SNVs.vcf",


        output_Dir + "/Phylogenies/fasttree_mpileupSNVs_NoFilt/MM2.mpileup.call.Merged.SNVs.min-100.fasttree.newick",   
        output_Dir + "/Phylogenies/fasttree_mpileupSNVs_10AmbFilt/MM2.mpileup.call.Merged.SNVs.10AmbThresh.NoMask.min-100.fasttree.newick",
        output_Dir + "/Phylogenies/iqtree_mpileupSNVs_PLCFilt/MM2.mpileup.call.Merged.SNVs.PLCMask.min-100.iq.treefile",


        expand(output_Dir + "/NucDiversity/NucDiv_SNVs_mpileup/MM2.mpileup.call.Merged.SNVs.NucDiv.{NucDiv_WindowSize}bp.windowed.pi", NucDiv_WindowSize = NucDiv_WindowSizes_bp),










#########################################################
#########################################################
## Analysis Versus H37rv (call variants against H37Rv) ##
#########################################################
#########################################################

##############################################################################
############ Minimap2: Assembly To H37rv Alignment & Variant Calling #########
##############################################################################

rule MM2_AsmToH37rv:
    input:
        i_Assembly_FA = lambda wildcards: SampleID_to_Asm_FA_Dict[wildcards.sampleID],
        H37rv_FA = refGenome_FA_PATH,
    output:
        MM2_AsmToH37rv_SAM = output_Dir + "/AsmAnalysis/{sampleID}/VariantCallingVersusH37Rv/MM2_AsmToH37rv/{sampleID}.mm2.AsmToH37Rv.sam",
        MM2_AsmToH37rv_PAF = output_Dir + "/AsmAnalysis/{sampleID}/VariantCallingVersusH37Rv/MM2_AsmToH37rv/{sampleID}.mm2.AsmToH37Rv.paf",
        MM2_AsmToH37rv_BAM = output_Dir + "/AsmAnalysis/{sampleID}/VariantCallingVersusH37Rv/MM2_AsmToH37rv/{sampleID}.mm2.AsmToH37Rv.bam",
        MM2_AsmToH37rv_bai = output_Dir + "/AsmAnalysis/{sampleID}/VariantCallingVersusH37Rv/MM2_AsmToH37rv/{sampleID}.mm2.AsmToH37Rv.bam.bai",
        MM2_AsmToH37rv_paftools_VCF = output_Dir + "/AsmAnalysis/{sampleID}/VariantCallingVersusH37Rv/MM2_AsmToH37rv/{sampleID}.mm2.AsmToH37Rv.paftools.vcf",
        MM2_AsmToH37rv_VarTSV = output_Dir + "/AsmAnalysis/{sampleID}/VariantCallingVersusH37Rv/MM2_AsmToH37rv/{sampleID}.mm2.AsmToH37Rv.var.tsv",
    conda: "CondaEnvs/mm2_v2_4_WiUtilities.V2.yml"
    threads: 1
    params:
        MM2_MinAlnLen_ForCoverage = 1000,
        MM2_MinAlnLen_ForVariantCalling = 1000,
        partition = 'short',
        slurm_log_dir = i_log_outdir,
    resources:
        mem_mb = lambda wildcards, attempt: (attempt * 1000) + 3000 , # Increment memory by 1000 mb per attempt 
        runtime_slurm = lambda wildcards, attempt: '0-00:' + str(attempt * 5) + ":00", #runtime = '0-00:10:00',
        runtime = lambda wildcards, attempt: (attempt * 5) , # minutes
    shell: 
        "minimap2 -ax asm10 --MD --cs {input.H37rv_FA} {input.i_Assembly_FA} | awk '$1 ~ /^@/ || ($5 == 60)' > {output.MM2_AsmToH37rv_SAM} \n"
        "samtools view -bS {output.MM2_AsmToH37rv_SAM} | samtools sort - > {output.MM2_AsmToH37rv_BAM} \n"
        "samtools index {output.MM2_AsmToH37rv_BAM} \n"
        "minimap2 -cx asm10 --cs {input.H37rv_FA} {input.i_Assembly_FA} | awk '$1 ~ /^R/ || ($12 == 60)' > {output.MM2_AsmToH37rv_PAF} \n"
        #"minimap2 -cx asm10 --cs {input.H37rv_FA} {input.i_Assembly_FA} | awk '$1 ~ /^R/ || ($12 == 60)' | sort -k6,6 -k8,8n | paftools.js call -s {wildcards.sampleID} -L {params.MM2_MinAlnLen_ForVariantCalling} -l {params.MM2_MinAlnLen_ForCoverage} -f {input.H37rv_FA} - > {output.MM2_AsmToH37rv_paftools_VCF} \n"
        #"minimap2 -cx asm10 --cs {input.H37rv_FA} {input.i_Assembly_FA} | awk '$1 ~ /^R/ || ($12 == 60)' | sort -k6,6 -k8,8n | paftools.js call -s {wildcards.sampleID} -L {params.MM2_MinAlnLen_ForVariantCalling} -l {params.MM2_MinAlnLen_ForCoverage} - > {output.MM2_AsmToH37rv_VarTSV} \n"
        " sort -k6,6 -k8,8n {output.MM2_AsmToH37rv_PAF} | /home/mm774/miniforge/envs/base2/bin/paftools.js call -s {wildcards.sampleID} -L {params.MM2_MinAlnLen_ForVariantCalling} -l {params.MM2_MinAlnLen_ForCoverage} -f {input.H37rv_FA} - > {output.MM2_AsmToH37rv_paftools_VCF} \n"
        " sort -k6,6 -k8,8n {output.MM2_AsmToH37rv_PAF} | /home/mm774/miniforge/envs/base2/bin/paftools.js call -s {wildcards.sampleID} -L {params.MM2_MinAlnLen_ForVariantCalling} -l {params.MM2_MinAlnLen_ForCoverage} - > {output.MM2_AsmToH37rv_VarTSV} \n"



rule filter_MM2_paftools_VCF_GC3_PP_AlignTo_H37rv_RemoveIndelsGreaterThan15bp:
    input:
        MM2_AsmToH37rv_paftools_VCF = output_Dir + "/AsmAnalysis/{sampleID}/VariantCallingVersusH37Rv/MM2_AsmToH37rv/{sampleID}.mm2.AsmToH37Rv.paftools.vcf",
    output:
        MM2_AsmToH37rv_paftools_VCF_SNPsAndINDELs_Lengths_1to15bp_Only = output_Dir + "/AsmAnalysis/{sampleID}/VariantCallingVersusH37Rv/MM2_AsmToH37rv/{sampleID}.mm2.AsmToH37rv.paftools.Lengths_1to15bp.vcf",    
        MM2_AsmToH37rv_paftools_VCF_SNPs_Only = output_Dir + "/AsmAnalysis/{sampleID}/VariantCallingVersusH37Rv/MM2_AsmToH37rv/{sampleID}.mm2.AsmToH37rv.paftools.SNPsOnly.vcf",    
    conda: "CondaEnvs/mm2_v2_4_WiUtilities.V2.yml"
    threads: 1
    params:
        partition = 'short',
        slurm_log_dir = i_log_outdir,
    resources:
        mem_mb = lambda wildcards, attempt: (attempt * 1000) + 3000 , # Increment memory by 1000 mb per attempt 
        runtime_slurm = lambda wildcards, attempt: '0-00:' + str(attempt * 5) + ":00", #runtime = '0-00:10:00',
        runtime = lambda wildcards, attempt: (attempt * 5) , # minutes
    shell:
        "bcftools view --types snps,indels -i 'abs(strlen(ALT)-strlen(REF))<=15' {input.MM2_AsmToH37rv_paftools_VCF} > {output.MM2_AsmToH37rv_paftools_VCF_SNPsAndINDELs_Lengths_1to15bp_Only} \n"
        "bcftools view --types snps {input.MM2_AsmToH37rv_paftools_VCF} > {output.MM2_AsmToH37rv_paftools_VCF_SNPs_Only} \n"








rule MM2_H37RvToAsm_asm20_ForLiftOver:
    input:
        i_Assembly_FA = lambda wildcards: SampleID_to_Asm_FA_Dict[wildcards.sampleID],
        H37rv_FA = refGenome_FA_PATH,
    output:
        MM2_H37rvToAsm_PAF = output_Dir + "/AsmAnalysis/{sampleID}/VariantCallingVersusH37Rv/MM2_H37RvToAsm_asm20_ForLiftOver/{sampleID}.mm2.H37RvToAsm.asm20.paf",
    conda: "CondaEnvs/mm2_v2_4_WiUtilities.V2.yml"
#    conda: "/home/mm774/conda3/envs/mm2_v2_4_WiUtilities/"
    threads: 1
    params:
        MM2_MinAlnLen_ForCoverage = 1000,
        MM2_MinAlnLen_ForVariantCalling = 1000,
        partition = 'short',
        slurm_log_dir = i_log_outdir,
    resources:
        mem_mb = lambda wildcards, attempt: (attempt * 1000) + 3000 , # Increment memory by 1000 mb per attempt 
        runtime_slurm = lambda wildcards, attempt: '0-00:' + str(attempt * 5) + ":00", #runtime = '0-00:10:00',
        runtime = lambda wildcards, attempt: (attempt * 5) , # minutes
    shell:
        #"conda activate /home/mm774/conda3/envs/mm2_v2_4_WiUtilities \n"
        "minimap2 -cx asm20 --cs {input.i_Assembly_FA} {input.H37rv_FA} | awk '$1 ~ /^R/ || ($12 == 60)' > {output.MM2_H37rvToAsm_PAF} \n"












##########################################################################
############ Subset SNPs from Mpileup for Phylogeny building #############


rule bcftools_mpileup_MM2_AsmToH37rv:
    input:
        MM2_AsmToH37rv_BAM = output_Dir + "/AsmAnalysis/{sampleID}/VariantCallingVersusH37Rv/MM2_AsmToH37rv/{sampleID}.mm2.AsmToH37Rv.bam",
        H37rv_FA = refGenome_FA_PATH,     
    output:
        MM2_AsmToH37rv_BAM_mpileup_out = output_Dir + "/AsmAnalysis/{sampleID}/VariantCallingVersusH37Rv/MM2_AsmToH37rv/{sampleID}.mm2.AsmToH37rv.bam.mpileup.txt.vcf",
        MM2_AsmToH37rv_BAM_mpileup_call_KeepAllPositions_VCF = output_Dir + "/AsmAnalysis/{sampleID}/VariantCallingVersusH37Rv/MM2_AsmToH37rv/{sampleID}.mm2.AsmToH37rv.bam.mpileup.call.KeepAllPositions.vcf"
    conda: "CondaEnvs/mm2_v2_4_WiUtilities.V2.yml"
    conda: "/home/mm774/conda3/envs/mm2_v2_4_WiUtilities/"
    threads: 1
    params:
        partition = 'short',
        slurm_log_dir = i_log_outdir,
    resources:
        mem_mb = lambda wildcards, attempt: (attempt * 1000) + 2000 , # Increment memory by 1000 mb per attempt 
        runtime_slurm = lambda wildcards, attempt: '0-00:' + str(attempt * 3) + ":00", #runtime = '0-00:10:00',
        runtime = lambda wildcards, attempt: (attempt * 3) , # minutes
    shell: # Run bcftools mpileup to summarize coverage and basepair outputs from the Minimap2 alignment (BAM)
        #"conda activate /home/mm774/conda3/envs/mm2_v2_4_WiUtilities \n"
        "bcftools mpileup -f {input.H37rv_FA} {input.MM2_AsmToH37rv_BAM} > {output.MM2_AsmToH37rv_BAM_mpileup_out} \n"
        "bcftools mpileup -f {input.H37rv_FA} {input.MM2_AsmToH37rv_BAM} | bcftools call -c -o {output.MM2_AsmToH37rv_BAM_mpileup_call_KeepAllPositions_VCF}"



##### Merge SNPs for Phylogeny generation (Based on Flye Asm) #####

rule getAll_SNPpositions_mpileup_MM2_AsmToH37rv:
    input:
        MM2_AsmToH37rv_BAM_mpileup_call_KeepAllPositions_VCF = output_Dir + "/AsmAnalysis/{sampleID}/VariantCallingVersusH37Rv/MM2_AsmToH37rv/{sampleID}.mm2.AsmToH37rv.bam.mpileup.call.KeepAllPositions.vcf",
    output:
        MM2_mpileup_VCF_AllSNPpositions_TSV = output_Dir + "/AsmAnalysis/{sampleID}/VariantCallingVersusH37Rv/MM2_AsmToH37rv/{sampleID}.mm2.AsmToH37rv.bam.mpileup.call.AllSNPpositions.tsv",
    conda: "CondaEnvs/mm2_v2_4_WiUtilities.V2.yml"
#    conda: "/home/mm774/conda3/envs/mm2_v2_4_WiUtilities/"
    threads: 1
    params:
        partition = 'short',
        slurm_log_dir = i_log_outdir,
    resources:
        mem_mb = lambda wildcards, attempt: (attempt * 1000) + 2000 , # Increment memory by 1000 mb per attempt 
        runtime_slurm = lambda wildcards, attempt: '0-00:' + str(attempt * 3) + ":00", #runtime = '0-00:10:00',
        runtime = lambda wildcards, attempt: (attempt * 3) , # minutes
    shell:
        #"conda activate /home/mm774/conda3/envs/mm2_v2_4_WiUtilities \n"
        "bcftools view --types snps {input.MM2_AsmToH37rv_BAM_mpileup_call_KeepAllPositions_VCF} | cut -f 1,2 | grep -v '#' > {output.MM2_mpileup_VCF_AllSNPpositions_TSV} \n"




rule combineAll_SNPpositions_mpileup_MM2_AsmToH37rv:
   input:
       expand(output_Dir + "/AsmAnalysis/{sampleID}/VariantCallingVersusH37Rv/MM2_AsmToH37rv/{sampleID}.mm2.AsmToH37rv.bam.mpileup.call.AllSNPpositions.tsv", sampleID = input_All_SampleIDs),           
   output:
       AllSample_AllSNPpositions_MM2_mpileup_TSV = output_Dir + "/Asm_MergeSNPs_mpileup/AllSNPpositions.MM2.mpileup.SNPs.Union.AllSamples.tsv",
   threads: 1
   params:
       partition = 'short',
       slurm_log_dir = i_log_outdir,
   resources:
       mem_mb = lambda wildcards, attempt: (attempt * 1000) + 2000 , # Increment memory by 1000 mb per attempt 
       runtime_slurm = lambda wildcards, attempt: '0-00:' + str(attempt * 3) + ":00", #runtime = '0-00:10:00',
       runtime = lambda wildcards, attempt: (attempt * 3) , # minutes
   shell:
       "cat {input} | sort -k 2n | uniq > {output.AllSample_AllSNPpositions_MM2_mpileup_TSV}"





rule Filter_MM2_mpileup_VarCalling_To_OnlySNPpositionsInUnionOfAllSNPs:
    input:
        MM2_AsmToH37rv_BAM_mpileup_call_KeepAllPositions_VCF = output_Dir + "/AsmAnalysis/{sampleID}/VariantCallingVersusH37Rv/MM2_AsmToH37rv/{sampleID}.mm2.AsmToH37rv.bam.mpileup.call.KeepAllPositions.vcf",
        AllSample_AllSNPpositions_MM2_mpileup_TSV = output_Dir + "/Asm_MergeSNPs_mpileup/AllSNPpositions.MM2.mpileup.SNPs.Union.AllSamples.tsv"
    output:
        MM2_AsmToRef_AllPositions_BCF_GZ = output_Dir + "/AsmAnalysis/{sampleID}/VariantCallingVersusH37Rv/MM2_AsmToH37rv/{sampleID}.mm2.AsmToH37rv.AllPositions.bcf.gz",
        MM2_AsmToRef_SNPs_BCF_GZ_SNPsInAllSamples = output_Dir + "/AsmAnalysis/{sampleID}/VariantCallingVersusH37Rv/MM2_AsmToH37rv/{sampleID}.mm2.AsmToH37rv.mpileup.call.SNPs.Union.AllSamples.bcf.gz",
        SampleID_TXT = output_Dir + "/AsmAnalysis/{sampleID}/VariantCallingVersusH37Rv/MM2_AsmToH37rv/{sampleID}.name.txt",
        MM2_AsmToRef_SNPs_BCF_GZ_SNPsInAllSamples_Renamed = output_Dir + "/AsmAnalysis/{sampleID}/VariantCallingVersusH37Rv/MM2_AsmToH37rv/{sampleID}.mm2.AsmToH37rv.mpileup.call.SNPs.Union.AllSamples.Renamed.bcf.gz",
        BCF_Renamed_PATH_TXT = output_Dir + "/AsmAnalysis/{sampleID}/VariantCallingVersusH37Rv/MM2_AsmToH37rv/{sampleID}.BCF_RenamedPATH.txt",
    conda: "CondaEnvs/mm2_v2_4_WiUtilities.V2.yml"
#    conda: "/home/mm774/conda3/envs/mm2_v2_4_WiUtilities/"
    threads: 1
    params:
        partition = 'short',
        slurm_log_dir = i_log_outdir,
    resources:
        mem_mb = lambda wildcards, attempt: (attempt * 1000) + 5000 , # Increment memory by 1000 mb per attempt 
        runtime_slurm = lambda wildcards, attempt: '0-00:' + str(attempt * 10) + ":00", #runtime = '0-00:10:00',
        runtime = lambda wildcards, attempt: (attempt * 10) , # minutes
    shell:
        "bcftools view {input.MM2_AsmToH37rv_BAM_mpileup_call_KeepAllPositions_VCF} -O b -o {output.MM2_AsmToRef_AllPositions_BCF_GZ} \n"
        
        "bcftools index {output.MM2_AsmToRef_AllPositions_BCF_GZ} \n"

        "bcftools view {output.MM2_AsmToRef_AllPositions_BCF_GZ} "
        " -R {input.AllSample_AllSNPpositions_MM2_mpileup_TSV} -e 'DP!=1' -O b -o {output.MM2_AsmToRef_SNPs_BCF_GZ_SNPsInAllSamples} \n"

        "bcftools index {output.MM2_AsmToRef_SNPs_BCF_GZ_SNPsInAllSamples} \n"

        "echo {wildcards.sampleID} > {output.SampleID_TXT} \n"
        
        "bcftools reheader -s {output.SampleID_TXT} {output.MM2_AsmToRef_SNPs_BCF_GZ_SNPsInAllSamples} -o {output.MM2_AsmToRef_SNPs_BCF_GZ_SNPsInAllSamples_Renamed} \n "
        
        "bcftools index {output.MM2_AsmToRef_SNPs_BCF_GZ_SNPsInAllSamples_Renamed} \n"
        
        "echo '{output.MM2_AsmToRef_SNPs_BCF_GZ_SNPsInAllSamples_Renamed}' > {output.BCF_Renamed_PATH_TXT} "





rule merge_BCFs_Renamed_PATHs_To_TXT:
    input:
        expand(output_Dir + "/AsmAnalysis/{sampleID}/VariantCallingVersusH37Rv/MM2_AsmToH37rv/{sampleID}.BCF_RenamedPATH.txt", sampleID = input_All_SampleIDs), 
    output:
        output_Dir + "/Asm_MergeSNPs_mpileup/ListOfAll_PATHs_BCFs_Renamed.txt"
    threads: 1
    params:
        partition = 'short',
        slurm_log_dir = i_log_outdir,
    resources:
        mem_mb = lambda wildcards, attempt: (attempt * 1000) + 5000 , # Increment memory by 1000 mb per attempt 
        runtime_slurm = lambda wildcards, attempt: '0-00:' + str(attempt * 5) + ":00", #runtime = '0-00:10:00',
        runtime = lambda wildcards, attempt: (attempt * 5) , # minutes
    shell:
        "cat {input} > {output}"






rule merge_All_BCFs_Renamed_To_VCF:
    input:
        output_Dir + "/Asm_MergeSNPs_mpileup/ListOfAll_PATHs_BCFs_Renamed.txt"
    output:
        Merged_VCF = output_Dir + "/Asm_MergeVar_mpileup/MM2.mpileup.call.Merged.vcf",
        Merged_VCF_SNVs = output_Dir + "/Asm_MergeVar_mpileup/MM2.mpileup.call.Merged.SNVs.vcf",
        Merged_VCF_SNVs_POS = output_Dir + "/Asm_MergeSNPs_mpileup/MM2.mpileup.call.MergedSNPs.snp.positions",
        Merged_VCF_SNVs_GZ = output_Dir + "/Asm_MergeVar_mpileup/MM2.mpileup.call.Merged.SNVs.vcf.gz",
    conda: "CondaEnvs/mm2_v2_4_WiUtilities.V2.yml"
#    conda: "/home/mm774/conda3/envs/mm2_v2_4_WiUtilities/"
    threads: 1
    params:
        partition = 'short',
        slurm_log_dir = i_log_outdir,
    resources:
        mem_mb = lambda wildcards, attempt: (attempt * 1000) + 5000 , # Increment memory by 1000 mb per attempt 
        runtime_slurm = lambda wildcards, attempt: '0-00:' + str(attempt * 5) + ":00", #runtime = '0-00:10:00',
        runtime = lambda wildcards, attempt: (attempt * 5) , # minutes
    shell:
        #"conda activate /home/mm774/conda3/envs/mm2_v2_4_WiUtilities \n"
        ' bcftools merge -i "-" -l {input} -o {output.Merged_VCF} -O v \n' # Does NOT fill AMB positions, ambigous calls stay ambigous
        ""
        " bcftools view {output.Merged_VCF} --types snps > {output.Merged_VCF_SNVs} \n"

        'grep -v "#" {output.Merged_VCF_SNVs} | cut -f 2  > {output.Merged_VCF_SNVs_POS} \n'
        ""
        " bgzip -c {output.Merged_VCF_SNVs} > {output.Merged_VCF_SNVs_GZ} \n"
        " tabix {output.Merged_VCF_SNVs_GZ} "



rule convert_MergedVCF_To_FASTA_ALN:
    input:
        Merged_VCF_SNVs = output_Dir + "/Asm_MergeVar_mpileup/MM2.mpileup.call.Merged.SNVs.vcf",
    output:
        MergedSNVs_FA = output_Dir + "/Asm_MergeVar_mpileup/MM2.mpileup.call.Merged.SNVs.min-100.fasta",
    #conda: "CondaEnvs/mm2_v2_4_WiUtilities.V2.yml"
    threads: 1
    params:
        min_Supporting = -100,
        partition = 'short',
        slurm_log_dir = i_log_outdir,
    resources:
        mem_mb = lambda wildcards, attempt: (attempt * 1000) + 5000 , # Increment memory by 1000 mb per attempt 
        runtime_slurm = lambda wildcards, attempt: '0-00:' + str(attempt * 5) + ":00", #runtime = '0-00:10:00',
        runtime = lambda wildcards, attempt: (attempt * 5) , # minutes
    shell:
        "Scripts/vcf2phylip/vcf2phylip.py -i {input} -f -m {params.min_Supporting} \n"


### Filtering of Merged-VCF ### 


### Filter by AMB threshold (Maximum % of AMB allowed a position)

rule filter_SNVs_10AmbThresh_MergeSNVs_mpileup:
    input:
        Merged_VCF_SNVs_GZ = output_Dir + "/Asm_MergeVar_mpileup/MM2.mpileup.call.Merged.SNVs.vcf.gz"
    output:
        Merged_VCF_SNVs_10AmbFilt = output_Dir + "/Asm_MergeVar_mpileup/MM2.mpileup.call.Merged.SNVs.10AmbThresh.NoMask.vcf",
        Merged_VCF_SNVs_10AmbFilt_GZ = output_Dir + "/Asm_MergeVar_mpileup/MM2.mpileup.call.Merged.SNVs.10AmbThresh.NoMask.vcf.gz",
        Merged_VCF_SNVs_10AmbFilt_POS = output_Dir + "/Asm_MergeVar_mpileup/MM2.mpileup.call.Merged.SNVs.10AmbThresh.NoMask.positions",
    conda: "CondaEnvs/mm2_v2_4_WiUtilities.V2.yml"
#    conda: "/home/mm774/conda3/envs/mm2_v2_4_WiUtilities/"
    threads: 1
    params:
        partition = 'short',
        slurm_log_dir = i_log_outdir,
    resources:
        mem_mb = lambda wildcards, attempt: (attempt * 1000) + 5000 , # Increment memory by 1000 mb per attempt 
        runtime_slurm = lambda wildcards, attempt: '0-00:' + str(attempt * 5) + ":00", #runtime = '0-00:10:00',
        runtime = lambda wildcards, attempt: (attempt * 5) , # minutes
    shell:
        #"conda activate /home/mm774/conda3/envs/mm2_v2_4_WiUtilities \n"
        'bcftools view {input.Merged_VCF_SNVs_GZ} --types snps -e "F_MISSING > 0.10"  > {output.Merged_VCF_SNVs_10AmbFilt} \n'

        'grep -v "#" {output.Merged_VCF_SNVs_10AmbFilt} | cut -f 2  > {output.Merged_VCF_SNVs_10AmbFilt_POS} \n'

        " bgzip -c {output.Merged_VCF_SNVs_10AmbFilt} > {output.Merged_VCF_SNVs_10AmbFilt_GZ} \n"
        " tabix {output.Merged_VCF_SNVs_10AmbFilt_GZ} "



rule convert_MergedVCF_10AmbThresh_To_FASTA_ALN:
    input:
        Merged_VCF_SNVs_10AmbFilt = output_Dir + "/Asm_MergeVar_mpileup/MM2.mpileup.call.Merged.SNVs.10AmbThresh.NoMask.vcf",
    output:
        MergedSNVs_10AmbFilt_FA = output_Dir + "/Asm_MergeVar_mpileup/MM2.mpileup.call.Merged.SNVs.10AmbThresh.NoMask.min-100.fasta",
    #conda: "CondaEnvs/mm2_v2_4_WiUtilities.V2.yml"
    threads: 1
    params:
        min_Supporting = -100,
        partition = 'short',
        slurm_log_dir = i_log_outdir,
    resources:
        mem_mb = lambda wildcards, attempt: (attempt * 1000) + 5000 , # Increment memory by 1000 mb per attempt 
        runtime_slurm = lambda wildcards, attempt: '0-00:' + str(attempt * 5) + ":00", #runtime = '0-00:10:00',
        runtime = lambda wildcards, attempt: (attempt * 5) , # minutes
    shell:
        "Scripts/vcf2phylip/vcf2phylip.py -i {input} -f -m {params.min_Supporting} \n"


########################################################################




rule filter_SNVs_RLC_Regions_MergeSNVs_mpileup:
    input:
        Merged_VCF_SNVs = output_Dir + "/Asm_MergeVar_mpileup/MM2.mpileup.call.Merged.SNVs.vcf",
        RLC_regions_BED = "References/Mtb_H37Rv_MaskingSchemes/RLC_Regions.H37Rv.bed"
    output:
        Merged_VCF_SNVs_RLCfilt = output_Dir + "/Asm_MergeVar_mpileup/MM2.mpileup.call.Merged.SNVs.RLCMask.vcf",
        Merged_VCF_SNVs_RLCfilt_GZ = output_Dir + "/Asm_MergeVar_mpileup/MM2.mpileup.call.Merged.SNVs.RLCMask.vcf.gz",
        Merged_VCF_SNVs_RLCfilt_POS = output_Dir + "/Asm_MergeVar_mpileup/MM2.mpileup.call.Merged.SNVs.RLCMask.positions",
    conda: "CondaEnvs/mm2_v2_4_WiUtilities.V2.yml"
#    conda: "/home/mm774/conda3/envs/mm2_v2_4_WiUtilities/"
    threads: 1
    params:
        partition = 'short',
        slurm_log_dir = i_log_outdir,
    resources:
        mem_mb = lambda wildcards, attempt: (attempt * 1000) + 2000 , # Increment memory by 1000 mb per attempt 
        runtime_slurm = lambda wildcards, attempt: '0-00:' + str(attempt * 2) + ":00", #runtime = '0-00:10:00',
        runtime = lambda wildcards, attempt: (attempt * 2), # minutes
    shell:
        #"conda activate /home/mm774/conda3/envs/mm2_v2_4_WiUtilities \n"
        'bedtools intersect -header -v -a {input.Merged_VCF_SNVs} -b {input.RLC_regions_BED} -wa > {output.Merged_VCF_SNVs_RLCfilt} \n'
        'grep -v "#" {output.Merged_VCF_SNVs_RLCfilt} | cut -f 2  > {output.Merged_VCF_SNVs_RLCfilt_POS} \n'
        ' bgzip -c {output.Merged_VCF_SNVs_RLCfilt} > {output.Merged_VCF_SNVs_RLCfilt_GZ} \n'
        ' tabix {output.Merged_VCF_SNVs_RLCfilt_GZ} '



rule filter_SNVs_PLC_Regions_MergeSNVs_mpileup:
    input:
        Merged_VCF_SNVs = output_Dir + "/Asm_MergeVar_mpileup/MM2.mpileup.call.Merged.SNVs.vcf",
        PLC_regions_BED = "References/Mtb_H37Rv_MaskingSchemes/PLC_Regions.H37Rv.bed"
    output:
        Merged_VCF_SNVs_PLCfilt = output_Dir + "/Asm_MergeVar_mpileup/MM2.mpileup.call.Merged.SNVs.PLCMask.vcf",
        Merged_VCF_SNVs_PLCfilt_GZ = output_Dir + "/Asm_MergeVar_mpileup/MM2.mpileup.call.Merged.SNVs.PLCMask.vcf.gz",
        Merged_VCF_SNVs_PLCfilt_POS = output_Dir + "/Asm_MergeVar_mpileup/MM2.mpileup.call.Merged.SNVs.PLCMask.positions",
    conda: "CondaEnvs/mm2_v2_4_WiUtilities.V2.yml"
    threads: 1
    params:
        partition = 'short',
        slurm_log_dir = i_log_outdir,
    resources:
        mem_mb = lambda wildcards, attempt: (attempt * 1000) + 2000 , # Increment memory by 1000 mb per attempt 
        runtime_slurm = lambda wildcards, attempt: '0-00:' + str(attempt * 2) + ":00", #runtime = '0-00:10:00',
        runtime = lambda wildcards, attempt: (attempt * 2) , # minutes
    shell:
        "bedtools intersect -header -v -a {input.Merged_VCF_SNVs} -b {input.PLC_regions_BED} -wa > {output.Merged_VCF_SNVs_PLCfilt} \n"
        'grep -v "#" {output.Merged_VCF_SNVs_PLCfilt} | cut -f 2  > {output.Merged_VCF_SNVs_PLCfilt_POS} \n'

        " bgzip -c {output.Merged_VCF_SNVs_PLCfilt} > {output.Merged_VCF_SNVs_PLCfilt_GZ} \n"
        " tabix {output.Merged_VCF_SNVs_PLCfilt_GZ} "


rule filter_SNVs_RLCandLowPMap_Regions_MergeSNVs_mpileup:
    input:
        Merged_VCF_SNVs_10AmbFilt = output_Dir + "/Asm_MergeVar_mpileup/MM2.mpileup.call.Merged.SNVs.10AmbThresh.NoMask.vcf",
        RLCandLowPmap_regions_BED = "References/Mtb_H37Rv_MaskingSchemes/RLC_Regions.Plus.LowPmapK50E4.H37Rv.bed"
    output:
        Merged_VCF_SNVs_RLCfilt = output_Dir + "/Asm_MergeVar_mpileup/MM2.mpileup.call.Merged.SNVs.10AmbThresh.RLCandLowPmapMask.vcf",
        Merged_VCF_SNVs_RLCfilt_GZ = output_Dir + "/Asm_MergeVar_mpileup/MM2.mpileup.call.Merged.SNVs.10AmbThresh.RLCandLowPmapMask.vcf.gz",
        Merged_VCF_SNVs_RLCfilt_POS = output_Dir + "/Asm_MergeVar_mpileup/MM2.mpileup.call.Merged.SNVs.10AmbThresh.RLCandLowPmapMask.positions",
    conda: "CondaEnvs/mm2_v2_4_WiUtilities.V2.yml"
    threads: 1
    params:
        partition = 'short',
        slurm_log_dir = i_log_outdir,
    resources:
        mem_mb = lambda wildcards, attempt: (attempt * 1000) + 2000 , # Increment memory by 1000 mb per attempt 
        runtime_slurm = lambda wildcards, attempt: '0-00:' + str(attempt * 2) + ":00", #runtime = '0-00:10:00',
        runtime = lambda wildcards, attempt: (attempt * 2) , # minutes
    shell:
        "bedtools intersect -header -v -a {input.Merged_VCF_SNVs_10AmbFilt} -b {input.RLCandLowPmap_regions_BED} -wa > {output.Merged_VCF_SNVs_RLCfilt} \n"
        'grep -v "#" {output.Merged_VCF_SNVs_RLCfilt} | cut -f 2  > {output.Merged_VCF_SNVs_RLCfilt_POS} \n'

        " bgzip -c {output.Merged_VCF_SNVs_RLCfilt} > {output.Merged_VCF_SNVs_RLCfilt_GZ} \n"
        " tabix {output.Merged_VCF_SNVs_RLCfilt_GZ} "




rule convert_MergedVCF_To_FASTA_ALN_RLC_removed:
    input:
        Merged_VCF_SNVs_RLCfilt = output_Dir + "/Asm_MergeVar_mpileup/MM2.mpileup.call.Merged.SNVs.RLCMask.vcf",
    output:
        MergedSNVs_RLCfilt_FA = output_Dir + "/Asm_MergeVar_mpileup/MM2.mpileup.call.Merged.SNVs.RLCMask.min-100.fasta",
    threads: 1
    params:
        min_Supporting = -100,
        partition = 'short',
        slurm_log_dir = i_log_outdir,
    resources:
        mem_mb = lambda wildcards, attempt: (attempt * 1000) + 5000 , # Increment memory by 1000 mb per attempt 
        runtime_slurm = lambda wildcards, attempt: '0-00:' + str(attempt * 5) + ":00", #runtime = '0-00:10:00',
        runtime = lambda wildcards, attempt: (attempt * 5) , # minutes
    shell:
        "Scripts/vcf2phylip/vcf2phylip.py -i {input} -f -m {params.min_Supporting} \n"



rule convert_MergedVCF_To_FASTA_ALN_PLC_removed:
    input:
        Merged_VCF_SNVs_RLCfilt = output_Dir + "/Asm_MergeVar_mpileup/MM2.mpileup.call.Merged.SNVs.PLCMask.vcf",
    output:
        MergedSNVs_RLCfilt_FA = output_Dir + "/Asm_MergeVar_mpileup/MM2.mpileup.call.Merged.SNVs.PLCMask.min-100.fasta",
    conda: "CondaEnvs/mm2_v2_4_WiUtilities.V2.yml"
    threads: 1
    params:
        min_Supporting = -100,
        partition = 'short',
        slurm_log_dir = i_log_outdir,
    resources:
        mem_mb = lambda wildcards, attempt: (attempt * 1000) + 5000 , # Increment memory by 1000 mb per attempt 
        runtime_slurm = lambda wildcards, attempt: '0-00:' + str(attempt * 5) + ":00", #runtime = '0-00:10:00',
        runtime = lambda wildcards, attempt: (attempt * 5) , # minutes
    shell:
        "Scripts/vcf2phylip/vcf2phylip.py -i {input} -f -m {params.min_Supporting} \n"


rule convert_MergedVCF_To_FASTA_ALN_RLCandLowPmap_removed:
    input:
        Merged_VCF_SNVs_RLCfilt = output_Dir + "/Asm_MergeVar_mpileup/MM2.mpileup.call.Merged.SNVs.10AmbThresh.RLCandLowPmapMask.vcf",
    output:
        MergedSNVs_RLCandPmap_filt_FA = output_Dir + "/Asm_MergeVar_mpileup/MM2.mpileup.call.Merged.SNVs.10AmbThresh.RLCandLowPmapMask.min-100.fasta",
    conda:"CondaEnvs/mm2_v2_4_WiUtilities.V2.yml"
    threads: 1
    params:
        min_Supporting = -100,
        partition = 'short',
        slurm_log_dir = i_log_outdir,
    resources:
        mem_mb = lambda wildcards, attempt: (attempt * 1000) + 5000 , # Increment memory by 1000 mb per attempt 
        runtime_slurm = lambda wildcards, attempt: '0-00:' + str(attempt * 5) + ":00", #runtime = '0-00:10:00',
        runtime = lambda wildcards, attempt: (attempt * 5) , # minutes
    shell:
        "Scripts/vcf2phylip/vcf2phylip.py -i {input} -f -m {params.min_Supporting} \n"







######## WGA - Nucleotide Diversity Analysis ########

rule calculate_NucDiversity_MergeSNVs_mpileup:
    input:
        Merged_VCF_SNVs = output_Dir + "/Asm_MergeVar_mpileup/MM2.mpileup.call.Merged.SNVs.vcf",
    output:
        output_Dir + "/NucDiversity/NucDiv_SNVs_mpileup/MM2.mpileup.call.Merged.SNVs.NucDiv.{NucDiv_WindowSize}bp.windowed.pi"
    conda: "CondaEnvs/mm2_v2_4_WiUtilities.V2.yml"
    threads: 1
    params:
        NucDiv_Output_Prefix = output_Dir + "/NucDiversity/NucDiv_SNVs_mpileup/MM2.mpileup.call.Merged.SNVs.NucDiv.{NucDiv_WindowSize}bp",
        partition = 'short',
        slurm_log_dir = i_log_outdir,
    resources:
        mem_mb = lambda wildcards, attempt: (attempt * 1000) + 5000 , # Increment memory by 1000 mb per attempt 
        runtime_slurm = lambda wildcards, attempt: '0-00:' + str(attempt * 5) + ":00", #runtime = '0-00:10:00',
        runtime = lambda wildcards, attempt: (attempt * 5) , # minutes
    shell:
        "vcftools --vcf {input} --window-pi {wildcards.NucDiv_WindowSize} --out {params.NucDiv_Output_Prefix}"

#####################################################




######## Phylogeny Building ########

rule fasttree_GTR_from_MergedSNPs:
    input:
        MergedSNVs_FA = output_Dir + "/Asm_MergeVar_mpileup/MM2.mpileup.call.Merged.SNVs.min-100.fasta",
    output:
        output_Dir + "/Phylogenies/fasttree_mpileupSNVs_NoFilt/MM2.mpileup.call.Merged.SNVs.min-100.fasttree.newick"   
    conda:
        "CondaEnvs/Gubbins_v3_2_1.yml"
    threads: 1
    params:
        partition = 'short',
        slurm_log_dir = i_log_outdir,
    resources:
        mem_mb = lambda wildcards, attempt: (attempt * 1000) + 8000 , # Increment memory by 1000 mb per attempt 
        runtime_slurm = lambda wildcards, attempt: '0-00:' + str(attempt * 30) + ":00", #runtime = '0-00:10:00',
        runtime = lambda wildcards, attempt: (attempt * 30) , # minutes
    shell:
        "time FastTree -nt -gtr {input} > {output}"


rule fasttree_GTR_from_MergedSNVs_10AmbFilt:
    input:
        MergedSNVs_10AmbFilt_FA = output_Dir + "/Asm_MergeVar_mpileup/MM2.mpileup.call.Merged.SNVs.10AmbThresh.NoMask.min-100.fasta",
    output:
        output_Dir + "/Phylogenies/fasttree_mpileupSNVs_10AmbFilt/MM2.mpileup.call.Merged.SNVs.10AmbThresh.NoMask.min-100.fasttree.newick",
    conda:
        "CondaEnvs/Gubbins_v3_2_1.yml"
    threads: 1
    params:
        partition = 'short',
        slurm_log_dir = i_log_outdir,
    resources:
        mem_mb = lambda wildcards, attempt: (attempt * 1000) + 8000 , # Increment memory by 1000 mb per attempt 
        runtime_slurm = lambda wildcards, attempt: '0-00:' + str(attempt * 30) + ":00", #runtime = '0-00:10:00',
        runtime = lambda wildcards, attempt: (attempt * 30) , # minutes
    shell:
        "time FastTree -nt -gtr {input} > {output}"




rule IQtree_GTR_from_MergedSNVs_10AmbFilt:
    input:
        MergedSNVs_10AmbFilt_FA = output_Dir + "/Asm_MergeVar_mpileup/MM2.mpileup.call.Merged.SNVs.10AmbThresh.NoMask.min-100.fasta",
    output:
        MergedSNVs_10AmbFilt_SNPSites_FA = output_Dir + "/Asm_MergeVar_mpileup/MM2.mpileup.call.Merged.SNVs.10AmbThresh.NoMask.min-100.VarSites.fasta",
        TREE = output_Dir + "/Phylogenies/iqtree_mpileupSNVs_10AmbFilt/MM2.mpileup.call.Merged.SNVs.10AmbThresh.NoMask.min-100.iq.treefile",
    threads: 1
    params:
        iqtree_outprefix = output_Dir + "/Phylogenies/iqtree_mpileupSNVs_10AmbFilt/MM2.mpileup.call.Merged.SNVs.10AmbThresh.NoMask.min-100.iq",
        bb = 10000,
        partition = 'short',
        slurm_log_dir = i_log_outdir,
    resources:
        mem_mb = lambda wildcards, attempt: (attempt * 1000) + 8000 , # Increment memory by 1000 mb per attempt 
        runtime_slurm = lambda wildcards, attempt: '0-00:' + str(attempt * 30) + ":00", #runtime = '0-00:10:00',
        runtime = lambda wildcards, attempt: (attempt * 30) , # minutes
    shell:
        "snp-sites {input.MergedSNVs_10AmbFilt_FA} > {output.MergedSNVs_10AmbFilt_SNPSites_FA} \n"
        "time iqtree -m GTR+ASC -nt 1 -bb {params.bb} -s {output.MergedSNVs_10AmbFilt_SNPSites_FA} -pre {params.iqtree_outprefix} " # -m # -redo


rule IQtree_GTR_from_MergedSNVs_PLCFilt:
    input:
        MergedSNVs_PLCfilt_FA = output_Dir + "/Asm_MergeVar_mpileup/MM2.mpileup.call.Merged.SNVs.PLCMask.min-100.fasta",
    output:
        MergedSNVs_PLCfilt_SNPSites_FA = output_Dir + "/Asm_MergeVar_mpileup/MM2.mpileup.call.Merged.SNVs.PLCMask.min-100.VarSites.fasta",
        TREE = output_Dir + "/Phylogenies/iqtree_mpileupSNVs_PLCFilt/MM2.mpileup.call.Merged.SNVs.PLCMask.min-100.iq.treefile",
    threads: 1
    params:
        iqtree_outprefix = output_Dir + "/Phylogenies/iqtree_mpileupSNVs_10AmbFilt/MM2.mpileup.call.Merged.SNVs.10AmbThresh.NoMask.min-100.iq",
        bb = 10000,
        partition = 'short',
        slurm_log_dir = i_log_outdir,
    resources:
        mem_mb = lambda wildcards, attempt: (attempt * 1000) + 8000 , # Increment memory by 1000 mb per attempt 
        runtime_slurm = lambda wildcards, attempt: '0-00:' + str(attempt * 30) + ":00", #runtime = '0-00:10:00',
        runtime = lambda wildcards, attempt: (attempt * 30) , # minutes
    shell:
        "snp-sites {input.MergedSNVs_PLCfilt_FA} > {output.MergedSNVs_PLCfilt_SNPSites_FA} \n"
        "time iqtree -m GTR+ASC -nt 1 -bb {params.bb} -s {output.MergedSNVs_PLCfilt_SNPSites_FA} -pre {params.iqtree_outprefix} " # -m # -redo





rule IQtree_GTR_from_MergedSNVs_10AmbFilt_RLCandLowPmap_Filt:
    input:
        MergedSNVs_RLCandPmap_filt_FA = output_Dir + "/Asm_MergeVar_mpileup/MM2.mpileup.call.Merged.SNVs.10AmbThresh.RLCandLowPmapMask.min-100.fasta",
    output:
        MergedSNVs_RLCandPmap_filt_SNPSites_FA = output_Dir + "/Asm_MergeVar_mpileup/MM2.mpileup.call.Merged.SNVs.10AmbThresh.RLCandLowPmapMask.min-100.VarSites.fasta",
        TREE = output_Dir + "/Phylogenies/iqtree_mpileupSNVs_RLCandLowPmapMask/MM2.mpileup.call.Merged.SNVs.10AmbThresh.RLCandLowPmapMask.min-100.iq.treefile",
    threads: 1
    params:
        iqtree_outprefix = output_Dir + "/Phylogenies/iqtree_mpileupSNVs_10AmbFilt/MM2.mpileup.call.Merged.SNVs.10AmbThresh.NoMask.min-100.iq",
        bb = 10000,
        partition = 'short',
        slurm_log_dir = i_log_outdir,
    resources:
        mem_mb = lambda wildcards, attempt: (attempt * 1000) + 8000 , # Increment memory by 1000 mb per attempt 
        runtime_slurm = lambda wildcards, attempt: '0-00:' + str(attempt * 30) + ":00", #runtime = '0-00:10:00',
        runtime = lambda wildcards, attempt: (attempt * 30) , # minutes
    shell:
        "snp-sites {input.MergedSNVs_RLCandPmap_filt_FA} > {output.MergedSNVs_RLCandPmap_filt_SNPSites_FA} \n"
        "time iqtree -m GTR+ASC -nt 1 -bb {params.bb} -s {output.MergedSNVs_RLCandPmap_filt_SNPSites_FA} -pre {params.iqtree_outprefix} " # -m # -redo



####################################






######## Gubbins Recombination Analysis ########


### Step 1: Create an full genome MSA w/ SNVs across all input genome (using H37Rv reference coordinates)

rule insert_SNVs_IntoRefGenome_PerSample_10AmbThresh:
    input:
        Merged_VCF_SNVs_10AmbFilt_GZ = output_Dir + "/Asm_MergeVar_mpileup/MM2.mpileup.call.Merged.SNVs.10AmbThresh.NoMask.vcf.gz",
        Merged_VCF_SNVs_GZ = output_Dir + "/Asm_MergeVar_mpileup/MM2.mpileup.call.Merged.SNVs.vcf.gz",
        H37rv_FA = refGenome_FA_PATH,
    output:
        inferredGenome_FromSNVs_FA = temp(output_Dir + "/RecombDetection/PopulationMSA_mpileup_SNVs_10AmbThresh/Indiv_ConsensusSeq_SNVs/{sampleID}.InferredGenome.WiSNVs.10AmbThresh.fasta"),
        inferredGenome_FromSNVs_Renamed_FA = output_Dir + "/RecombDetection/PopulationMSA_mpileup_SNVs_10AmbThresh/Indiv_ConsensusSeq_SNVs/{sampleID}.InferredGenome.WiSNVs.10AmbThresh.renamed.fasta",
    conda: "CondaEnvs/mm2_v2_4_WiUtilities.V2.yml"
    threads: 1
    params:
        partition = 'short',
        slurm_log_dir = i_log_outdir,
    resources:
        mem_mb = lambda wildcards, attempt: (attempt * 1000) + 2000 , # Increment memory by 1000 mb per attempt 
        runtime_slurm = lambda wildcards, attempt: '0-00:' + str(attempt * 3) + ":00", #runtime = '0-00:10:00',
        runtime = lambda wildcards, attempt: (attempt * 3) , # minutes
    shell:
        "bcftools consensus -f {input.H37rv_FA} --missing '-' -o {output.inferredGenome_FromSNVs_FA} -s {wildcards.sampleID} {input.Merged_VCF_SNVs_10AmbFilt_GZ} \n"   
        ""
        " bioawk -c fastx '{{ print \">{wildcards.sampleID}\" \"\\n\" $seq }}' {output.inferredGenome_FromSNVs_FA} > {output.inferredGenome_FromSNVs_Renamed_FA}"


rule create_FullGenome_MSA_SNVs_10AmbThresh:
    input:
        expand(output_Dir + "/RecombDetection/PopulationMSA_mpileup_SNVs_10AmbThresh/Indiv_ConsensusSeq_SNVs/{sampleID}.InferredGenome.WiSNVs.10AmbThresh.renamed.fasta", sampleID = input_All_SampleIDs), 
    output:
        Inferred_FullGenomeMSA_FA = output_Dir + "/RecombDetection/PopulationMSA_mpileup_SNVs_10AmbThresh/InferredFullGenomeMSA.FromSNVs.10AmbThresh.fasta"
    threads: 1
    params:
        partition = 'short',
        slurm_log_dir = i_log_outdir,
    resources:
        mem_mb = lambda wildcards, attempt: (attempt * 1000) + 3000 , # Increment memory by 1000 mb per attempt 
        runtime_slurm = lambda wildcards, attempt: '0-00:' + str(attempt * 5) + ":00", #runtime = '0-00:10:00',
        runtime = lambda wildcards, attempt: (attempt * 5) , # minutes
    shell:
        "cat {input} > {output.Inferred_FullGenomeMSA_FA}"




### Step 2: Run Gubbins on the constructed MSA with inserted SNPs

rule run_Gubbins_v321_ExtensiveSearch_SmallWin_MS4_Wi_mpileup_SNVs_10AmbThresh:
    input:
        Inferred_FullGenomeMSA_FA = output_Dir + "/RecombDetection/PopulationMSA_mpileup_SNVs_10AmbThresh/InferredFullGenomeMSA.FromSNVs.10AmbThresh.fasta",
        phylogeny_mpileup_SNVs_10Amb = output_Dir + "/Phylogenies/fasttree_mpileupSNVs_10AmbFilt/MM2.mpileup.call.Merged.SNVs.10AmbThresh.NoMask.min-100.fasttree.newick",
    output:
        NodeLabelled_Tree = output_Dir + "/RecombDetection/Gubbins_v321_ExtSearch_SW_MS4_mpileup_SNVs_10AmbThresh/Gubbins.node_labelled.final_tree.tre",
        BranchStats_CSV = output_Dir + "/RecombDetection/Gubbins_v321_ExtSearch_SW_MS4_mpileup_SNVs_10AmbThresh/Gubbins.per_branch_statistics.csv",
        RecombPreds_GFF = output_Dir + "/RecombDetection/Gubbins_v321_ExtSearch_SW_MS4_mpileup_SNVs_10AmbThresh/Gubbins.recombination_predictions.gff",
    conda: "CondaEnvs/Gubbins_v3_2_1.yml"
    threads: 1
    params:
        Gubbins_Output_Prefix_PATH = output_Dir + "/RecombDetection/Gubbins_v321_ExtSearch_SW_MS4_mpileup_SNVs_10AmbThresh/Gubbins",
        Gubbins_OutputDir_PATH     = output_Dir + "/RecombDetection/Gubbins_v321_ExtSearch_SW_MS4_mpileup_SNVs_10AmbThresh/",
        min_SNPs = 4, # Min SNPs to identify a recombination block (default: 3)
        minWinSize = 25,
        maxWinSize = 1000,
        partition = 'short',
        slurm_log_dir = i_log_outdir,
    resources:
        mem_mb = lambda wildcards, attempt: (attempt * 1000) + 3000 , # Increment memory by 1000 mb per attempt 
        runtime_slurm = lambda wildcards, attempt: '0-' + str(attempt * 2) + ":00:00", #runtime = '0-00:10:00',
        runtime = lambda wildcards, attempt: (attempt * 120) , # minutes
    shell:
        #" mkdir {params.Gubbins_OutputDir_PATH}/ \n"
        " cd {params.Gubbins_OutputDir_PATH} \n"
        "time run_gubbins.py --extensive-search "
        " --min-window-size {params.minWinSize} --max-window-size {params.maxWinSize} --min-snps {params.min_SNPs} "
        " --starting-tree {input.phylogeny_mpileup_SNVs_10Amb} "
        " --prefix {params.Gubbins_Output_Prefix_PATH} {input.Inferred_FullGenomeMSA_FA} "


SED_Update_Chr_To_H37Rv = 's/SEQUENCE/NC_000962.3/g'

rule reformat_GubbinsEvents_v321_ExtensiveSearch_SmallWin_MS4_Wi_SNV_10AmbThresh:
    input:
        RecombPreds_GFF = output_Dir + "/RecombDetection/Gubbins_v321_ExtSearch_SW_MS4_mpileup_SNVs_10AmbThresh/Gubbins.recombination_predictions.gff",
    output:
        RecombPreds_RenamedChr_GFF = output_Dir + "/RecombDetection/Gubbins_v321_ExtSearch_SW_MS4_mpileup_SNVs_10AmbThresh/Gubbins.recombination_predictions.RenamedCHR.gff",
        RecombPreds_RenamedChr_BED = output_Dir + "/RecombDetection/Gubbins_v321_ExtSearch_SW_MS4_mpileup_SNVs_10AmbThresh/Gubbins.recombination_predictions.RenamedCHR.bed",
    conda: "CondaEnvs/mm2_v2_4_WiUtilities.V2.yml"
    threads: 1
    params:
        partition = 'short',
        slurm_log_dir = i_log_outdir,
    resources:
        mem_mb = lambda wildcards, attempt: (attempt * 1000) + 3000 , # Increment memory by 1000 mb per attempt 
        runtime_slurm = lambda wildcards, attempt: '0-00:' + str(attempt * 5) + ":00", #runtime = '0-00:10:00',
        runtime = lambda wildcards, attempt: (attempt * 5) , # minutes
    shell:
        "sed {SED_Update_Chr_To_H37Rv} {input.RecombPreds_GFF} | gff2bed > {output.RecombPreds_RenamedChr_BED} \n"

        "sed {SED_Update_Chr_To_H37Rv} {input.RecombPreds_GFF} > {output.RecombPreds_RenamedChr_GFF} "










###############################################
############ Lineage Calling ##################
###############################################


rule FastLinCaller_AsmToH37Rv_MM2:
    input:
        MM2_AsmToH37rv_paftools_VCF = output_Dir + "/AsmAnalysis/{sampleID}/VariantCallingVersusH37Rv/MM2_AsmToH37rv/{sampleID}.mm2.AsmToH37Rv.paftools.vcf",
    output:
        output_Dir + "/AsmAnalysis/{sampleID}/LineageCalling/LinCall_Paftools_AsmToH37Rv/{sampleID}.AsmToH37Rv.lineage_call.tsv"
    conda:
        "CondaEnvs/fastlincaller_v032.yml"
    params:
        partition = 'short',
        slurm_log_dir = i_log_outdir,
    resources:
        mem_mb = 3000,
        runtime_slurm = lambda wildcards, attempt: '0-00:' + str(attempt * 5) + ":00", #runtime = '0-00:10:00',
        runtime = lambda wildcards, attempt: (attempt * 5) , # minutes
    shell:
        "fast-lineage-caller {input} --out {output}"

###############################################







###############################################
############ Homlogy Mapping Analysis #########
###############################################

refGenome_FA_PATH = config["RefGenome_FA_PATH"]
refGenome_GFF_PATH = config["RefGenome_GFF_PATH"]



###############################################


###############################################
############ Homlogy Mapping Analysis #########
###############################################


AWK_STR_PAF_HmMap_Update = '{if ( ($0 ~ "s2:i:0") && ($0 ~ "zd:i:2") ) print $0; else print $0"\\ts2:i:0\\tzd:i:2"}'

SED_STR_Update_SecToPri_Aln = 's/tp:A:S/tp:A:P/g'

AWK_STR_FiltBy_RefLen = ' $4 - $3 <= 10000 ' # 4 and 5th columns signify the target (reference) coordinates of the PAF alignment


rule HomologyMapping_Ref_k19w19Param:
    input:
        H37rv_FA = refGenome_FA_PATH,
    output:
        H37Rv_Hmap_SAM = temp(output_Dir + "/HomologyMapping/H37Rv_Ref_HomologyMapping_k19w19Param/H37Rv.HomologyMap.sam"),
        H37Rv_Hmap_Viz_SAM = output_Dir + "/HomologyMapping/H37Rv_Ref_HomologyMapping_k19w19Param/H37Rv.HomologyMap.VizUpdated.sam",
        H37Rv_Hmap_Viz_BAM = output_Dir + "/HomologyMapping/H37Rv_Ref_HomologyMapping_k19w19Param/H37Rv.HomologyMap.VizUpdated.bam",
        H37Rv_Hmap_Viz_BAI = output_Dir + "/HomologyMapping/H37Rv_Ref_HomologyMapping_k19w19Param/H37Rv.HomologyMap.VizUpdated.bam.bai",
        H37Rv_Hmap_PAF = output_Dir + "/HomologyMapping/H37Rv_Ref_HomologyMapping_k19w19Param/H37Rv.HomologyMap.paf",
        H37Rv_Hmap_TagMod_PAF = output_Dir + "/HomologyMapping/H37Rv_Ref_HomologyMapping_k19w19Param/H37Rv.HomologyMap.TagMod.paf",
        H37Rv_Hmap_TagMod_VarTSV = output_Dir + "/HomologyMapping/H37Rv_Ref_HomologyMapping_k19w19Param/H37Rv.HomologyMap.TagMod.var.tsv",
        H37Rv_Hmap_TagMod_noR_VarTSV = output_Dir + "/HomologyMapping/H37Rv_Ref_HomologyMapping_k19w19Param/H37Rv.HomologyMap.TagMod.noR.var.tsv",
        H37Rv_Hmap_PAF_Trimmed = output_Dir + "/HomologyMapping/H37Rv_Ref_HomologyMapping_k19w19Param/H37Rv.HomologyMap.Col1to12.paf",
        H37Rv_Hmap_MergedRegions_BED = output_Dir + "/HomologyMapping/H37Rv_Ref_HomologyMapping_k19w19Param/H37Rv.HomologyMap.merged.WiCounts.bed"
    conda: "CondaEnvs/mm2_v2_4_WiUtilities.V2.yml"
    params:
        scripts_dir = "./Scripts",
        mm2_aln_params = "-k19 -w19", # Set to "" if you do not want to set any parameters
        partition = 'short',
        slurm_log_dir = i_log_outdir,
    resources:
        mem_mb = 3000,
        runtime_slurm = lambda wildcards, attempt: '0-00:' + str(attempt * 5) + ":00", #runtime = '0-00:10:00',
        runtime = lambda wildcards, attempt: (attempt * 5) , # minutes
    shell:
        "minimap2 --MD -DP {params.mm2_aln_params} -a --cs {input.H37rv_FA} {input.H37rv_FA} > {output.H37Rv_Hmap_SAM} \n"

        "{params.scripts_dir}/updateSAM_ForMM2SelfAlign.py --input_sam {output.H37Rv_Hmap_SAM} "
        "--aligned_sequence_fasta {input.H37rv_FA} --output_sam {output.H37Rv_Hmap_Viz_SAM} \n"

        "samtools view -bS {output.H37Rv_Hmap_Viz_SAM} | samtools sort - > {output.H37Rv_Hmap_Viz_BAM} \n"
        "samtools index {output.H37Rv_Hmap_Viz_BAM} \n"

        "minimap2 -DP {params.mm2_aln_params} -c --cs {input.H37rv_FA} {input.H37rv_FA} > {output.H37Rv_Hmap_PAF} \n"

        "awk '{AWK_STR_PAF_HmMap_Update}' {output.H37Rv_Hmap_PAF} | sed '{SED_STR_Update_SecToPri_Aln}' | sort -k3,3n -k4,4n > {output.H37Rv_Hmap_TagMod_PAF} \n"

        "paftools.js call -q 0 -l 100 -L 100 {output.H37Rv_Hmap_TagMod_PAF} | sort -k3n -k4n > {output.H37Rv_Hmap_TagMod_VarTSV}  \n"

        "cut -f 1-12 {output.H37Rv_Hmap_PAF} > {output.H37Rv_Hmap_PAF_Trimmed} \n"

        "bedtools bamtobed -i {output.H37Rv_Hmap_SAM} | cut -f 1-3 | sort -k2,2n -k3,3n | bedtools merge -c 1 -o count > {output.H37Rv_Hmap_MergedRegions_BED} \n"

        "grep -v ^'R' {output.H37Rv_Hmap_TagMod_VarTSV} > {output.H37Rv_Hmap_TagMod_noR_VarTSV} "





