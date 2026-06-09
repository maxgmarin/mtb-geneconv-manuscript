# pgqc/module1.py

import pandas as pd
import numpy as np

import bioframe as bf

import io 
import pathlib

from tqdm import tqdm

import re



## Compute genomic distance - have to figure out whether positive or negative direction is closer on the genome
def Rv_dist(p1,p2):
    """
    computes distance accounting for circular genome
    valid for MTB H37Rv only (relies on knowing total number of nucleotides)
    
    p1, p2: int
        positions 1 and 2 on the genome
    """
    # accounts for circularization of genome
    
    total_nucleotides = 4411532    

    max_dist = int(total_nucleotides / 2)
    
    d1 = np.abs(p1 - p2)

    #  distance cannot be farther than the maximum distance
    if d1 > max_dist:
        # sort so that "i" is the smaller of the two positions
        i,j=sorted([p1,p2])        
        # add the genomic distance to i
        i = i + total_nucleotides
        d1 = i-j
        
    return d1




# regex for cs tag parsing
CS_OP_RE = re.compile(r'([:=*+-])(\d+|[A-Za-z]+)')

def count_substitutions_in_cs(cs: str) -> int:
    """
    Count the number of substitution events (*xy) in a minimap2 cs tag.

    Parameters
    ----------
    cs : str
        A single cs tag string (e.g. ':6-ata:10+gtc:4*at:3').

    Returns
    -------
    int
        Number of substitutions detected in this cs string.
    """
    if cs is None or (isinstance(cs, float) and pd.isna(cs)):
        return 0

    n_sub = 0
    for op, arg in CS_OP_RE.findall(cs):
        if op == '*':
            # arg is two bases: ref, query
            ref_base = arg[0]
            alt_base = arg[1]
            # mirror your variant logic: skip ambiguous Ns
            if ref_base.lower() != "n" and alt_base.lower() != "n":
                n_sub += 1

    return n_sub



def calculate_gc_content_perwindow(genome_seq, df,
                                   start_col = "Target_Start",
                                   end_col   = "Target_End"):
    """
    Compute GC content for each genomic window in the DataFrame.
    
    Parameters:
    - genome (str): The full genome sequence as a string.
    - df (pd.DataFrame): DataFrame containing "Target_Start" and "Target_End" columns.
    
    Returns:
    - pd.DataFrame: Updated DataFrame with a new column "GC_Content".
    """
    gc_contents = []
    
    for _, row in df.iterrows():
        start = int(row[start_col]) - 1  # Convert to 0-based indexing
        end = int(row[end_col])
        
        # Extract sequence from genome string
        subseq = genome_seq[start:end]
        
        # Calculate GC content
        gc_count = subseq.count('G') + subseq.count('C')
        total_bases = len(subseq)
        gc_content = (gc_count / total_bases) * 100 if total_bases > 0 else None
        
        gc_contents.append(gc_content)

    # Add results to DataFrame
    df["GC_Content"] = gc_contents
    df["GC_Content"] = df["GC_Content"].astype(float)

    return df




def compute_Paralog_EdgeToEdge_Distance(row):
    """
    Given a Pandas row with Target_Start, Target_End, Query_Start, Query_End,
    compute the minimal circular distance between any pair of edges
    using Rv_dist().
    """
    t_start = row["Target_Start"]
    t_end   = row["Target_End"]
    q_start = row["Query_Start"]
    q_end   = row["Query_End"]
    
    d1 = Rv_dist(t_start, q_start)
    d2 = Rv_dist(t_start, q_end)
    d3 = Rv_dist(t_end,   q_start)
    d4 = Rv_dist(t_end,   q_end)
    
    return min(d1, d2, d3, d4)



def label_DF_ByOvrLapGenes(i_DF,
                           i_GenomeAnno_Genes_DF,
                           GeneNameColumn = "Symbol",
                           i_cols1 = ["Chrom", "Start", "End"],
                           i_cols2 = ["Chrom", "Start", "End"],):
    
    listOf_Overlap_Genes = []

    for i, row in i_DF.iterrows():
        
        # a) Target overlapping genes
        i_Chrom = row[ i_cols1[0] ]
        i_Start, i_End = int(row[i_cols1[1]]), int(row[i_cols1[2]])

        Target_Range = f"{i_Chrom}:{i_Start}-{i_End}"

        sub_DF_Overlap_Genes = bf.select(i_GenomeAnno_Genes_DF, Target_Range, cols = i_cols2)

        listOf_Overlap_Genes.append( ",".join(list(sub_DF_Overlap_Genes[ GeneNameColumn ].values)) )


    i_Anno_DF = i_DF
    i_Anno_DF["Overlap_Genes"] = listOf_Overlap_Genes

    # For all cases where there are NO OVERLAPPING GENES, simply put "_"
    i_Anno_DF["Overlap_Genes"] = i_Anno_DF["Overlap_Genes"].replace("", np.nan).fillna("_")

    return i_Anno_DF


#### Variant Parsing and Annotation (For minimap2 + paftools.js) ####


def infer_NT_var_type(row):
    
    """Return 'SNP'/'INS'/'DEL'/'INDEL'/'UNKNOWN' from available columns."""

    ref = row["Ref"].replace("-", "")
    alt = row["Alt"].replace("-", "")
    
    if row["Ref"] and row["Alt"]:
        if len(ref) == len(alt):
            if len(ref) == 1:
                return "SNP"
            if len(ref) > 1:
                return "MNP"
                
        if len(alt) > len(ref): return "INS"
        if len(alt) < len(ref): return "DEL"

    return "UNKNOWN"


# From the paftools.js source, the fields after V are:

# [0] Ref_Name
# [1] Ref_Start
# [2] Ref_End
# [3] Cov
# [4] MapQ
# [5] Ref_Allele
# [6] Alt_Allele
# [7] Query_Name
# [8] Query_Start
# [9] Query_End
# [10] Strand

def parse_PAF_VarTSV(i_Var_TSV):

    content = "\n".join( line for line in pathlib.Path(i_Var_TSV).read_text().split("\n") if line.startswith("V"))
    IO_Content = io.StringIO(content)

    
    i_Var_DF = pd.read_csv(IO_Content, sep = "\t", 
                            header=None)

    # PafVar_Cols_1to12 = ["VariantTag",
    #                      "Query_Name", "Query_Start", "Query_End",
    #                      "Cov", "MapQ", "Ref", "Alt",
    #                      "Target_Name", "Target_Start", "Target_End",
    #                      "Strand"]
    PafVar_Cols_1to12 = ["VariantTag",
                         "Target_Name", "Target_Start", "Target_End",
                         "Cov", "MapQ", "Ref", "Alt",
                         "Query_Name", "Query_Start", "Query_End",
                         "Strand"]


    i_Var_DF.columns = PafVar_Cols_1to12

    # Remove the R rows of the "paftools.js" var.tsv output
    i_Var_DF = i_Var_DF.query("VariantTag == 'V'").copy()

    i_Var_DF = i_Var_DF.sort_values(["Target_Start", "Target_End", "Query_Start", "Query_End", "Strand"], ascending=True)

    # Capitilize nucleotide letters
    i_Var_DF['Ref'] = i_Var_DF['Ref'].str.upper()
    i_Var_DF['Alt'] = i_Var_DF['Alt'].str.upper()

    # Set type of coordinate columns to 'int64'
    for col in ["Target_Start", "Target_End", "Query_Start", "Query_End"]:
        i_Var_DF[col] = i_Var_DF[col].astype("int64")

    i_Var_DF['SNP']  = (i_Var_DF["Alt"].str.len() == 1) & (i_Var_DF["Ref"].str.len() == 1) & (i_Var_DF["Ref"] != '-') & (i_Var_DF["Alt"] != '-')
    i_Var_DF["Type"] = i_Var_DF.apply(lambda r: infer_NT_var_type(r), axis=1)

    return i_Var_DF
    

#### Define function for adding Gene Codon info to Variant info DF ####


def addCodonInfo_To_Var(i_Var_DF, i_H37Rv_Genes_DF):
    listOfRows = []

    for i, row in tqdm(i_Var_DF.iterrows()):
        i_pos_0 = row["Start_0"]
        i_pos_1 = i_pos_0 + 1
        
        # Define SNP range
        SNP_Range = f"NC_000962.3:{i_pos_0}-{i_pos_1}"

        # 
        i_Gene_Info = bf.select(i_H37Rv_Genes_DF, SNP_Range, cols=["Chrom", "Start", "End"])

        row["N_OverlapGenes"] = i_Gene_Info.shape[0]

        if i_Gene_Info.shape[0] > 0:
            Gene_Genomic_Start_0based = i_Gene_Info["Start"].values[0]
            Gene_Genomic_End = i_Gene_Info["End"].values[0]
            Gene_Genomic_Strand = i_Gene_Info["Strand"].values[0]
            Gene_Symbol = i_Gene_Info["Symbol"].values[0]

            if Gene_Genomic_Strand == "+":
                WithinGene_Pos_0 = i_pos_0 - Gene_Genomic_Start_0based
            elif Gene_Genomic_Strand == "-":
                WithinGene_Pos_0 = Gene_Genomic_End - i_pos_1

            CodonNum = (WithinGene_Pos_0 // 3) + 1
            CodonPos = WithinGene_Pos_0 % 3 + 1

            row["Symbol"] = Gene_Symbol
            row["Strand"] = Gene_Genomic_Strand
            row["Gene_Pos_0"] = WithinGene_Pos_0
            row["Codon"] = CodonNum
            row["Codon_Pos"] = CodonPos

        listOfRows.append(row)

    i_Var_WiCodon_DF = pd.DataFrame(listOfRows)

    i_Var_WiCodon_DF["Symbol"] = i_Var_WiCodon_DF["Symbol"].fillna("None")

    return i_Var_WiCodon_DF

from Bio.Seq import Seq

def reverse_complement_base(base):
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return complement[base.upper()]


def infer_SNP_CDS_Consequences(i_SNPs_WiCodonInfo_DF,
                               dictOf_Rv_Gene_Seq,
                               i_Symbol_To_RvID_Dict):

    
    # Group by EventID and Symbol
    grouped = i_SNPs_WiCodonInfo_DF.query("Symbol != 'None' & SNP == True").groupby(['SampleID', 'Symbol'])
    
    # Dictionary to store the results
    mutation_consequences = {}

    for event_id, group in tqdm(grouped):
        i_GeneSymbol = group["Symbol"].values[0]
        i_RvID = i_Symbol_To_RvID_Dict[i_GeneSymbol]

        try:
            ref_dna_seq = Seq(dictOf_Rv_Gene_Seq[i_RvID])
            ref_protein_seq = ref_dna_seq.translate()
    
            # Convert Seq object to a list to make it mutable
            mutated_dna_seq = list(ref_dna_seq)

            # Apply each mutation
            for index, row in group.iterrows():
                
                position = int(row['Gene_Pos_0'])  # Use 0-based index in gene
                
                alt_Allele_PlusStrand = row['Alt']
                mutated_dna_seq[position] = alt_Allele_PlusStrand 
                
                if row["Strand"] == '+':
                    mutated_dna_seq[position] = alt_Allele_PlusStrand  # Apply mutation to gene sequence (+ strand)

                elif row["Strand"] == '-':
                    mutated_dna_seq[position] = reverse_complement_base(alt_Allele_PlusStrand)  # Apply mutation to gene sequence (- strand)
                else:
                    print(f"Error inserting mutations for  {event_id} - {i_GeneSymbol} - {i_RvID} - {group.shape[0]}. NO STRAND INFO for variant that is causing ISSUE!")


            # Convert the mutated list back to a Seq object
            mutated_dna_seq = Seq("".join(mutated_dna_seq))
            
            # Translate the mutated DNA sequence
            mutated_protein_seq = mutated_dna_seq.translate(table="Bacterial")

            # Compare with reference protein sequence
            differences = [(i_GeneSymbol, i+1, ref_aa, mut_aa) for i, (ref_aa, mut_aa) in enumerate(zip(ref_protein_seq, mutated_protein_seq)) if ref_aa != mut_aa]
            
            # Store the consequences
            mutation_consequences[event_id] = differences
                        
        except:
            print(f"Error translating and inferring AA changes for {event_id} - {i_GeneSymbol} - {i_RvID} - {group.shape[0]} mutations to be inserted")

    # Initialize a list to store the data
    MutCons_data = []

    # Loop through the dictionary and unpack each mutation
    for event_tuple, mutations in mutation_consequences.items():
        for Symbol, position, ref_aa, mut_aa in mutations:
            MutCons_data.append([event_tuple[0], event_tuple[1], position, ref_aa, mut_aa])

    # Create a DataFrame
    MutCon_DF = pd.DataFrame(MutCons_data, columns=['SampleID', 'Symbol', 'Codon', 'Ref_AA', 'Mut_AA'])

    return MutCon_DF









#### Functions for comparing lists of genes 
def check_overlap_with_gene_group(overlap_str, group_gene_ids):
    """
    Checks whether any gene in a comma-separated string overlaps with a given group.

    Parameters:
    -----------
    overlap_str : str
        Comma-separated gene IDs, e.g. "Rv0010c,Rv0011,Rv0020".
    group_gene_ids : set
        A set of gene IDs in the group, e.g. {"Rv0010c", "Rv1234"}.

    Returns:
    --------
    bool
        True if any gene in overlap_str is in group_gene_ids.
    """
    if pd.isna(overlap_str) or not overlap_str:
        return False
    
    group_gene_ids_Set = set(group_gene_ids)

    overlap_genes = set(overlap_str.split(","))

    return not overlap_genes.isdisjoint(group_gene_ids_Set)






def gini_coefficient_FromPDSeries(x: pd.Series) -> float:
    """
    Compute the Gini coefficient of a pandas Series.

    Assumes non-negative values.
    Returns np.nan if the Series is empty or sums to zero.
    """
    x = x.dropna().to_numpy()

    if len(x) == 0:
        return np.nan

    if np.any(x < 0):
        raise ValueError("Gini coefficient is not defined for negative values.")

    total = x.sum()
    if total == 0:
        return 0.0

    x_sorted = np.sort(x)
    n = len(x)

    index = np.arange(1, n + 1)
    gini = (2 * np.sum(index * x_sorted)) / (n * total) - (n + 1) / n

    return gini


