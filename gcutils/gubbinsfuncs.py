# pgqc/module1.py

import pandas as pd
import numpy as np

import bioframe as bf

from io import StringIO
from Bio import SeqIO
from Bio.Seq import Seq
from tqdm import tqdm

from .general import check_overlap_with_gene_group


def get_RecombEvents_From_Gubbins_GFF(input_GFF_PATH):
    listOf_Tuples = []

    counter = 1

    with open(input_GFF_PATH,'r') as fin:        
        for line in fin:

            line_Columns = line.rstrip("\n").split("\t")

            #print(counter, line_Columns)

            counter += 1
            #if counter > 10: break


            if line.startswith("#"): continue

            if len(line_Columns) < 2: continue


            #elif (True): continue
            elif (line_Columns[2] == "CDS" ): 

                NEW_Line_Columns = line_Columns[:7]

                i_Attribute_Column = line_Columns[8]

                ###For each line of your GTF, create a dictionnary with this array key : info, value : value of this info
                dict_Attr = {}
                for i in i_Attribute_Column.split(";"):

                    #print(i)

                    ###Just looking for line with " " character (as key = value)
                    if "=" in i:

                        key = i.strip().split("=")[0].strip('"')

                        value = i.strip().split("=")[1].strip('"')

                        ###Put them in a dictionnary
                        dict_Attr[key]=value


                i_nodes_Text = dict_Attr["node"]

                From_Node = i_nodes_Text.split("->")[0]
                To_Node = i_nodes_Text.split("->")[1]


                i_neg_log_likelihood = dict_Attr["neg_log_likelihood"]

                i_taxa_List = dict_Attr["taxa"].split(" ")
                i_taxa_List = ' '.join(i_taxa_List).split()

                i_snp_count = dict_Attr["snp_count"]

                NEW_Line_Columns.append(From_Node)
                NEW_Line_Columns.append(To_Node)
                NEW_Line_Columns.append(i_neg_log_likelihood)
                NEW_Line_Columns.append(i_snp_count)
                NEW_Line_Columns.append(i_taxa_List)


                NEW_Line_Columns_Tuple = tuple(NEW_Line_Columns)

                listOf_Tuples.append( NEW_Line_Columns_Tuple  )

    GFF_DF = pd.DataFrame(listOf_Tuples)

    GFF_DF.columns = ["seqname",
                      "source",
                      "feature",
                      "start_1based",
                      "end_1based",
                      "score",
                      "strand",
                      "Parent_Node",
                      "Child_Node",
                      "neg_log_likelihood",
                      "snp_count",
                      "taxa_List"] 

    return GFF_DF





def get_RecombEvents_From_Gubbins_GFF_V2(input_GFF_PATH):
    listOf_Tuples = []

    counter = 1

    with open(input_GFF_PATH,'r') as fin:        
        for line in fin:

            line_Columns = line.rstrip("\n").split("\t")

            #print(counter, line_Columns)

            counter += 1
            #if counter > 10: break


            if line.startswith("#"): continue

            if len(line_Columns) < 2: continue


            #elif (True): continue
            elif (line_Columns[2] == "CDS" ): 

                NEW_Line_Columns = line_Columns[:7]

                i_Attribute_Column = line_Columns[8]

                ###For each line of your GTF, create a dictionnary with this array key : info, value : value of this info
                dict_Attr = {}
                for i in i_Attribute_Column.split(";"):

                    #print(i)

                    ###Just looking for line with " " character (as key = value)
                    if "=" in i:

                        key = i.strip().split("=")[0].strip('"')

                        value = i.strip().split("=")[1].strip('"')

                        ###Put them in a dictionnary
                        dict_Attr[key]=value


                i_nodes_Text = dict_Attr["node"]

                From_Node = i_nodes_Text.split("->")[0]
                To_Node = i_nodes_Text.split("->")[1]


                i_neg_log_likelihood = dict_Attr["neg_log_likelihood"]

                i_taxa_List = dict_Attr["taxa"].split(" ")
                i_taxa_List = ' '.join(i_taxa_List).split()

                i_snp_count = dict_Attr["snp_count"]

                NEW_Line_Columns.append(From_Node)
                NEW_Line_Columns.append(To_Node)
                NEW_Line_Columns.append(i_neg_log_likelihood)
                NEW_Line_Columns.append(i_snp_count)
                NEW_Line_Columns.append(i_taxa_List)


                NEW_Line_Columns_Tuple = tuple(NEW_Line_Columns)

                listOf_Tuples.append( NEW_Line_Columns_Tuple  )

    i_GRE_DF = pd.DataFrame(listOf_Tuples)

    i_GRE_DF.columns = ["seqname",
                      "source",
                      "feature",
                      "start_1based",
                      "end_1based",
                      "score",
                      "strand",
                      "Parent_Node",
                      "Child_Node",
                      "neg_log_likelihood",
                      "snp_count",
                      "taxa_List"] 

    i_GRE_DF["start_1based"] = i_GRE_DF["start_1based"].astype(int)
    i_GRE_DF["end_1based"] = i_GRE_DF["end_1based"].astype(int)

    i_GRE_DF["start_0based"] = i_GRE_DF["start_1based"] - 1

    i_GRE_DF["CenterOfRegion"] = ((i_GRE_DF["start_1based"] + i_GRE_DF["end_1based"] - 1) / 2)
    i_GRE_DF["EventLen"] = i_GRE_DF["end_1based"] - i_GRE_DF["start_1based"] + 1

    i_GRE_DF["snp_count"] = i_GRE_DF["snp_count"].astype(int)


    # Replace very large Neg Log Likelihoods
    #### Replace "-nan" with a very large #, this is for cases where the neg log likelihood was too large to calculate
    i_GRE_DF["neg_log_likelihood"] = i_GRE_DF["neg_log_likelihood"].replace("-nan", 999999)  
    i_GRE_DF["neg_log_likelihood"] = i_GRE_DF["neg_log_likelihood"].astype(float)  

    i_GRE_DF = i_GRE_DF.drop(columns=["source", "feature", "score"])

    return i_GRE_DF





def subset_EventsDF_ForGenes(i_GRE_DF, gene_set):
    # Normalize input: if a single string is passed, wrap it in a list
    if isinstance(gene_set, str):
        gene_set = [gene_set]

    # Build boolean mask
    bool_array_gene_set = i_GRE_DF["Overlap_Genes"].str.contains('|'.join(gene_set))

    # Subset
    return i_GRE_DF[bool_array_gene_set]







def label_Gubbins_Events_DF_ByOvrLap_H37RvGenes(i_Gubbins_Events_DF,
                                                i_H37Rv_Genes_DF, 
                                                Start_Column="start_0based",
                                                End_Column="end_1based",
                                                GeneNameColumn = "Symbol"):
    
    overlap_genes = []
    overlap_rv_ids = []

    for _, row in i_Gubbins_Events_DF.iterrows():
        event_start = int(row[Start_Column])
        event_end = int(row[End_Column])
        event_range = f"NC_000962.3:{event_start}-{event_end}"

        overlapping_genes_df = bf.select(i_H37Rv_Genes_DF, event_range, cols=("Chrom", "Start", "End"))

        if not overlapping_genes_df.empty:
            genes = ",".join(overlapping_genes_df["Symbol"].astype(str))
            rv_ids = ",".join(overlapping_genes_df["H37rv_GeneID"].astype(str))
        else:
            genes = "_"
            rv_ids = "_"

        overlap_genes.append(genes)
        overlap_rv_ids.append(rv_ids)

    i_Gubbins_Events_WiGenes_DF = i_Gubbins_Events_DF.copy()
    i_Gubbins_Events_WiGenes_DF["Overlap_Genes"] = overlap_genes
    i_Gubbins_Events_WiGenes_DF["Overlap_Gene_RvIDs"] = overlap_rv_ids

    i_Gubbins_Events_WiGenes_DF["Overlap_Genes"] = i_Gubbins_Events_WiGenes_DF["Overlap_Genes"].fillna("_")
    i_Gubbins_Events_WiGenes_DF["Overlap_Gene_RvIDs"] = i_Gubbins_Events_WiGenes_DF["Overlap_Gene_RvIDs"].fillna("_")

    return i_Gubbins_Events_WiGenes_DF



def label_Gubbins_Events_DF_ByOvrLap_HHRs(i_GRE_DF, i_HHR_DF,
                                          HHR_ID_Column = "HmRegionID",
                                          GCE_CoordCols = ("seqname", "start_0based", "end_1based"),
                                          HHR_CoordCols = ("Chr", "Start", "End")):
    
    listOf_Overlap_HHRs = []

    for i, row in i_GRE_DF.iterrows():
        
        # a) Target overlapping genes
        i_Chrom = row[ GCE_CoordCols[0] ]
        i_Start = int(row[ GCE_CoordCols[1] ])
        i_End =  int(row[ GCE_CoordCols[2] ])

        Target_Range = f"{i_Chrom}:{i_Start}-{i_End}"

        sub_DF_Overlap_HHRs = bf.select(i_HHR_DF, Target_Range, cols = HHR_CoordCols)

        listOf_Overlap_HHRs.append( ",".join(list(sub_DF_Overlap_HHRs[ HHR_ID_Column ].values)) )


    i_GRE_Anno_DF = i_GRE_DF.copy()
    i_GRE_Anno_DF["Overlap_HHRs"] = listOf_Overlap_HHRs

    # For all cases where there are NO OVERLAPPING GENES, simply put "None"
    i_GRE_Anno_DF["Overlap_HHRs"] = i_GRE_Anno_DF["Overlap_HHRs"].fillna("_") #.replace("", "None")

    return i_GRE_Anno_DF





def add_HighHomologyAndRepeat_OvrlapInfo_ToGubbinsEvents(i_Events_DF,
                                                         i_DictOf_HighHomologyAndRepeat_InfoDFs ):


    col_to_drop = set(i_Events_DF.columns) & {
        "OvrlapWi_HmMap_ParalogReg",
        "OvrlapWi_HmMap_LocalRepeat",
        "OvrlapWi_LowComplexityRegion",
        "OvrlapWi_LowPmapRegion",
        "N_HmMapAln_PR_Ovrlap",
        "N_HmMapAln_LR_Ovrlap",
        "N_LCR_Ovrlap",
        "N_LowPmap_Ovrlap",
    }

    i_Events_DF = i_Events_DF.drop(columns=list(col_to_drop))
    
    RE_CoordCols = ("seqname", "start_0based", "end_1based")
    HmAlnPAF_CoordCols = ("Query_Name", "Query_Start", "Query_End")
    HHR_CoordCols = ("Chr", "Start", "End")

    
    i_HmMap_LocalRepeats_PairwiseAln_DF = i_DictOf_HighHomologyAndRepeat_InfoDFs["HmMap_LocalRepeat_Aln_DF"]
    
    i_HmMap_Paralogs_PairwiseAln_DF     = i_DictOf_HighHomologyAndRepeat_InfoDFs["HmMap_Paralogous_Aln_DF"]

    i_LowComp_Regions_DF = i_DictOf_HighHomologyAndRepeat_InfoDFs["LowComplexity_Regions_DF"]
    i_LowPMap_Regions_DF = i_DictOf_HighHomologyAndRepeat_InfoDFs["LowPmap_Regions_DF"]

    i_HmMap_Paralogous_Regions_DF   = i_DictOf_HighHomologyAndRepeat_InfoDFs["HmMap_Paralogous_Regions"]
    i_HmMap_LocalRepeats_Regions_DF = i_DictOf_HighHomologyAndRepeat_InfoDFs["HmMap_LocalRepeat_Regions"]

    i_HmMap_LRandPR_Regions_DF = pd.concat([i_HmMap_Paralogous_Regions_DF,
                                            i_HmMap_LocalRepeats_Regions_DF])


    #### Step 1: Add a column which uses the unique IDs for the Local Repeat and ParalogRegions (From Homology-Map)
    
    i_Events_DF = label_Gubbins_Events_DF_ByOvrLap_HHRs(i_Events_DF,
                                                        i_HmMap_LRandPR_Regions_DF) 



    #### Step 2:
    i_Events_DF = bf.count_overlaps(i_Events_DF,
                                    i_HmMap_Paralogs_PairwiseAln_DF, 
                                    cols1 = RE_CoordCols,
                                    cols2 = HmAlnPAF_CoordCols ).rename(columns={'count': 'N_HmMapAln_PR_Ovrlap'})
    
    
    i_Events_DF["OvrlapWi_HmMap_ParalogReg"] = np.where(i_Events_DF['N_HmMapAln_PR_Ovrlap'] > 0, 1, 0)
    
    
    #### Step 3:
    i_Events_DF = bf.count_overlaps(i_Events_DF,
                                    i_HmMap_LocalRepeats_PairwiseAln_DF, 
                                    cols1 = RE_CoordCols,
                                    cols2 = HmAlnPAF_CoordCols ).rename(columns={'count': 'N_HmMapAln_LR_Ovrlap'})
    
    
    i_Events_DF["OvrlapWi_HmMap_LocalRepeat"] = np.where(i_Events_DF['N_HmMapAln_LR_Ovrlap'] > 0, 1, 0)
    
    
    ## Annotate each event by overlaps w/ `LCR` or `LowPmap` Regions of H37Rv
    
    #### Step 4: Annotate each event by overlap w/ `LCR` (low complexity regions identified with longdust)
    
    i_Events_DF = bf.count_overlaps(i_Events_DF,
                                    i_LowComp_Regions_DF, 
                                    cols1 = RE_CoordCols,
                                    cols2 = ('chrom', 'start', 'end') ).rename(columns={'count': 'N_LCR_Ovrlap'})
    
    i_Events_DF["OvrlapWi_LowComplexityRegion"] = np.where(i_Events_DF['N_LCR_Ovrlap'] > 0, 1, 0)
    
    #### Step 5: Annotate each event by overlap with Low Pileup Mappability scores (Based on genmap and pupmapper)
    
    i_Events_DF = bf.count_overlaps(i_Events_DF,
                                                 i_LowPMap_Regions_DF, 
                                                 cols1 = RE_CoordCols,
                                                 cols2 = ('chrom', 'start', 'end')).rename(columns={'count': 'N_LowPmap_Ovrlap'})
    
    i_Events_DF["OvrlapWi_LowPmapRegion"] = np.where(i_Events_DF['N_LowPmap_Ovrlap'] > 0, 1, 0)
    

    return i_Events_DF

    


def Parse_Process_Gubbins_Events_Standard(in_Gubbins_Events_GFF_PATH,
                                          in_TreeNode_LineageLabels_Dict,
                                          in_TreeNode_NumDescendants_Dict,
                                          in_GenomeAnno_Genes_DF,
                                          in_Esx_PEPPE_REP13E12_GeneLists_Dict,
                                          i_DictOf_HighHomologyAndRepeat_InfoDFs,
                                          verbose = False):
    
    # Step 1: Parse the events predicted by Gubbins from GFF file
    i_Events_DF = get_RecombEvents_From_Gubbins_GFF_V2(in_Gubbins_Events_GFF_PATH)


    # Step 2: Label Gubbins Events by lineage of tree branch it occured on
    i_Events_DF["Lineage"] = i_Events_DF["Child_Node"].map(in_TreeNode_LineageLabels_Dict).fillna("None")


    i_Events_DF["Num_Tips_Downstream"] = i_Events_DF["Child_Node"].map(in_TreeNode_NumDescendants_Dict).fillna("None")

    
    # Step 3: Label Gubbins Events by overlapping genes
    i_Events_DF = label_Gubbins_Events_DF_ByOvrLap_H37RvGenes(i_Events_DF, 
                                                              in_GenomeAnno_Genes_DF)


    # Step 4: Sort all events by coordinates, phylogenetic branch nodeIDs, + overlap genes
    
    ## - Let's sort by the following columns: ["start_1based", "end_1based", "Parent_Node", "Child_Node", "Lineage", "Overlap_Genes"]
    Gub_ColumnsToSortBy = ["start_1based", "end_1based", "Parent_Node", "Child_Node", "Lineage", "Overlap_Genes"]
    
    i_Events_DF = i_Events_DF.sort_values(Gub_ColumnsToSortBy, kind="mergesort").reset_index(drop=True)
    
    i_Events_DF["EventNum"] = i_Events_DF.index + 1
    i_Events_DF["EventID"] = "Event_" + i_Events_DF["EventNum"].astype(str).str.rjust(3, '0')
    
    i_Events_DF = i_Events_DF.drop(["EventNum"], axis = 1)

    
    # Step 5: Let's label each putative recombination event by whether it overlaps with one of the 3 gene-groups (Esx, PE/PPE, REP13E12 repeats, or None) 

    i_ListOf_Esx_RvIDs = in_Esx_PEPPE_REP13E12_GeneLists_Dict["Esx"]
    i_ListOf_PEPPE_RvIDs = in_Esx_PEPPE_REP13E12_GeneLists_Dict["PEPPE"]
    i_ListOf_REP13E12_RvIDs = in_Esx_PEPPE_REP13E12_GeneLists_Dict["REP13E12"]

    G = i_Events_DF.copy()
    G["Contains_Esx"] = G["Overlap_Gene_RvIDs"].apply(lambda x: check_overlap_with_gene_group(x, i_ListOf_Esx_RvIDs))
    G["Contains_PEPPE"] = G["Overlap_Gene_RvIDs"].apply(lambda x: check_overlap_with_gene_group(x, i_ListOf_PEPPE_RvIDs))
    G["Contains_REP13E12"] = G["Overlap_Gene_RvIDs"].apply(lambda x: check_overlap_with_gene_group(x, i_ListOf_REP13E12_RvIDs))
    G["NoOverlap_Wi_PEPPE_Esx_13E12_Genes"] = ~(G["Contains_Esx"] | G["Contains_PEPPE"] | G["Contains_REP13E12"] )
    
    i_Events_DF = G.copy()


    # Step 6: 
    i_Events_DF = add_HighHomologyAndRepeat_OvrlapInfo_ToGubbinsEvents(i_Events_DF,
                                                                       i_DictOf_HighHomologyAndRepeat_InfoDFs )

    
    return i_Events_DF





















########################################################################################################################










def annotate_HHR_with_GCE_overlaps(
    i_HHR_DF,
    i_GCE_DF,
    HHR_CoordCols=("Chr", "Start", "End"),
    GCE_CoordCols=("seqname", "start_0based", "end_1based"),
    GCE_ID_Col="EventID"
):
    """
    Annotate each row in i_HHR_DF with comma-separated GC EventIDs from i_GCE_DF that overlap.

    Parameters:
    - i_HHR_DF: DataFrame of high-homology regions.
    - i_GCE_DF: DataFrame of gene conversion events.
    - HHR_CoordCols: Tuple of (chrom, start, end) column names in i_HHR_DF.
    - GCE_CoordCols: Tuple of (chrom, start, end) column names in i_GCE_DF.
    - GCE_ID_Col: Column name in i_GCE_DF that contains unique GC Event IDs.

    Returns:
    - Annotated DataFrame with an "Overlap_GC_EventIDs" column.
    """

    listOf_Overlap_GCEs = []

    for _, row in i_HHR_DF.iterrows():
        i_Chrom = row[HHR_CoordCols[0]]
        i_Start = int(row[HHR_CoordCols[1]])
        i_End = int(row[HHR_CoordCols[2]])
        Target_Range = f"{i_Chrom}:{i_Start}-{i_End}"

        overlapping = bf.select(i_GCE_DF, Target_Range, cols=GCE_CoordCols)
        overlap_ids = ",".join(overlapping[GCE_ID_Col].astype(str)) if not overlapping.empty else "_"

        listOf_Overlap_GCEs.append(overlap_ids)

    annotated_df = i_HHR_DF.copy()
    annotated_df["Overlap_GC_EventIDs"] = listOf_Overlap_GCEs
    
    annotated_df["Overlap_GC_EventIDs"] = annotated_df["Overlap_GC_EventIDs"].fillna("_")


    return annotated_df






##### Define functions for parsing and analyzing the SNP mutation events from the Gubbins Ancestral State Reconstruction ######


def parse_BaseReconstruction_EMBL_Gubbins(input_BaseReconstr_EMBL_PATH):
    """ Parse the {Name}.branch_base_reconstruction.embl from Gubbins"""
    
    # read as string, modify header and footer
    data = open(input_BaseReconstr_EMBL_PATH, 'r').read()
    head = 'ID   X56734; standard; linear; DNA; STD; PLN; 999 BP.\nFH   Key             Location/Qualifiers\n'
    foot = '\nXX\nSQ   Sequence 999 BP;\n//'
    data_w_headers = f'{head}{data}{foot}'

    # read into seqio object
    file_io = StringIO(data_w_headers)
    records = list(SeqIO.parse(file_io, 'embl'))

    results = []
    for record in records:
        for feature in record.features:
            if feature.type == 'variation':
                #print(feature)
                i_Pos_0 = feature.location.start
                i_Pos_1 = i_Pos_0 + 1

                i_Parent_Call = feature.qualifiers['parent_base'][0]
                i_Child_Call = feature.qualifiers['replace'][0]

                i_Taxa = feature.qualifiers['taxa'][0]


                i_AffectedNodes = feature.qualifiers['node'][0]
                i_Parent_Node = i_AffectedNodes.replace("->", " ").split(" ")[0]

                i_Child_Node = i_AffectedNodes.replace("->", " ").split(" ")[1]


                results.append([i_Pos_1, i_Parent_Node, i_Child_Node, i_Parent_Call, i_Child_Call, i_Taxa])

    i_BR_DF = pd.DataFrame(results)

    i_BR_DF.columns = ["Pos_1based", "Parent_Node", "Child_Node", "Parent_Call", "Child_Call", "taxa_List"]
    
    return i_BR_DF



def annotate_Gubbins_SNP_Events_By_RecombEventID(i_Gubbins_SNPs_DF, i_GRE_DF):
    """
    Annotate SNP events with overlapping recombination event IDs from Gubbins.

    Parameters:
    -----------
    i_Gubbins_SNPs_DF : pd.DataFrame
        DataFrame of all SNP mutation events inferred by Gubbins.
        Must include columns: 'Parent_Node', 'Child_Node', 'Pos_1based'.

    i_GRE_DF : pd.DataFrame
        DataFrame of recombination events inferred by Gubbins.
        Must include columns: 'EventID', 'start_1based', 'end_1based',
                              'Parent_Node', 'Child_Node'.

    Returns:
    --------
    pd.DataFrame
        Updated SNP DataFrame with an added 'EventID' column indicating
        the event that each SNP overlaps with (or "None" if no overlap).
    """
    listOf_Labeled_DFs = []

    for _, event_row in i_GRE_DF.iterrows():
        # Extract recombination event details
        event_id = event_row["EventID"]
        start = event_row["start_1based"]
        end = event_row["end_1based"]
        parent = event_row["Parent_Node"]
        child = event_row["Child_Node"]

        # Find SNPs that match both the node pair and coordinate range
        overlapping_snps = i_Gubbins_SNPs_DF.query(
            f"(Parent_Node == '{parent}') & (Child_Node == '{child}') & "
            f"(Pos_1based >= {start}) & (Pos_1based <= {end})"
        ).copy()

        overlapping_snps["EventID"] = event_id
        listOf_Labeled_DFs.append(overlapping_snps)

    # Concatenate all labeled SNPs
    labeled_snps_df = pd.concat(listOf_Labeled_DFs, ignore_index=True) if listOf_Labeled_DFs else pd.DataFrame()

    # Merge back with the full SNPs DF to include unmatched SNPs
    SNPs_full_annotated_df = pd.merge(
        i_Gubbins_SNPs_DF,
        labeled_snps_df[["Parent_Node", "Child_Node", "Pos_1based", "EventID"]],
        on=["Parent_Node", "Child_Node", "Pos_1based"],
        how="left"
    )

    SNPs_full_annotated_df["Pos_0based"] = SNPs_full_annotated_df["Pos_1based"] - 1
    SNPs_full_annotated_df["Chrom"] = "NC_000962.3"


    # Fill non-overlapping SNPs with "None"
    SNPs_full_annotated_df["EventID"] = SNPs_full_annotated_df["EventID"].fillna("None")

    # Sort for consistent output
    sorting_cols = ["Parent_Node", "Child_Node", "EventID", "Pos_1based"]
    SNPs_full_annotated_df = SNPs_full_annotated_df.sort_values(sorting_cols, kind="mergesort").reset_index(drop=True)

    return SNPs_full_annotated_df

########################################################################################################################




def annotate_Gubbins_SNP_Events_By_EventID_And_ParalogRegion(i_Gubbins_SNPs_DF,
                                                             i_GRE_DF,
                                                             i_HmMap_Paralogous_Aln_DF):
    """
    Annotate SNP events with overlapping recombination event IDs from Gubbins.

    Parameters:
    -----------
    i_Gubbins_SNPs_DF : pd.DataFrame
        DataFrame of all SNP mutation events inferred by Gubbins.
        Must include columns: 'Parent_Node', 'Child_Node', 'Pos_1based'.

    i_GRE_DF : pd.DataFrame
        DataFrame of recombination events inferred by Gubbins.
        Must include columns: 'EventID', 'start_1based', 'end_1based',
                              'Parent_Node', 'Child_Node'.

    Returns:
    --------
    pd.DataFrame
        Updated SNP DataFrame with an added 'EventID' column indicating
        the event that each SNP overlaps with (or "None" if no overlap).
    """
    listOf_Labeled_DFs = []

    for _, event_row in i_GRE_DF.iterrows():
        # Extract recombination event details
        event_id = event_row["EventID"]
        start = event_row["start_1based"]
        end = event_row["end_1based"]
        parent = event_row["Parent_Node"]
        child = event_row["Child_Node"]

        # Find SNPs that match both the node pair and coordinate range
        overlapping_snps = i_Gubbins_SNPs_DF.query(
            f"(Parent_Node == '{parent}') & (Child_Node == '{child}') & "
            f"(Pos_1based >= {start}) & (Pos_1based <= {end})"
        ).copy()

        overlapping_snps["EventID"] = event_id
        listOf_Labeled_DFs.append(overlapping_snps)

    # Concatenate all labeled SNPs
    labeled_snps_df = pd.concat(listOf_Labeled_DFs, ignore_index=True) if listOf_Labeled_DFs else pd.DataFrame()

    # Merge back with the full SNPs DF to include unmatched SNPs
    SNPs_full_annotated_df = pd.merge(
        i_Gubbins_SNPs_DF,
        labeled_snps_df[["Parent_Node", "Child_Node", "Pos_1based", "EventID"]],
        on=["Parent_Node", "Child_Node", "Pos_1based"],
        how="left"
    )

    SNPs_full_annotated_df["Pos_0based"] = SNPs_full_annotated_df["Pos_1based"] - 1
    SNPs_full_annotated_df["Chrom"] = "NC_000962.3"


    # Fill non-overlapping SNPs with "None"
    SNPs_full_annotated_df["EventID"] = SNPs_full_annotated_df["EventID"].fillna("None")

    # Sort for consistent output
    sorting_cols = ["Parent_Node", "Child_Node", "EventID", "Pos_1based"]
    SNPs_full_annotated_df = SNPs_full_annotated_df.sort_values(sorting_cols, kind="mergesort").reset_index(drop=True)


    Gub_SNP_CoordCols = ("Chrom", "Pos_0based", "Pos_1based")
    Hm_CoordCols = ("Query_Name", "Query_Start", "Query_End")

    SNPs_full_annotated_df = bf.count_overlaps(SNPs_full_annotated_df,
                                               i_HmMap_Paralogous_Aln_DF,
                                               cols1 = Gub_SNP_CoordCols,
                                                cols2 = Hm_CoordCols )

    SNPs_full_annotated_df.rename(columns={'count': 'Num_ParalogAln_Ovrlap'}, inplace=True)

    SNPs_full_annotated_df["PR_HmOvrlap"]  = np.where(SNPs_full_annotated_df['Num_ParalogAln_Ovrlap'] > 0, 1, 0).astype(int)

    SNPs_full_annotated_df["RegionType"] = np.where(SNPs_full_annotated_df['Num_ParalogAln_Ovrlap'] > 0, "PR", "Unq").astype(str)


    return SNPs_full_annotated_df






























############ Define functions for annotation of Gubbins SNP events by CDS consequences/effect ############

def reverse_complement_base(base):
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return complement[base.upper()]




#### Define function for adding Gene Codon info to Gubbins Var DF
def addCDSInfo_Gubbins_ASR_SNP(i_Gubbins_SNPs_DF,
                                 i_H37Rv_Genes_DF):
    listOfRows = []

    for i, row in tqdm(i_Gubbins_SNPs_DF.iterrows()):
        i_pos_1 = row["Pos_1based"]
        i_pos_0 = i_pos_1 - 1

        # Define SNP range
        SNP_Range = f"NC_000962.3:{i_pos_0}-{i_pos_1}"

        # Step 1: Get the list of all BEST Donor matches for each event
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

    GRE_SNPs_DF = pd.DataFrame(listOfRows)

    GRE_SNPs_DF["Symbol"] = GRE_SNPs_DF["Symbol"].fillna("None")
    GRE_SNPs_DF = GRE_SNPs_DF.drop('taxa_List', axis=1, errors='ignore')

    return GRE_SNPs_DF




def calculate_mutation_consequences_Gubbins_SNPs(i_SNPs_WiCodonInfo_DF, dictOf_Rv_Gene_Seq, i_Symbol_To_RvID_Dict):
    # Group by EventID and Symbol
    grouped = i_SNPs_WiCodonInfo_DF.query("Symbol != 'None'").groupby(['Child_Node', 'Symbol'])
    print(i_SNPs_WiCodonInfo_DF.columns)
    #print(grouped.columns)

    # Dictionary to store the results
    mutation_consequences = {}

    for event_id, group in tqdm(grouped):
        #print(event_id)

        i_GeneSymbol = group["Symbol"].values[0]
        i_RvID = i_Symbol_To_RvID_Dict[i_GeneSymbol]
        try:
            ref_dna_seq = Seq(dictOf_Rv_Gene_Seq[i_RvID])
            ref_protein_seq = ref_dna_seq.translate()
    
            # Convert Seq object to a list to make it mutable
            mutated_dna_seq = list(ref_dna_seq)
            #print(i_GeneSymbol, i_RvID)
            # Apply each mutation
            for index, row in group.iterrows():

                position = int(row['Gene_Pos_0'])  # Use 0-based index in gene

                Mut_Allele_OnPlusStrand = row['Child_Call']

                if row["Strand"] == '+':
                    mutated_dna_seq[position] = Mut_Allele_OnPlusStrand  # Apply mutation to gene sequence (+ strand)

                elif row["Strand"] == '-':
                    mutated_dna_seq[position] = reverse_complement_base(Mut_Allele_OnPlusStrand)  # Apply mutation to gene sequence (- strand)
                else:
                    print(f"Error inserting mutations for  {event_id} - {i_GeneSymbol} - {i_RvID} - {group.shape[0]}. NO STRAND INFO for variant that is causing ISSUE!")

                # mutated_dna_seq[position] = row['Child_Call']  # Apply mutation
    
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
    MutCon_DF = pd.DataFrame(MutCons_data, columns=['Child_Node', 'Symbol', 'Codon', 'Ref_AA', 'Mut_AA'])

    return MutCon_DF



def Add_AA_Consequences_ToGubbinsSNPs_DF(i_Gubbins_SNPs_DF,
                                         i_H37Rv_Genes_DF,
                                         dictOf_Rv_Gene_Seq,
                                         i_Symbol_To_RvID_Dict):

    i_SNPs_WiCodonInfo_DF = addCDSInfo_Gubbins_ASR_SNP(i_Gubbins_SNPs_DF,
                                                         i_H37Rv_Genes_DF)
    
    i_SNPs_AAChanges_DF = calculate_mutation_consequences_Gubbins_SNPs(i_SNPs_WiCodonInfo_DF,
                                                                       dictOf_Rv_Gene_Seq,
                                                                       i_Symbol_To_RvID_Dict)
    i_SNPs_AA_Anno_DF = pd.merge(i_SNPs_WiCodonInfo_DF,
                                         i_SNPs_AAChanges_DF,
                                         on = ["Child_Node", "Symbol", "Codon"], how = "left")
    
    i_SNPs_AA_Anno_DF["Chrom"] = "NC_000962.3"
    i_SNPs_AA_Anno_DF["Start"] = i_SNPs_AA_Anno_DF["Pos_1based"] - 1

    i_SNPs_AA_Anno_DF["End"] = i_SNPs_AA_Anno_DF["Pos_1based"]
    
    i_SNPs_AA_Anno_DF["MissenseMut"] = (~i_SNPs_AA_Anno_DF["Mut_AA"].isna() ) & (i_SNPs_AA_Anno_DF["Ref_AA"] != i_SNPs_AA_Anno_DF["Mut_AA"])

    return i_SNPs_AA_Anno_DF


from Bio.Seq import Seq


##########################################################################################################







###### Define functions for comparison/overlap of Gubbins SNP events DFs between different analyses ######





def compare_Gubbins_ASR_SNP_DFs(df_a, df_b):
    """
    Compare SNPs between two DataFrames and return:
    1. Number of intersecting SNPs with matching parent and child calls
    2. DataFrame of intersecting SNPs
    3. Outer-merged DataFrame of both SNP sets

    Parameters:
    - df_a: DataFrame with columns ["Pos_1based", "Parent_Call", "Child_Call", "EventID"]
    - df_b: DataFrame with columns ["Pos_1based", "Parent_Call", "Child_Call", "EventID"]

    Returns:
    - N_Intersect_SNPs: int
    - intersect_df: pd.DataFrame
    - outer_merged_df: pd.DataFrame
    """
    # Subset relevant columns
    df_a_trim = df_a[["Pos_1based", "Parent_Call", "Child_Call", "EventID"]]
    df_b_trim = df_b[["Pos_1based", "Parent_Call", "Child_Call", "EventID"]]

    # Outer merge on position
    outer_merged_df = pd.merge(df_a_trim, df_b_trim, how='outer',
                               on='Pos_1based',
                               suffixes=("_A", "_B"))

    # Filter for matching SNPs
    intersect_df = outer_merged_df.query(
        "(Parent_Call_A == Parent_Call_B) & (Child_Call_A == Child_Call_B)"
    )

    N_Intersect_SNPs = intersect_df.shape[0]

    return N_Intersect_SNPs, intersect_df, outer_merged_df

##########################################################################################################






############# Gubbins Phylogeny analysis Functions ##############

def get_BranchLengths_and_DescendantCounts_FromTree(i_ete3_tree):
    """
    Extract branch lengths and number of descendant leaf nodes for each node in an ETE3 tree.

    Parameters:
    -----------
    i_ete3_tree : ete3.Tree
        The input ETE3 tree object.

    Returns:
    --------
    branch_len_dict : dict
        Dictionary mapping node name to branch length from its parent.
    
    descendant_count_dict : dict
        Dictionary mapping node name to number of descendant leaf nodes.
    """
    branch_len_dict = {}
    descendant_count_dict = {}

    for node in i_ete3_tree.iter_descendants("postorder"):
        # Count descendant leaf nodes
        descendant_leaves = [desc.name for desc in node.get_descendants() if desc.is_leaf()]
        descendant_count_dict[node.name] = len(descendant_leaves)

        # Record branch length from parent
        branch_len_dict[node.name] = node.dist

    return branch_len_dict, descendant_count_dict







def annotate_tree_with_lineages(
    i_ete3_tree,
    i_ID_To_PrimLineage_Dict,
):
    """
    Annotate leaves of an ETE3 tree with lineage info, and infer primary
    lineage for internal nodes that are monophyletic for a single lineage.

    Parameters
    ----------
    i_ete3_tree : ete3.Tree
        The ETE3 tree object whose nodes will be annotated.

    i_ID_To_PrimLineage_Dict : dict
        Maps leaf names -> primary lineage.

    i_ID_To_SubLineage_Dict : dict
        Maps leaf names -> sublineage.

    Returns
    -------
    i_ete3_tree : ete3.Tree
        The same tree object, modified in place.

    node_To_PrimaryLin_Dict : dict
        Dictionary mapping node name -> inferred primary lineage.
        Includes both internal nodes (when monophyletic) and leaves
        (from i_ID_To_PrimLineage_Dict).
    """

    # 1) Annotate leaves with provided lineage info
    leaf_count = 0
    for n in i_ete3_tree.get_leaves():
        n.add_feature(
            "Primary_lineage",
            i_ID_To_PrimLineage_Dict.get(n.name, "Unknown Lineage")
        )

        leaf_count += 1

    i_ete3_tree.sort_descendants(attr='Primary_lineage')

    print("Number of leaves processed:", leaf_count)

    # 4) Sort tree by primary lineage (uses the feature we just added)
    i_ete3_tree.sort_descendants(attr="Primary_lineage")

    return i_ete3_tree



def infer_InternalNode_Lineages(i_ete3_tree,
                                i_ID_To_PrimLineage_Dict):

    node_To_PrimaryLin_Dict = {}

    for node in i_ete3_tree.iter_descendants("postorder"):
        
        listOf_ChildLineages = []
        
        for child_node in node.get_descendants():
            if child_node.is_leaf():
                listOf_ChildLineages.append(  (child_node.Primary_lineage) )
                    #print(node.name, listOf_ChildLineages)
            
        set_Of_ChildLineages = list(set(listOf_ChildLineages))
        
        if len(set_Of_ChildLineages) == 1:
            OnlyOneLineage = True
        else:
            OnlyOneLineage = False
        
        if OnlyOneLineage:
            node_To_PrimaryLin_Dict[node.name] = set_Of_ChildLineages[0]
    
    node_To_PrimaryLin_Dict.update(i_ID_To_PrimLineage_Dict)

    
    return i_ete3_tree, node_To_PrimaryLin_Dict














from typing import Tuple
from ete3 import Tree
import pandas as pd

def get_parsimony_stats_AllVariants_Per_GC_Event(
    tree: Tree,
    i_pGCE_DF: pd.DataFrame,
    i_GubSNPs_All_DF: pd.DataFrame,
    i_PerAsm_AllVar_DF: pd.DataFrame,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    For each GCE EventID, compute per-site and per-event parsimony statistics.
    Always returns:
        - per_site_all_df       (one row per EventID × site)
        - event_summary_all_df  (one row per EventID)
    """

    per_site_list = []
    event_stats_list = []

    # Loop over all EventIDs in the provided pGCE dataframe
    for i_EventID in i_pGCE_DF["EventID"].unique():

        # 1) Extract metadata for this event
        i_row = i_pGCE_DF.loc[i_pGCE_DF["EventID"] == i_EventID].iloc[0]
        i_ParentNode_Name = i_row["Parent_Node"]
        i_ChildNode_Name  = i_row["Child_Node"]

        # 2) Extract subtree rooted at Child_Node
        i_GCE_SubTree_AtChildNode = subset_ETE_Tree(i_ChildNode_Name, tree)

        # 3) Collect SNP positions associated with this EventID
        i_GubSNPs_TargetGCE_DF = i_GubSNPs_All_DF.loc[
            i_GubSNPs_All_DF["EventID"] == i_EventID
        ]
        i_GCE_SNP_Pos_0based = list(i_GubSNPs_TargetGCE_DF["Pos_0based"].values)

        # 4) Compute per-site parsimony for these positions
        i_per_snp_parsimony_df = parsimony_per_snp(
            tree=i_GCE_SubTree_AtChildNode,
            variants_df=i_PerAsm_AllVar_DF,
            leaf_col="SampleID",
            state_col="Alt",
            pos_col="Start_0",
            ref_col="Ref",
            alt_col="Alt",
            snp_col="SNP",
            snp_only=True,
            positions_of_interest=i_GCE_SNP_Pos_0based,
        )

        # 5) Annotate per-site DF with EventID + node metadata
        if not i_per_snp_parsimony_df.empty:
            i_per_snp_parsimony_df = i_per_snp_parsimony_df.copy()
            i_per_snp_parsimony_df["EventID"]     = i_EventID
            i_per_snp_parsimony_df["Parent_Node"] = i_ParentNode_Name
            i_per_snp_parsimony_df["Child_Node"]  = i_ChildNode_Name

            # Reorder columns: EventID, Parent_Node, Child_Node first
            first_cols = ["EventID", "Parent_Node", "Child_Node"]
            other_cols = [c for c in i_per_snp_parsimony_df.columns if c not in first_cols]
            i_per_snp_parsimony_df = i_per_snp_parsimony_df[first_cols + other_cols]

        per_site_list.append(i_per_snp_parsimony_df)

        # 6) Per-event summary
        i_GCE_ParsimonyScore_Stats_DF = summarize_parsimony_across_sites(
            per_site_df=i_per_snp_parsimony_df,
            SetID=i_EventID,   # This returns VariantSetID, but we will rename it
        )

        # Rename VariantSetID → EventID
        i_GCE_ParsimonyScore_Stats_DF = i_GCE_ParsimonyScore_Stats_DF.rename(
            columns={"VariantSetID": "EventID"}
        )

        # Add nodes
        i_GCE_ParsimonyScore_Stats_DF["Parent_Node"] = i_ParentNode_Name
        i_GCE_ParsimonyScore_Stats_DF["Child_Node"]  = i_ChildNode_Name

        # Reorder columns
        first_cols = ["EventID", "Parent_Node", "Child_Node"]
        other_cols = [c for c in i_GCE_ParsimonyScore_Stats_DF.columns if c not in first_cols]
        i_GCE_ParsimonyScore_Stats_DF = i_GCE_ParsimonyScore_Stats_DF[first_cols + other_cols]

        event_stats_list.append(i_GCE_ParsimonyScore_Stats_DF)

    # 7) Concatenate per-site tables
    if per_site_list:
        per_site_all_df = pd.concat(per_site_list, ignore_index=True)
    else:
        per_site_all_df = pd.DataFrame(
            columns=[
                "EventID",
                "Parent_Node",
                "Child_Node",
                "Start_0",
                "Ref",
                "Alt_Alleles",
                "AltAllele_Count",
                "RefAllele_Count",
                "Total_Leaves",
                "Parsimony_Score",
            ]
        )

    # 8) Concatenate event summaries
    if event_stats_list:
        event_summary_all_df = pd.concat(event_stats_list, ignore_index=True)
    else:
        event_summary_all_df = pd.DataFrame(
            columns=[
                "EventID",
                "Parent_Node",
                "Child_Node",
                "N_Sites",
                "Total_Parsimony_Score",
                "Mean_Parsimony_Score",
                "Median_Parsimony_Score",
                "Total_Leaves",
                "Mutations_per_Leaf",
                "Mutations_per_Leaf_per_Site",
            ]
        )

    return per_site_all_df, event_summary_all_df





    



    