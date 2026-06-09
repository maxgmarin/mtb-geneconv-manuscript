
import pandas as pd
import numpy as np

from typing import Optional, Sequence, List, Dict

import bioframe as bf

import re

from .general import label_DF_ByOvrLapGenes, Rv_dist


def parse_HomologyMap_PAF_V1(i_HmMap_PAF):
    
    Mtb_HM_PAF_DF = pd.read_csv(i_HmMap_PAF, sep = "\t", header=None, usecols=np.arange(0,12) )

    # https://github.com/lh3/miniasm/blob/master/PAF.md
    Paf_Cols_1to12 = ["Query_Name", "Query_Len", "Query_Start", "Query_End",
                      "Strand", "Target_Name", "Target_Len", "Target_Start", "Target_End",
                      "Num_ResidueMatches", "Aln_BlockLength", "MapQual"]

    Mtb_HM_PAF_DF.columns = Paf_Cols_1to12

    #Mtb_HM_PAF_DF = Mtb_HM_PAF_DF.drop(["Query_Name", "Query_Len", "Target_Name", "Target_Len", "MapQual"], axis = 1)
    Mtb_HM_PAF_DF = Mtb_HM_PAF_DF.drop(["Query_Len",  "Target_Len", "MapQual"], axis = 1)

    Mtb_HM_PAF_DF["Prop_Match"] = Mtb_HM_PAF_DF["Num_ResidueMatches"] / Mtb_HM_PAF_DF["Aln_BlockLength"]


    Mtb_HM_PAF_DF["Target_Middle"] = (Mtb_HM_PAF_DF["Target_End"] + Mtb_HM_PAF_DF["Target_Start"]) / 2
    Mtb_HM_PAF_DF["Query_Middle"] = (Mtb_HM_PAF_DF["Query_End"] + Mtb_HM_PAF_DF["Query_Start"]) / 2
    Mtb_HM_PAF_DF["Dist_Middles"] = abs(Mtb_HM_PAF_DF["Query_Middle"] - Mtb_HM_PAF_DF["Target_Middle"])

    Mtb_HM_PAF_DF["Query_Length"] = (Mtb_HM_PAF_DF["Query_End"] - Mtb_HM_PAF_DF["Query_Start"])
    Mtb_HM_PAF_DF["Target_Length"] = (Mtb_HM_PAF_DF["Target_End"] - Mtb_HM_PAF_DF["Target_Start"])


    Mtb_HM_PAF_DF = Mtb_HM_PAF_DF.sort_values(["Query_Start", "Query_End", "Target_Start", "Target_End", "Strand"], ascending=True)

    return Mtb_HM_PAF_DF




PAF_COLS_1TO12 = [
    "Query_Name", "Query_Len", "Query_Start", "Query_End",
    "Strand", "Target_Name", "Target_Len", "Target_Start", "Target_End",
    "Num_ResidueMatches", "Aln_BlockLength", "MapQual"
]

def read_mm2_homology_map_paf(
    paf_path: str,
    *,
    drop_cols: Optional[Sequence[str]] = ["Query_Len", "Target_Len", "MapQual"],
    sort_cols: Optional[Sequence[str]] = ["Query_Start", "Query_End", "Target_Start", "Target_End", "Strand"],
) -> pd.DataFrame:
    """
    Read a PAF file and compute common homology-mapping metrics.

    Parameters
    ----------
    paf_path : str
        Path to a (tab-delimited) PAF file (minimap2).
    drop_cols : sequence of str or None, default ("Query_Len", "Target_Len", "MapQual")
        Columns to drop after loading. Set to None to keep all columns.
    sort_cols : sequence of str or None, default (...)
        Columns to sort by. Set to None to skip sorting.

    Returns
    -------
    pandas.DataFrame
        DataFrame with the first 12 PAF fields named per spec, optional drops,
        and extra computed columns:
            - SeqID          = Num_ResidueMatches / Aln_BlockLength
            - Target_Middle       = (Target_End + Target_Start) / 2
            - Query_Middle        = (Query_End + Query_Start) / 2
            - Dist_Middles        = |Query_Middle - Target_Middle|
            - Query_Length        = Query_End - Query_Start
            - Target_Length       = Target_End - Target_Start
    """
    # Dtypes for memory efficiency and correctness
    dtype_map = {
        0: "string",     # Query_Name
        1: "Int64",      # Query_Len
        2: "Int64",      # Query_Start
        3: "Int64",      # Query_End
        4: "string",     # Strand
        5: "string",     # Target_Name
        6: "Int64",      # Target_Len
        7: "Int64",      # Target_Start
        8: "Int64",      # Target_End
        9: "float64",      # Num_ResidueMatches
        10: "float64",     # Aln_BlockLength
        11: "Int64",     # MapQual
    }

    hm_df = pd.read_csv(paf_path,
                        sep="\t",
                        header=None,
                        usecols=range(12),           # read the standard first 12 PAF fields
                        dtype=dtype_map,
                        comment="#",                 # ignore commented lines if any
                    )

    hm_df.columns = PAF_COLS_1TO12
    #print(hm_df.columns)

    # Drop unwanted columns
    hm_df = hm_df.drop(drop_cols, axis = 1)


    hm_df["SeqID"] = hm_df["Num_ResidueMatches"] / hm_df["Aln_BlockLength"]

    q_start, q_end = hm_df["Query_Start"].astype("float64"), hm_df["Query_End"].astype("float64")
    t_start, t_end = hm_df["Target_Start"].astype("float64"), hm_df["Target_End"].astype("float64")

    hm_df["Target_Middle"] = (t_end + t_start) / 2.0
    hm_df["Query_Middle"]  = (q_end + q_start) / 2.0
    hm_df["Dist_Middles"]  = np.abs(hm_df["Query_Middle"] - hm_df["Target_Middle"])

    hm_df["Query_Length"]  = (hm_df["Query_End"] - hm_df["Query_Start"]).astype("Int64")
    hm_df["Target_Length"] = (hm_df["Target_End"] - hm_df["Target_Start"]).astype("Int64")

    # Optional sort
    if sort_cols: hm_df = hm_df.sort_values(list(sort_cols), ascending=True, kind="mergesort")

    return hm_df

# --- Example ---
# paf_df = read_mm2_homology_map_paf("HomMap_k19w19_Trimmed.paf")



PAF_COLS_1TO12 = [
    "Query_Name", "Query_Len", "Query_Start", "Query_End",
    "Strand", "Target_Name", "Target_Len", "Target_Start", "Target_End",
    "Num_ResidueMatches", "Aln_BlockLength", "MapQual"
]

def read_mm2_homology_map_paf_wics(
    paf_path: str,
    *,
    drop_cols: Optional[Sequence[str]] = ("Query_Len", "Target_Len", "MapQual"),
    sort_cols: Optional[Sequence[str]] = ("Query_Start", "Query_End",
                                          "Target_Start", "Target_End", "Strand"),
) -> pd.DataFrame:
    """
    Read a minimap2 PAF homology-map file and compute common metrics,
    while also extracting the cs tag (if present) into a 'cs' column.

    Parameters
    ----------
    paf_path : str
        Path to a (tab-delimited) PAF file.
    drop_cols : sequence of str or None, default ("Query_Len", "Target_Len", "MapQual")
        Columns to drop after loading. Set to None to keep all columns.
    sort_cols : sequence of str or None, default (...)
        Columns to sort by. Set to None to skip sorting.

    Returns
    -------
    pandas.DataFrame
        DataFrame with:
          - First 12 PAF fields named per spec.
          - New column:
              - cs : str or <NA>
                The value portion of the cs tag (everything after "cs:Z:").
          - Extra computed columns:
              - SeqID          = Num_ResidueMatches / Aln_BlockLength
              - Target_Middle  = (Target_End + Target_Start) / 2
              - Query_Middle   = (Query_End + Query_Start) / 2
              - Dist_Middles   = |Query_Middle - Target_Middle|
              - Query_Length   = Query_End - Query_Start
              - Target_Length  = Target_End - Target_Start
    """
    rows: List[Dict] = []

    with open(paf_path, "r") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue

            fields = line.split("\t")
            if len(fields) < 12:
                # Skip malformed lines
                continue

            core = fields[:12]
            aux  = fields[12:]

            # Extract cs tag (if present)
            cs_val = None
            for tag in aux:
                # Optional tags are TAG:TYPE:VALUE, e.g. "cs:Z::50*ac+TT"
                if tag.startswith("cs:"):
                    parts = tag.split(":", 2)
                    if len(parts) == 3:
                        # parts = ["cs", "Z", "<value>"]
                        cs_val = parts[2]
                    else:
                        # Fallback: strip "cs:" and keep whatever remains
                        cs_val = tag[3:]
                    break

            row = dict(zip(PAF_COLS_1TO12, core))
            row["cs"] = cs_val
            rows.append(row)

    if not rows:
        return pd.DataFrame(columns=PAF_COLS_1TO12 + ["cs"])

    hm_df = pd.DataFrame(rows)

    # Dtypes for memory and correctness (mirrors your original function)
    hm_df["Query_Name"]        = hm_df["Query_Name"].astype("string")
    hm_df["Query_Len"]         = hm_df["Query_Len"].astype("Int64")
    hm_df["Query_Start"]       = hm_df["Query_Start"].astype("Int64")
    hm_df["Query_End"]         = hm_df["Query_End"].astype("Int64")
    hm_df["Strand"]            = hm_df["Strand"].astype("string")
    hm_df["Target_Name"]       = hm_df["Target_Name"].astype("string")
    hm_df["Target_Len"]        = hm_df["Target_Len"].astype("Int64")
    hm_df["Target_Start"]      = hm_df["Target_Start"].astype("Int64")
    hm_df["Target_End"]        = hm_df["Target_End"].astype("Int64")
    hm_df["Num_ResidueMatches"] = hm_df["Num_ResidueMatches"].astype("float64")
    hm_df["Aln_BlockLength"]   = hm_df["Aln_BlockLength"].astype("float64")
    hm_df["MapQual"]           = hm_df["MapQual"].astype("Int64")
    hm_df["cs"]                = hm_df["cs"].astype("string")

    # Derived metrics
    hm_df["SeqID"] = hm_df["Num_ResidueMatches"] / hm_df["Aln_BlockLength"]

    q_start = hm_df["Query_Start"].astype("float64")
    q_end   = hm_df["Query_End"].astype("float64")
    t_start = hm_df["Target_Start"].astype("float64")
    t_end   = hm_df["Target_End"].astype("float64")

    hm_df["Target_Middle"] = (t_end + t_start) / 2.0
    hm_df["Query_Middle"]  = (q_end + q_start) / 2.0

    hm_df["Dist_Middles"] = hm_df.apply(
        lambda r: Rv_dist(r["Query_Middle"], r["Target_Middle"]),
        axis=1)
    
    #hm_df["Dist_Middles"]  = (hm_df["Query_Middle"] - hm_df["Target_Middle"]).abs()

    hm_df["Query_Length"]  = (hm_df["Query_End"] - hm_df["Query_Start"]).astype("Int64")
    hm_df["Target_Length"] = (hm_df["Target_End"] - hm_df["Target_Start"]).astype("Int64")

    # Drop requested columns (cs is never dropped by default)
    if drop_cols:
        hm_df = hm_df.drop(list(drop_cols), axis=1, errors="ignore")

    # Optional sort
    if sort_cols:
        hm_df = hm_df.sort_values(list(sort_cols), ascending=True, kind="mergesort")

    return hm_df







def label_HmMap_PAF_DF_ByOvrLapGenes(i_HM_PAF_DF, i_GenomeAnno_Genes_DF):
    
    listOf_TargetOverlap_Genes = []
    listOf_QueryOverlap_Genes = []

    for i, row in i_HM_PAF_DF.iterrows():
        # a) Target overlapping genes
        tar_Start, tar_End = int(row["Target_Start"]), int(row["Target_End"])
        Target_Range = f"NC_000962.3:{tar_Start}-{tar_End}"

        sub_DF_Overlap_Target_Genes = bf.select(i_GenomeAnno_Genes_DF, Target_Range, cols = ("Chrom", "Start", "End"))

        listOf_TargetOverlap_Genes.append( ",".join(list(sub_DF_Overlap_Target_Genes["Symbol"].values)) )

        # b) query overlapping genes
        Query_Start, Query_End = int(row["Query_Start"]), int(row["Query_End"])
        Query_Range = f"NC_000962.3:{Query_Start}-{Query_End}"

        sub_DF_Overlap_Query_Genes = bf.select(i_GenomeAnno_Genes_DF, Query_Range, cols = ("Chrom", "Start", "End"))    

        listOf_QueryOverlap_Genes.append( ",".join(list(sub_DF_Overlap_Query_Genes["Symbol"].values)) )
        
        
    i_HM_PAF_Anno_DF = i_HM_PAF_DF.copy()
    i_HM_PAF_Anno_DF["QueryOverlap_Genes"] = listOf_QueryOverlap_Genes
    i_HM_PAF_Anno_DF["TargetOverlap_Genes"] = listOf_TargetOverlap_Genes


    # For all cases where there are NO OVERLAPPING GENES, simply put "_"
    i_HM_PAF_Anno_DF["QueryOverlap_Genes"] = i_HM_PAF_Anno_DF["QueryOverlap_Genes"].fillna("_")
    i_HM_PAF_Anno_DF["TargetOverlap_Genes"] = i_HM_PAF_Anno_DF["TargetOverlap_Genes"].fillna("_")

    return i_HM_PAF_Anno_DF






def add_self_overlap_column(df: pd.DataFrame) -> pd.DataFrame:
    """
    Add a boolean column 'AlnOverlapsWithSelf' that indicates whether
    the query and target intervals overlap at all.

    Parameters
    ----------
    df : pandas.DataFrame
        DataFrame with PAF-style columns: Query_Start, Query_End,
        Target_Start, Target_End.

    Returns
    -------
    pandas.DataFrame
        Copy of df with new boolean column 'AlnOverlapsWithSelf'.
    """
    df = df.copy()

    q_start, q_end = df["Query_Start"], df["Query_End"]
    t_start, t_end = df["Target_Start"], df["Target_End"]

    # Overlap condition: intervals intersect if neither ends before the other starts
    df["AlnOverlapsWithSelf"] = ~((t_end < q_start) | (t_start > q_end))

    return df

# Example usage
# paf_df = load_paf_with_homology_metrics("HomMap_k19w19_Trimmed.paf")
# paf_df = add_self_overlap_column(paf_df)



def Annotate_HmMap_Aln_DF(i_HmMap_Aln_DF, i_GenomeAnno_Genes_DF, i_Transposase_Genes_DF):
    """
    Process a parsed homology map PAF dataframe by adding gene overlaps,
    self-overlap classification, perfect repeat labeling, and transposase overlap annotation.
    
    Parameters
    ----------
    HmMap_Aln_DF : pd.DataFrame
        Parsed homology map alignment DataFrame (from load_mm2_homology_map_paf).
    H37Rv_GenomeAnno_Genes_DF : pd.DataFrame
        Genome annotation DataFrame for H37Rv.
    Transposase_Genes_DF : pd.DataFrame
        Annotation DataFrame of transposase genes.

    Returns
    -------
    pd.DataFrame
        Processed alignment DataFrame with additional annotations.
    """

    PAF_Query_CoordCols = ("Query_Name", "Query_Start", "Query_End")
    PAF_Target_CoordCols = ("Target_Name", "Target_Start", "Target_End")
    GenomeAnno_CoordCols = ("Chrom", "Start", "End")

    A = i_HmMap_Aln_DF.copy()

    # Step 2: Add info about overlapping genes
    A = label_HmMap_PAF_DF_ByOvrLapGenes(A, i_GenomeAnno_Genes_DF)

    # Step 3: Add column for self overlap
    A = add_self_overlap_column(A)

    # Step 4: Classify perfect repeats
    A["PerfectRepeat"] = np.where(A['SeqID'] == 1.0, True, False)

    # Step 5: Count overlaps with transposases
    A = bf.count_overlaps(
        A, i_Transposase_Genes_DF,
        cols1 = PAF_Query_CoordCols,
        cols2 = GenomeAnno_CoordCols
    ).rename(columns={'count': 'Query_NOvrlap_Transposase'})

    A = bf.count_overlaps(
        A, i_Transposase_Genes_DF,
        cols1 = PAF_Target_CoordCols,
        cols2 = GenomeAnno_CoordCols
    ).rename(columns={'count': 'Target_NOvrlap_Transposase'})
    
    # Step 6: infer query and target coords and then create query-to-target ID
    A["QueryCoords"]      = A["Query_Name"] + ":" + A["Query_Start"].astype(str) + "-" + A["Query_End"].astype(str)
    A["TargetCoords"]     = A["Target_Name"] + ":" + A["Target_Start"].astype(str) + "-" + A["Target_End"].astype(str)    
    A["QueryToTarget_ID"] = A["QueryCoords"] + ";" + A["TargetCoords"] 

    return A



def load_and_process_mm2_homology_map_paf(paf_path,
                                          i_GenomeAnno_Genes_DF,
                                          i_Transposase_Genes_DF):


    # i_HmMap_Aln_DF = read_mm2_homology_map_paf(paf_path)
    i_HmMap_Aln_DF = read_mm2_homology_map_paf_wics(paf_path)


    i_HmMap_Aln_WiAnno_DF = Annotate_HmMap_Aln_DF(i_HmMap_Aln_DF,
                                                  i_GenomeAnno_Genes_DF,
                                                  i_Transposase_Genes_DF)





    return i_HmMap_Aln_WiAnno_DF





def annotate_homology_aln_DF_by_HmRegionNum(
    hm_pair_df,
    hm_regions_df,
    verbose = True
) -> pd.DataFrame:
    """
    Annotate homologous region pairs with region IDs and overlapping gene info.

    Parameters
    ----------
    hm_pair_df : pd.DataFrame
        DataFrame containing homology pairs with columns 'QueryCoords' and 'TargetCoords'.
    hm_regions_df : pd.DataFrame
        DataFrame containing homologous region information (must include
        'HmRegion_Num' and 'Overlap_Genes').
    bf : module or object
        Object exposing a `select(df, region, cols)` function.
    verbose : bool, optional
        If True, prints warnings when multiple matches are found.

    Returns
    -------
    pd.DataFrame
        Updated DataFrame with 'Query_HmRegionID', 'Target_HmRegionID',
        'Query_Overlap_Genes', and 'Target_Overlap_Genes' added.
    """
    updated_rows = []

    for idx, row in hm_pair_df.iterrows():
        query_region = row["QueryCoords"]
        target_region = row["TargetCoords"]

        q_regions_df = bf.select(hm_regions_df, query_region,
                        cols=["Chr", "Start", "End"])
        t_regions_df = bf.select(hm_regions_df, target_region,
                        cols=["Chr", "Start", "End" ])

        if q_regions_df.shape[0] != 1 and verbose:
            print(f"BIG ERROR: Query region {query_region} matched {q_regions_df.shape[0]} rows")
        if t_regions_df.shape[0] != 1 and verbose:
            print(f"BIG ERROR: Target region {target_region} matched {t_regions_df.shape[0]} rows")

        row = row.copy()  # avoid mutating the iterator view
        row["Query_HmRegionID"]    = q_regions_df["HmRegionID"].iat[0]
        row["Target_HmRegionID"]   = t_regions_df["HmRegionID"].iat[0]
        row["Query_HmRegionNum"]    = q_regions_df["HmRegion_Num"].iat[0]
        row["Target_HmRegionNum"]   = t_regions_df["HmRegion_Num"].iat[0]
#        row["Target_Overlap_Genes"] = t_regions_df["Overlap_Genes"].iat[0]
#        row["Query_Overlap_Genes"] = q_regions_df["Overlap_Genes"].iat[0]

        updated_rows.append(row)

    return pd.DataFrame(updated_rows)



def build_ungapped_homology_regions(i_HmPair_DF, i_GenomeAnno_Genes_DF, i_Transposase_Genes_DF, genome_chr="NC_000962.3", regionPrefix = ""):
    """
    Build merged ungapped homology regions from paired homology-map alignments.

    Parameters
    ----------
    i_HmPair_DF : pd.DataFrame
        Input homology-map pairs (e.g., HmMap_Aln_k19w19_NoOverlap_DF).
    i_GenomeAnno_Genes_DF : pd.DataFrame
        Genome annotation table (genes) for overlap labeling.
    i_Transposase_Genes_DF : pd.DataFrame
        Transposase annotations for overlap counting.
    genome_chr : str, optional
        Chromosome name to assign in output, by default "NC_000962.3".

    Returns
    -------
    pd.DataFrame
        UngHmRegions_Merged_DF with merged ranges, gene overlap labels,
        transposase overlap counts, and region IDs.
    """

    # Step 7: Cluster all OVERLAPPING query and target alignments into a single region
    HmCluster_Query_Cols =  ["Query_Name", "Query_Start", "Query_End"]
    HmCluster_Target_Cols = ["Target_Name", "Target_Start", "Target_End"]
    GenomeAnno_CoordCols = ("Chrom", "Start", "End")

    i_HmPair_DF = bf.cluster(i_HmPair_DF,
                   cols = HmCluster_Query_Cols,
                   min_dist=0).rename(columns={'cluster': 'Query_Cluster',
                                               'cluster_start': 'Query_Cluster_Start',
                                               'cluster_end': 'Query_Cluster_End'})

    i_HmPair_DF = bf.cluster(i_HmPair_DF,
                                cols = HmCluster_Target_Cols,
                                min_dist=0).rename(columns={'cluster': 'Target_Cluster',
                                                            'cluster_start': 'Target_Cluster_Start',
                                                            'cluster_end': 'Target_Cluster_End'})

    
    # --- Derive unique target and query ranges ---
    i_Hm_All_TargetRanges_DF = (
        i_HmPair_DF[["Target_Cluster", "Target_Cluster_Start", "Target_Cluster_End"]]
        .drop_duplicates()
        .sort_values(["Target_Cluster", "Target_Cluster_Start", "Target_Cluster_End"])
        .reset_index(drop=True).rename(columns={'Target_Cluster_Start': 'R1_cluster_start',
                                                'Target_Cluster_End': 'R1_cluster_end'})
    )
    
    i_Hm_All_TargetRanges_DF["Chr"] = genome_chr
    
    i_Hm_All_QueryRanges_DF = (
        i_HmPair_DF[["Query_Cluster", "Query_Cluster_Start", "Query_Cluster_End"]]
        .drop_duplicates()
        .sort_values(["Query_Cluster", "Query_Cluster_Start", "Query_Cluster_End"])
        .reset_index(drop=True).rename(columns={'Query_Cluster_Start': 'R1_cluster_start',
                                                'Query_Cluster_End': 'R1_cluster_end'})
    )
    
    i_Hm_All_QueryRanges_DF["Chr"] = genome_chr
    
    
    #print("Shape of i_Hm_All_TargetRanges_DF:", i_Hm_All_TargetRanges_DF.shape)
    #print("Shape of i_Hm_All_QueryRanges_DF:",  i_Hm_All_QueryRanges_DF.shape)
    
    
    # --- Merge query & target ranges on cluster id ---
    Hm_MergedQamdT_Ranges_DF = pd.concat([i_Hm_All_QueryRanges_DF, i_Hm_All_TargetRanges_DF],
                                         ignore_index=True)
    
    #print("shape of Hm_MergedQamdT_Ranges_DF:", Hm_MergedQamdT_Ranges_DF.shape)
    #print("columns of Hm_MergedQamdT_Ranges_DF:", Hm_MergedQamdT_Ranges_DF.columns)
    
    
    Hm_MergedQamdT_Ranges_V2_DF = bf.cluster(Hm_MergedQamdT_Ranges_DF,
                             cols = ["Chr", "R1_cluster_start", 'R1_cluster_end'],
                             min_dist=0).rename(columns={'cluster': 'HmRegion_Num',
                                                         'cluster_start': 'Start',
                                                         'cluster_end': 'End'})
    
    
    #print("shape of Hm_MergedQamdT_Ranges_V2_DF:", Hm_MergedQamdT_Ranges_DF.shape)
    #print("columns of Hm_MergedQamdT_Ranges_V2_DF:", Hm_MergedQamdT_Ranges_DF.columns)
    

    
    UngHmRegions_Merged_DF = Hm_MergedQamdT_Ranges_V2_DF[["HmRegion_Num", "Chr", "Start", "End"]].drop_duplicates(
                                                         ).sort_values(["HmRegion_Num", "Start", "End"]).reset_index(drop=True)
    
    
    #print("shape of UngHmRegions_Merged_DF:", UngHmRegions_Merged_DF.shape)
    #print("columns of UngHmRegions_Merged_DF:", UngHmRegions_Merged_DF.columns)
    
    
    # --- Build the Ungapped Homology Regions table ---
    #UngHmRegions_Merged_DF["Chr"] = genome_chr
    UngHmRegions_Merged_DF["Center"] = (UngHmRegions_Merged_DF["Start"] + UngHmRegions_Merged_DF["End"]) / 2
    UngHmRegions_Merged_DF["Length"] = UngHmRegions_Merged_DF["End"] - UngHmRegions_Merged_DF["Start"]
    
    
    UngHmRegions_Merged_DF = UngHmRegions_Merged_DF.sort_values(["Chr", "Start", "End"])
    


    # --- Build the Ungapped Homology Regions table ---
    #UngHmRegions_Merged_DF["Chr"] = genome_chr
    UngHmRegions_Merged_DF["Center"] = (UngHmRegions_Merged_DF["Start"] + UngHmRegions_Merged_DF["End"]) / 2
    UngHmRegions_Merged_DF["Length"] = UngHmRegions_Merged_DF["End"] - UngHmRegions_Merged_DF["Start"]

    # Label by overlapping genes
    UngHmRegions_Merged_DF = label_DF_ByOvrLapGenes(
        UngHmRegions_Merged_DF,
        i_GenomeAnno_Genes_DF,
        i_cols1=("Chr", "Start", "End")
    )

    # Count overlaps with transposases
    UngHmRegions_Merged_DF = bf.count_overlaps(
        UngHmRegions_Merged_DF,
        i_Transposase_Genes_DF,
        cols1=("Chr", "Start", "End"),
        cols2=GenomeAnno_CoordCols
    ).rename(columns={"count": "NOvrlap_Transposase"})


    i_HmPair_NoPerfAln_DF = i_HmPair_DF.query("SeqID != 1.0")

    i_HmPair_MaxSeqID99_DF = i_HmPair_DF.query("SeqID <= 0.99")


    UngHmRegions_Merged_DF = bf.count_overlaps(UngHmRegions_Merged_DF,
                                               i_HmPair_DF, # ALL PAIRS of homology alignments, overlap and non-overlap
                                               cols1 = ("Chr", "Start", "End"),
                                               cols2 = HmCluster_Query_Cols ).rename(columns={'count': 'num_HmAln_All'})


    UngHmRegions_Merged_DF = bf.count_overlaps(UngHmRegions_Merged_DF,
                                               i_HmPair_NoPerfAln_DF,
                                               cols1 = ("Chr", "Start", "End"),
                                               cols2 = HmCluster_Query_Cols ).rename(columns={'count': 'num_HmAln_PR_BelowSeqID100'})

    UngHmRegions_Merged_DF = bf.count_overlaps(UngHmRegions_Merged_DF,
                                               i_HmPair_MaxSeqID99_DF, 
                                               cols1 = ("Chr", "Start", "End"),
                                               cols2 = HmCluster_Query_Cols ).rename(columns={'count': 'num_HmAln_PR_MinSeqID99'})



    # Assign unique IDs
    UngHmRegions_Merged_DF["HmRegionNum"] = UngHmRegions_Merged_DF.index  # + 1 if you want 1-based
    UngHmRegions_Merged_DF["HmRegionID"] = regionPrefix + "HmRegion_" + UngHmRegions_Merged_DF["HmRegionNum"].astype(str).str.rjust(3, "0")



    # Add unique Region IDs that link each pair of aligned regions (Query and Target region)
    i_HmPair_DF = annotate_homology_aln_DF_by_HmRegionNum(i_HmPair_DF, UngHmRegions_Merged_DF)

    return UngHmRegions_Merged_DF, i_HmPair_DF










###### Functions for variant inference from cs tags associated w/ each homology-map alignment ######

CS_OP_RE = re.compile(r'([:=*+-])(\d+|[A-Za-z]+)')


def infer_variants_from_PAF_per_row(row: pd.Series) -> pd.DataFrame:
    """
    Decode the minimap2 'cs' string for a single homology-map alignment and
    return a per-variant DataFrame with coordinates in both target (reference)
    and query.

    Expected columns in `row`:
      - Query_Name, Query_Start, Query_End
      - Target_Name, Target_Start, Target_End
      - Strand ('+' or '-')
      - cs (string)

    Returns
    -------
    DataFrame with columns (in this order):
      - Target_Name, Target_Start, Target_End
      - Strand, Ref, Alt, SNP, Type
      - Query_Name, Query_Start, Query_End
      - HmMap_Aln_ID, Aln_Query_Start, Aln_Query_End
      - Aln_Target_Start, Aln_Target_End, Aln_Strand
    """
    cols_order = [
        "Target_Name","Target_Start","Target_End",
        "Strand","Ref","Alt","SNP","Type",
        "Query_Name","Query_Start","Query_End",
        "HmMap_Aln_ID","Aln_Query_Start","Aln_Query_End",
        "Aln_Target_Start","Aln_Target_End","Aln_Strand"
    ]

    cs = row["cs"]
    # Gracefully handle missing cs (no variants can be decoded)
    if pd.isna(cs) or cs is None:
        return pd.DataFrame(columns=cols_order)

    ref_name   = row["Target_Name"]
    ref_start  = int(row["Target_Start"])
    ref_end    = int(row["Target_End"])
    query_name = row["Query_Name"]
    query_start = int(row["Query_Start"])
    query_end   = int(row["Query_End"])
    strand = "+" if str(row["Strand"]) == "+" else "-"

    # Walk coordinates along target (x) and query (y)
    x = ref_start
    if strand == "+":
        y = query_start
    else:
        # For reverse alignments, we traverse query backwards
        y = query_end

    #aln_id = f"{query_name}:{query_start}-{query_end}:{strand}->{ref_name}:{ref_start}-{ref_end}"
    aln_id = f"{query_name}:{query_start}-{query_end};{ref_name}:{ref_start}-{ref_end}"


    variants = []

    for m in CS_OP_RE.finditer(cs):
        op  = m.group(1)
        arg = m.group(2)

        # Matches (no variant)
        if op == '=':
            l = len(arg)
            if strand == "+":
                y += l
            else:
                y -= l
            x += l

        elif op == ':':
            l = int(arg)
            if strand == "+":
                y += l
            else:
                y -= l
            x += l

        # Substitution: *xy where x = ref base, y = query base
        elif op == '*':
            ref_base = arg[0]
            alt_base = arg[1]

            # Query coords depend on strand
            if strand == "+":
                qs = y
                qe = y + 1
                y += 1
            else:
                qs = y - 1
                qe = y
                y -= 1

            rs = x
            re_ = x + 1
            x += 1

            # Skip ambiguous bases but still advance coords
            if ref_base.lower() == 'n' or alt_base.lower() == 'n':
                continue

            variants.append({
                "Target_Name": ref_name,
                "Target_Start": rs,
                "Target_End": re_,
                "Strand": strand,
                "Ref": ref_base.upper(),
                "Alt": alt_base.upper(),
                "Query_Name": query_name,
                "Query_Start": qs,
                "Query_End": qe,
                "HmMap_Aln_ID": aln_id,
                "Aln_Query_Start": query_start,
                "Aln_Query_End": query_end,
                "Aln_Target_Start": ref_start,
                "Aln_Target_End": ref_end,
                "Aln_Strand": strand,
            })

        # Insertion in query (relative to target): +seq
        elif op == '+':
            ins_seq = arg.upper()
            l = len(ins_seq)

            if strand == "+":
                qs = y
                qe = y + l
                y += l
            else:
                qs = y - l
                qe = y
                y -= l

            rs = x
            re_ = x  # insertion sits between x-1 and x

            variants.append({
                "Target_Name": ref_name,
                "Target_Start": rs,
                "Target_End": re_,
                "Strand": strand,
                "Ref": "-",
                "Alt": ins_seq,
                "Query_Name": query_name,
                "Query_Start": qs,
                "Query_End": qe,
                "HmMap_Aln_ID": aln_id,
                "Aln_Query_Start": query_start,
                "Aln_Query_End": query_end,
                "Aln_Target_Start": ref_start,
                "Aln_Target_End": ref_end,
                "Aln_Strand": strand,
            })

        # Deletion in query (bases present in ref only): -seq
        elif op == '-':
            del_seq = arg.upper()
            l = len(del_seq)

            rs = x
            re_ = x + l
            x += l

            qs = y
            qe = y  # no query bases

            variants.append({
                "Target_Name": ref_name,
                "Target_Start": rs,
                "Target_End": re_,
                "Strand": strand,
                "Ref": del_seq,
                "Alt": "-",
                "Query_Name": query_name,
                "Query_Start": qs,
                "Query_End": qe,
                "HmMap_Aln_ID": aln_id,
                "Aln_Query_Start": query_start,
                "Aln_Query_End": query_end,
                "Aln_Target_Start": ref_start,
                "Aln_Target_End": ref_end,
                "Aln_Strand": strand,
            })

        else:
            raise ValueError(f"Unknown cs op {op!r} in cs={cs!r}")

    if not variants:
        return pd.DataFrame(columns=cols_order)

    var_df = pd.DataFrame(variants)

    # SNP flag, then a simple type
    var_df["SNP"] = (
        (var_df["Ref"].str.len() == 1) &
        (var_df["Alt"].str.len() == 1) &
        (var_df["Ref"] != '-') &
        (var_df["Alt"] != '-')
    )

    def classify(r):
        if r["SNP"]:
            return "SNP"
        if r["Ref"] == "-" and r["Alt"] != "-":
            return "INS"
        if r["Alt"] == "-" and r["Ref"] != "-":
            return "DEL"
        return "OTHER"

    var_df["Type"] = var_df.apply(classify, axis=1)

    # Reorder columns as requested
    var_df = var_df[cols_order]
    return var_df


def get_variants_from_PAF_Aln_df(hm_df: pd.DataFrame) -> pd.DataFrame:
    """
    Apply infer_variants_from_PAF_per_row() to every alignment in a homology
    map DataFrame and concatenate the per-alignment variant tables.

    Returns a DataFrame with the same column order as
    infer_variants_from_PAF_per_row().
    """
    cols_order = [
        "Target_Name","Target_Start","Target_End",
        "Strand","Ref","Alt","SNP","Type",
        "Query_Name","Query_Start","Query_End",
        "HmMap_Aln_ID","Aln_Query_Start","Aln_Query_End",
        "Aln_Target_Start","Aln_Target_End","Aln_Strand"
    ]

    dfs = []
    for _, row in hm_df.iterrows():
        vdf = infer_variants_from_PAF_per_row(row)
        if not vdf.empty:
            dfs.append(vdf)

    if not dfs:
        return pd.DataFrame(columns=cols_order)

    out = pd.concat(dfs, ignore_index=True)
    # Just in case, enforce column order
    out = out[cols_order]
    return out


def label_HmMap_Variants_DF_ByOvrLapGenes(i_HmAln_Variants_DF, i_GenomeAnno_Genes_DF):
    
    listOf_TargetOverlap_Genes = []
    listOf_QueryOverlap_Genes = []

    for i, row in i_HmAln_Variants_DF.iterrows():
        # a) Target overlapping genes
        tar_Start, tar_End = int(row["Aln_Target_Start"]), int(row["Aln_Target_End"])
        Target_Range = f"NC_000962.3:{tar_Start}-{tar_End}"

        sub_DF_Overlap_Target_Genes = bf.select(i_GenomeAnno_Genes_DF, Target_Range, cols = ("Chrom", "Start", "End"))

        listOf_TargetOverlap_Genes.append( ",".join(list(sub_DF_Overlap_Target_Genes["Symbol"].values)) )

        # b) query overlapping genes
        Query_Start, Query_End = int(row["Aln_Query_Start"]), int(row["Aln_Query_End"])
        Query_Range = f"NC_000962.3:{Query_Start}-{Query_End}"

        sub_DF_Overlap_Query_Genes = bf.select(i_GenomeAnno_Genes_DF, Query_Range, cols = ("Chrom", "Start", "End"))    

        listOf_QueryOverlap_Genes.append( ",".join(list(sub_DF_Overlap_Query_Genes["Symbol"].values)) )
        
        
    HmAln_Var_V2_DF = i_HmAln_Variants_DF.copy()
    HmAln_Var_V2_DF["QueryOverlap_Genes"] = listOf_QueryOverlap_Genes
    HmAln_Var_V2_DF["TargetOverlap_Genes"] = listOf_TargetOverlap_Genes

    # For all cases where there are NO OVERLAPPING GENES, simply put "_"
    HmAln_Var_V2_DF["QueryOverlap_Genes"] = HmAln_Var_V2_DF["QueryOverlap_Genes"].fillna("_")
    HmAln_Var_V2_DF["TargetOverlap_Genes"] = HmAln_Var_V2_DF["TargetOverlap_Genes"].fillna("_")

    HmAln_Var_V2_DF["QueryParalog_RegionID"] = HmAln_Var_V2_DF["QueryOverlap_Genes"] +"-"+ HmAln_Var_V2_DF["Query_Name"] +":"+ HmAln_Var_V2_DF["Aln_Query_Start"].astype(str) +"-"+ HmAln_Var_V2_DF["Aln_Query_End"].astype(str)

    return HmAln_Var_V2_DF


################################################################################################







