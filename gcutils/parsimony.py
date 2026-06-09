# gcutils/parsimony.py

from __future__ import annotations

from typing import Dict, Set, Tuple, Any, Iterable
from typing import Any, Iterable, Optional

from ete3 import Tree
import pandas as pd

from .treeviz_ete_utils import subset_ETE_Tree


# ---------------------------
# Traversal helper
# ---------------------------

def _postorder_nodes(tree: Tree) -> Iterable[Any]:
    """
    Yield nodes in post-order using ETE3.

    Children are visited before their parent, which is what
    the Fitch algorithm requires for the down-pass.
    """
    return tree.traverse("postorder")


# ---------------------------
# Core Fitch scoring (one site)
# ---------------------------

def _fitch_sets_and_score(
    tree: Tree,
    leaf_state: Dict[str, Any],
    *,
    default_state: Any = None,
) -> Tuple[int, Dict[Any, Set[Any]]]:
    """
    Core Fitch algorithm for a single site.

    Parameters
    ----------
    tree : ete3.Tree
        Tree topology (rooting arbitrary; score is root-independent).
    leaf_state : dict
        {leaf_name -> state} mapping for this site.
        States can be any hashable Python objects (str, int, tuple, ...).
    default_state : Any, optional
        State used for leaves missing from leaf_state.

    Returns
    -------
    score : int
        Minimal number of state changes (Fitch parsimony score) for this site.
    R : dict
        {node_obj -> set_of_possible_states} for each node.
        This is mainly useful if you ever want to add reconstruction later.
    """
    R: Dict[Any, Set[Any]] = {}
    score = 0

    for node in _postorder_nodes(tree):
        if node.is_leaf():
            s = leaf_state.get(node.name, default_state)
            R[node] = {s}
        else:
            child_sets = [R[ch] for ch in node.children]
            if not child_sets:
                # Degenerate case (shouldn't happen in a proper tree)
                R[node] = {default_state}
                continue

            running = set(child_sets[0])
            for S in child_sets[1:]:
                inter = running.intersection(S)
                if inter:
                    # intersection non-empty: no extra step
                    running = inter
                else:
                    # intersection empty: we must incur at least one change
                    running = running.union(S)
                    score += 1
            R[node] = running

    return score, R


def fitch_parsimony_score(
    tree: Tree,
    leaf_state: Dict[str, Any],
    *,
    default_state: Any = None,
) -> int:
    """
    Compute the Fitch parsimony score (minimum number of changes) for one site.

    Parameters
    ----------
    tree : ete3.Tree
        Tree topology (rooting arbitrary; score is root-independent).
    leaf_state : dict
        {leaf_name -> state} mapping for this site.
        States can be any hashable Python objects (str, int, tuple, ...).
    default_state : Any, optional
        State used for leaves missing from leaf_state.

    Returns
    -------
    int
        Fitch parsimony score for this site (minimal number of state changes).
    """
    score, _ = _fitch_sets_and_score(tree, leaf_state, default_state=default_state)
    return score





def parsimony_from_df_for_site(
    tree: Tree,
    site_df: pd.DataFrame,
    leaf_col: str,
    state_col: str,
    default_state: Any = 0,
) -> int:
    """
    Compute parsimony score for a single site using a subset DataFrame.

    Parameters
    ----------
    tree : ete3.Tree
        The phylogeny.
    site_df : pd.DataFrame
        Subset of the full variants DF corresponding to a single site.
        Must contain `leaf_col` and `state_col`.
    leaf_col : str
        Column with leaf IDs (matching leaf.name in the tree).
    state_col : str
        Column with the state to use in parsimony (can be string/int/etc.).
    default_state : Any, optional
        State to assign to leaves not present in site_df.

    Returns
    -------
    score : int
        Fitch parsimony score for this site.
    """
    leaf_state_map = dict(zip(site_df[leaf_col], site_df[state_col]))
    return fitch_parsimony_score(tree, leaf_state_map, default_state=default_state)







def parsimony_per_snp(
    tree: Tree,
    variants_df: pd.DataFrame,
    leaf_col: str,
    state_col: str,
    *,
    pos_col: str = "Start_0",
    ref_col: str = "Ref",
    alt_col: str = "Alt",
    snp_col: str = "SNP",
    snp_only: bool = True,
    default_state: Any = 0,
    positions_of_interest: Optional[Iterable[Any]] = None,
) -> pd.DataFrame:
    """
    Compute Fitch parsimony scores for many sites in a variants DataFrame,
    optionally restricted to:
      - leaves present in `tree`
      - a user-provided list of genomic positions.

    Sites are defined by (pos_col, ref_col), which may have multiple
    alternative alleles.

    For each site, defined by a position in `positions_of_interest` (if
    provided) or by positions observed in the filtered variants_df:

      1. Subset rows for that site (site_df) within the subtree leaves.
      2. Build leaf → state map from (leaf_col, state_col).
      3. Compute Fitch parsimony score (0 if no variants => all default_state).
      4. Collect all unique Alt alleles and count alt observations.
      5. Infer Ref from the full variants_df, falling back to site_df.

    Returns
    -------
    DataFrame
        One row per site with columns:
            pos_col
            ref_col
            "Alt_Alleles"
            "AltAllele_Count"
            "RefAllele_Count"
            "Total_Leaves"
            "Parsimony_Score"
    """
    # --- 1) Reference lookup: one ref per position ---
    # Group by position and take the first ref allele so that each pos maps
    # to a *scalar* ref, not a Series.
    if ref_col in variants_df.columns:
        ref_lookup = (
            variants_df[[pos_col, ref_col]]
            .dropna()
            .drop_duplicates()
            .groupby(pos_col)[ref_col]
            .first()
        )
    else:
        ref_lookup = pd.Series(dtype=object)

    # --- 2) Subset variants_df to leaves present in the tree ---
    tree_leaf_ids = set(tree.get_leaf_names())
    df = variants_df[variants_df[leaf_col].isin(tree_leaf_ids)].copy()

    # --- 3) Optionally restrict to SNP rows ---
    if snp_only and snp_col in df.columns:
        df = df[df[snp_col].astype(bool)]

    # --- 4) Optionally filter by positions_of_interest ---
    if positions_of_interest is not None:
        pos_set = set(positions_of_interest)
        df = df[df[pos_col].isin(pos_set)]
        site_positions = sorted(pos_set)
    else:
        site_positions = sorted(df[pos_col].unique())

    # Total number of leaves in the tree (same for all sites)
    total_leaves = len(tree_leaf_ids)

    results = []

    # If there are no positions at all, return an empty DF with expected columns
    if len(site_positions) == 0:
        return pd.DataFrame(
            columns=[
                pos_col,
                ref_col,
                "Alt_Alleles",
                "AltAllele_Count",
                "RefAllele_Count",
                "Total_Leaves",
                "Parsimony_Score",
            ]
        )

    # --- 5) Iterate over each site position ---
    for pos in site_positions:
        # Rows for this position in the filtered DF (subtree leaves only)
        site_df = df[df[pos_col] == pos]

        # Infer reference allele from global lookup or site_df
        if pos in ref_lookup.index:
            ref = ref_lookup.loc[pos]
        elif not site_df.empty and ref_col in site_df.columns:
            ref = site_df[ref_col].iloc[0]
        else:
            ref = pd.NA

        # Safety: if ref somehow ends up as a Series, take the first element
        if isinstance(ref, pd.Series):
            ref = ref.iloc[0]

        # Compute parsimony for this site.
        # If site_df is empty, parsimony_from_df_for_site sees no alt states
        # and all leaves are default_state → score = 0.
        score = parsimony_from_df_for_site(
            tree=tree,
            site_df=site_df,
            leaf_col=leaf_col,
            state_col=state_col,
            default_state=default_state,
        )

        # Collect all unique Alt alleles at this site among target leaves
        if not site_df.empty and alt_col in site_df.columns:
            unique_alts = (
                site_df[alt_col]
                .dropna()
                .astype(str)
                .unique()
                .tolist()
            )
        else:
            unique_alts = []

        alt_alleles_str = ",".join(sorted(unique_alts)) if unique_alts else ""

        AltAllele_Count = site_df.shape[0]
        RefAllele_Count = total_leaves - AltAllele_Count

        results.append(
            {
                pos_col: pos,
                ref_col: ref,
                "Alt_Alleles": alt_alleles_str,
                "AltAllele_Count": AltAllele_Count,
                "RefAllele_Count": RefAllele_Count,
                "Total_Leaves": total_leaves,
                "Parsimony_Score": score,
            }
        )

    result_df = pd.DataFrame(results)
    if not result_df.empty:
        result_df = (
            result_df
            .sort_values([pos_col, ref_col])
            .reset_index(drop=True)
        )

    return result_df




def summarize_parsimony_across_sites(
    per_site_df: pd.DataFrame,
    SetID: str ,
) -> pd.DataFrame:
    """
    Summarize parsimony statistics across many sites.

    Parameters
    ----------
    per_site_df : pd.DataFrame
        Output of `parsimony_per_snp()`. Must contain:
            - "Parsimony_Score"
            - "Total_Leaves"
    SetID : str
        Identifier for the set of variants being summarized. Used to populate
        the `VariantSetID` column in the output.

    Returns
    -------
    summary_df : pd.DataFrame
        A single-row DataFrame with:
            EventID
            N_Sites
            Total_Parsimony_Score
            Mean_Parsimony_Score
            Median_Parsimony_Score
            Total_Leaves
            Score_per_Leaf
            Score_per_Site
            Socre_per_Leaf_per_Site
    """

    if per_site_df.empty:
        cols = [
            "EventID",
            "N_Sites",
            "Total_Parsimony_Score",
            "Mean_Score",
            "Median_Score",
            "Total_Leaves",
            "Score_per_Leaf",
            "Score_per_Site",
            "Socre_per_Leaf_per_Site",
        ]
        return pd.DataFrame([{c: pd.NA for c in cols}])

    # Number of sites
    n_sites = per_site_df.shape[0]

    # Total inferred changes
    total_parsimony = per_site_df["Parsimony_Score"].sum()

    # Mean & median parsimony per site
    mean_parsimony = per_site_df["Parsimony_Score"].mean()
    median_parsimony = per_site_df["Parsimony_Score"].median()

    # Number of leaves (assumed identical for all rows)
    total_leaves = int(per_site_df["Total_Leaves"].unique()[0])

    # Normalized statistics
    mutations_per_leaf = (
        total_parsimony / total_leaves if total_leaves > 0 else pd.NA
    )

    mutations_per_site = (
        total_parsimony / n_sites if n_sites > 0 else pd.NA
    )
    
    mutations_per_leaf_per_site = (
        total_parsimony / (total_leaves * n_sites)
        if (total_leaves > 0 and n_sites > 0)
        else pd.NA
    )

    # 5) number of sites with parsimony score == 0
    # n_sites_parsimony_0 = (per_site_df["Parsimony_Score"] == 0).sum()

    summary = {
        "EventID": SetID,
        "N_Sites": n_sites,
        "Total_Leaves": total_leaves,
        "Total_Score_AllSites": total_parsimony,
        "Mean_Score": mean_parsimony,
        "Median_Score": median_parsimony,
        "Score_per_Leaf": mutations_per_leaf,
        "Score_per_Site": mutations_per_site,
        "Score_per_Leaf_per_Site": mutations_per_leaf_per_site,
    }

    return pd.DataFrame([summary])







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



