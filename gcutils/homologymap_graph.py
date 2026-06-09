
import pandas as pd
import numpy as np
import networkx as nx

# from typing import Optional, Sequence
# import bioframe as bf

# from .general import label_DF_ByOvrLapGenes


def _to_float(x, default=0.0):
    try:
        if pd.isna(x):
            return default
        return float(x)
    except Exception:
        return default



## Define functions for building a graph all of Paralogous Regions in the Mtb genome
def CreatePR_DiGraph(
    hm_regions_df,
    hm_pairs_df,
    edge_keep_mode="keep_all",  # default behavior
):
    """
    Create a directed graph of homologous regions (PRs).

    Nodes:
        One node per region in hm_regions_df; isolated nodes removed at end.

    Edges:
        Directed edges from Query_HmRegionNum → Target_HmRegionNum.

    edge_keep_mode:
        - "keep_all": keep every edge in hm_pairs_df (even if counts are 0)
        - "keep_edges_Any_GCE": keep edges where N_Ovrlap_pGCE_ToTarget > 0
        - "keep_edges_Mapped_GCE": keep edges where EventsMapped > 0
    """

    valid_modes = {"keep_all", "keep_edges_Any_GCE", "keep_edges_Mapped_GCE"}
    if edge_keep_mode not in valid_modes:
        raise ValueError(f"edge_keep_mode must be one of {sorted(valid_modes)}; got {edge_keep_mode!r}")

    # Build mapping RegionNum → gene(s)
    RegionNum_To_Genes = (
        hm_regions_df
        .set_index("HmRegion_Num")["Overlap_Genes"]
        .astype(str)
        .to_dict()
    )

    G = nx.DiGraph()

    has_pr       = "PR_SetID" in hm_regions_df.columns
    has_mapped   = "N_Events_Mapped" in hm_regions_df.columns
    has_putative = "N_Events_Putative" in hm_regions_df.columns

    # ------------------
    # Add nodes
    # ------------------
    for _, row in hm_regions_df.iterrows():
        region_num        = row["HmRegion_Num"]
        overlapping_genes = row["Overlap_Genes"]

        num_events_mGCE = _to_float(row["N_Events_Mapped"], default=0.0) if has_mapped else 0.0
        num_events_pGCE = _to_float(row["N_Events_Putative"], default=0.0) if has_putative else 0.0

        node_name = f"{region_num} - {overlapping_genes}"

        node_attrs = dict(
            region_idnum=region_num,
            genes=overlapping_genes,
            label=overlapping_genes,
            counts_mGCE=num_events_mGCE,
            counts_pGCE=num_events_pGCE,
        )
        if has_pr:
            node_attrs["PR_SetID"] = row["PR_SetID"]

        G.add_node(node_name, **node_attrs)

    # ------------------
    # Add directed edges
    # ------------------
    for _, row in hm_pairs_df.iterrows():
        q = row["Query_HmRegionNum"]
        t = row["Target_HmRegionNum"]

        qnode = f"{q} - {RegionNum_To_Genes.get(q, 'None')}"
        tnode = f"{t} - {RegionNum_To_Genes.get(t, 'None')}"

        if qnode not in G.nodes or tnode not in G.nodes:
            print(f"WARNING: missing node(s) for edge {qnode} -> {tnode}")
            # continue  # uncomment if you want to skip missing endpoints

        edge_N_GCEs_Mapped = _to_float(row.get("EventsMapped", 0.0), default=0.0)
        edge_N_GCEs_Total  = _to_float(row.get("N_Ovrlap_pGCE_ToTarget", 0.0), default=0.0)

        # decide whether to keep this edge
        if edge_keep_mode == "keep_all":
            keep_edge = True
        elif edge_keep_mode == "keep_edges_Any_GCE":
            keep_edge = (edge_N_GCEs_Total > 0)
        elif edge_keep_mode == "keep_edges_Mapped_GCE":
            keep_edge = (edge_N_GCEs_Mapped > 0)

        if keep_edge:
            G.add_edge(
                qnode, tnode,
                weight=edge_N_GCEs_Mapped,
                label=str(edge_N_GCEs_Mapped),
                GCEs_Mapped=edge_N_GCEs_Mapped,
                GCEs_Total=edge_N_GCEs_Total,
            )

    # Remove isolated nodes (same behavior as before)
    G.remove_nodes_from(list(nx.isolates(G)))

    return G



def CreatePR_UGraph(
    hm_regions_df,
    hm_pairs_df,
    edge_keep_mode="keep_all",
):
    """
    Create an undirected graph of homologous regions (PRs).

    Nodes come from hm_regions_df. Edges are region-pairs from hm_pairs_df.
    Multiple rows between the same two nodes are collapsed into a single edge
    by summing the edge attributes.

    edge_keep_mode:
        - "keep_all": keep every edge row
        - "keep_edges_Any_GCE": keep rows where N_Ovrlap_pGCE_ToTarget > 0
        - "keep_edges_Mapped_GCE": keep rows where EventsMapped > 0

    Returns:
        G : nx.Graph (undirected, no isolates)
    """
    valid_modes = {"keep_all", "keep_edges_Any_GCE", "keep_edges_Mapped_GCE"}
    if edge_keep_mode not in valid_modes:
        raise ValueError(f"edge_keep_mode must be one of {sorted(valid_modes)}; got {edge_keep_mode!r}")

    # Build mapping RegionNum → GeneName(s)
    RegionNum_To_Genes = (
        hm_regions_df
        .set_index("HmRegion_Num")["Overlap_Genes"]
        .astype(str)
        .to_dict()
    )

    G = nx.Graph()

    has_pr       = "PR_SetID" in hm_regions_df.columns
    has_mapped   = "N_Events_Mapped" in hm_regions_df.columns
    has_putative = "N_Events_Putative" in hm_regions_df.columns

    # ------------------
    # Add nodes
    # ------------------
    for _, row in hm_regions_df.iterrows():
        region_num        = row["HmRegion_Num"]
        overlapping_genes = row["Overlap_Genes"]

        num_events_mGCE = _to_float(row["N_Events_Mapped"], default=0.0) if has_mapped else 0.0
        num_events_pGCE = _to_float(row["N_Events_Putative"], default=0.0) if has_putative else 0.0

        node_name = f"{region_num} - {overlapping_genes}"

        node_attrs = dict(
            region_idnum=region_num,
            genes=overlapping_genes,
            label=overlapping_genes,
            counts_mGCE=num_events_mGCE,
            counts_pGCE=num_events_pGCE,
        )
        if has_pr:
            node_attrs["PR_SetID"] = row["PR_SetID"]

        G.add_node(node_name, **node_attrs)

    # ------------------
    # Add edges (collapse multi-rows by summing)
    # ------------------
    for _, row in hm_pairs_df.iterrows():
        q = row["Query_HmRegionNum"]
        t = row["Target_HmRegionNum"]

        qnode = f"{q} - {RegionNum_To_Genes.get(q, 'None')}"
        tnode = f"{t} - {RegionNum_To_Genes.get(t, 'None')}"

        if qnode not in G.nodes or tnode not in G.nodes:
            print(f"WARNING: missing node(s) for edge {qnode} -- {tnode}")
            # continue

        GCEs_Mapped = _to_float(row.get("EventsMapped", 0.0), default=0.0)
        GCEs_Total  = _to_float(row.get("N_Ovrlap_pGCE_ToTarget", 0.0), default=0.0)

        # decide whether to keep this *row* as evidence for an undirected edge
        if edge_keep_mode == "keep_all":
            keep_row = True
        elif edge_keep_mode == "keep_edges_Any_GCE":
            keep_row = (GCEs_Total > 0)
        else:  # "keep_edges_Mapped_GCE"
            keep_row = (GCEs_Mapped > 0)

        if not keep_row:
            continue

        # collapse into an undirected edge by summing attributes
        if G.has_edge(qnode, tnode):
            d = G[qnode][tnode]
            d["GCEs_Mapped"] = d.get("GCEs_Mapped", 0.0) + GCEs_Mapped
            d["GCEs_Total"]  = d.get("GCEs_Total",  0.0) + GCEs_Total

            # choose what "weight" means for the undirected graph:
            # here: weight tracks mapped GCEs (consistent w/ DiGraph)
            d["weight"] = d["GCEs_Mapped"]
            d["label"]  = str(d["weight"])
        else:
            G.add_edge(
                qnode,
                tnode,
                GCEs_Mapped=GCEs_Mapped,
                GCEs_Total=GCEs_Total,
                weight=GCEs_Mapped,
                label=str(GCEs_Mapped),
            )

    # Remove isolated nodes
    G.remove_nodes_from(list(nx.isolates(G)))

    return G




# Let's create a function that maps each HmRegion ID to its connected set of nodes (i.e., components) in the graph.
# We will then produce both a table and a dictionary showing the mapping of each HmRegion ID to its connected component (PR_SetID).

def infer_HmMap_ParalogSet_Labels(graph):
    """
    Return:
        pr_set_dict: {region_idnum → PR_Set_#}
        pr_set_df:   DataFrame with columns ['HmRegion_ID', 'PR_SetID']
    """

    undirected = graph.to_undirected()
    components = list(nx.connected_components(undirected))

    pr_set_dict = {}

    for i, component in enumerate(components, start=1):
        pr_set_id = f"PR_Set_{i}"

        for node in component:
            region_id = graph.nodes[node].get("region_idnum")

            if region_id is None:
                raise ValueError(
                    f"Node {node} missing required attribute 'region_idnum'."
                )

            pr_set_dict[region_id] = pr_set_id

    pr_set_df = (
        pd.DataFrame
        .from_dict(pr_set_dict, orient='index', columns=["PR_SetID"])
        .rename_axis("HmRegion_ID")
        .reset_index()
    )

    return pr_set_dict, pr_set_df





def add_PR_SetID_to_graph(graph, pr_set_dict):
    """
    Return a COPY of the graph where each node has a 'PR_SetID' attribute
    assigned based on the region_idnum → PR_SetID mapping.

    Args:
        graph (nx.Graph or nx.DiGraph):
            The input graph whose nodes contain a 'region_idnum' attribute.
        pr_set_dict (dict):
            Mapping {region_idnum → PR_SetID} returned by
            infer_HmMap_ParalogSet_Labels().

    Returns:
        new_graph (nx.Graph or nx.DiGraph):
            A copy of the input graph with PR_SetID added to node attributes.
    """

    # Copy the graph so original is not modified
    new_graph = graph.copy()

    for node, attrs in new_graph.nodes(data=True):
        region_id = attrs.get("region_idnum", None)

        # Assign PR_SetID (or None if missing)
        new_graph.nodes[node]["PR_SetID"] = pr_set_dict.get(region_id, None)

    return new_graph



def subset_graph_by_PR_SetIDs(graph, pr_set_list):
    """
    Return a new graph containing only nodes whose PR_SetID attribute
    is in the provided list of PR_SetIDs.

    Args:
        graph: NetworkX graph (directed or undirected)
        pr_set_list: list or set of PR_SetIDs to keep
                     e.g. ["PR_Set_1", "PR_Set_3"]

    Returns:
        new_graph: A new NetworkX graph object containing only the
                   nodes (and their induced edges) belonging to the
                   specified paralog sets.
    """

    # Cast to set for fast membership lookups
    pr_set_list = set(pr_set_list)

    # Identify nodes to keep
    selected_nodes = [
        node for node, attrs in graph.nodes(data=True)
        if attrs.get("PR_SetID") in pr_set_list
    ]

    # Build and return the induced subgraph (copy=True gives new object)
    new_graph = graph.subgraph(selected_nodes).copy()

    return new_graph





def summarize_graph_stats(G):
    """
    Compute general descriptive statistics for a NetworkX graph (directed or undirected).

    Returns:
        stats (dict): dictionary of summary statistics
        component_sizes_df (DataFrame): two-column table:
            ['ComponentID', 'NumNodes']
    """

    is_directed = G.is_directed()

    # Basic stats
    num_nodes = G.number_of_nodes()
    num_edges = G.number_of_edges()

    # Degrees
    degrees = dict(G.degree())
    degree_vals = np.array(list(degrees.values()))

    # Component structure
    if is_directed:
        weak_components  = list(nx.weakly_connected_components(G))
        strong_components = list(nx.strongly_connected_components(G))
        components = weak_components  # report weak as default
    else:
        components = list(nx.connected_components(G))

    component_sizes = [len(c) for c in components]

    component_sizes_df = pd.DataFrame({
        "ComponentID": [f"Comp_{i+1}" for i in range(len(component_sizes))],
        "NumNodes": component_sizes
    })

    # Edge weight stats (if weights exist)
    weights = []
    for _, _, data in G.edges(data=True):
        if "weight" in data:
            weights.append(data["weight"])

    if len(weights) > 0:
        weights = np.array(weights)
        weight_stats = {
            "edge_weight_min": float(weights.min()),
            "edge_weight_max": float(weights.max()),
            "edge_weight_mean": float(weights.mean()),
            "edge_weight_median": float(np.median(weights)),
        }
    else:
        weight_stats = {
            "edge_weight_min": None,
            "edge_weight_max": None,
            "edge_weight_mean": None,
            "edge_weight_median": None,
        }

    stats = {
        "num_nodes": num_nodes,
        "num_edges": num_edges,

        "is_directed": is_directed,

        # Connected components
        "num_connected_components": len(components),
        "connected_component_sizes": component_sizes,

        # Degree stats
        "degree_min": float(degree_vals.min()) if len(degree_vals) else None,
        "degree_max": float(degree_vals.max()) if len(degree_vals) else None,
        "degree_mean": float(degree_vals.mean()) if len(degree_vals) else None,
        "degree_median": float(np.median(degree_vals)) if len(degree_vals) else None,

        # Edge weights
        **weight_stats,
    }

    # If directed, also include strong/weak component info
    if is_directed:
        stats.update({
            "num_weak_components": len(weak_components),
            "num_strong_components": len(strong_components),
        })

    return stats, component_sizes_df





















def CreatePR_UGraph_OLD(hm_regions_df,
                    hm_pairs_df,
                    add_nonzero_edges=True):
    """
    Create an undirected graph of homologous regions (PRs).
    
    Nodes come from hm_regions_df. Edges are paralogous region-pairs from
    hm_pairs_df. Multiple edges between the same two nodes will have their
    EventsMapped values summed into a single edge weight.

    Args:
        hm_regions_df : DataFrame with:
            Required:
                - 'HmRegion_Num'
                - 'Overlap_Genes'
            Optional:
                - 'N_Events_Mapped'
                - 'N_Events_Putative'
                - 'PR_SetID'
            Missing or NaN event columns are treated as 0.

        hm_pairs_df   : DataFrame with:
            Required:
                - 'Query_HmRegionNum'
                - 'Target_HmRegionNum'
            Optional:
                - 'EventsMapped'
            Missing or NaN EventsMapped is treated as 0.

        add_nonzero_edges : bool
            If True, edges with EventsMapped == 0 are also added.
            If False, edges with EventsMapped == 0 are ignored.

    Returns:
        G : nx.Graph (undirected, no isolates)
    """

    # Build mapping RegionNum → GeneName(s)
    RegionNum_To_Genes = (
        hm_regions_df
        .set_index("HmRegion_Num")["Overlap_Genes"]
        .astype(str)
        .to_dict()
    )

    # New undirected graph
    G = nx.Graph()

    has_pr        = "PR_SetID" in hm_regions_df.columns
    has_mapped    = "N_Events_Mapped" in hm_regions_df.columns
    has_putative  = "N_Events_Putative" in hm_regions_df.columns
    has_edge_evts = "EventsMapped" in hm_pairs_df.columns

    # ------------------
    # Add nodes
    # ------------------
    for _, row in hm_regions_df.iterrows():
        region_num        = row["HmRegion_Num"]
        overlapping_genes = row["Overlap_Genes"]

        # Safely get event counts, default to 0 if missing or NaN
        if has_mapped:
            num_events_mGCE = row["N_Events_Mapped"]
            if pd.isna(num_events_mGCE):
                num_events_mGCE = 0
        else:
            num_events_mGCE = 0

        if has_putative:
            num_events_pGCE = row["N_Events_Putative"]
            if pd.isna(num_events_pGCE):
                num_events_pGCE = 0
        else:
            num_events_pGCE = 0

        node_name = f"{region_num} - {overlapping_genes}"

        node_attrs = dict(
            region_idnum = region_num,
            genes        = overlapping_genes,
            label        = overlapping_genes,
            counts_mGCE  = num_events_mGCE,
            counts_pGCE  = num_events_pGCE,
        )

        if has_pr:
            node_attrs["PR_SetID"] = row["PR_SetID"]

        G.add_node(node_name, **node_attrs)

    # ------------------
    # Add edges (summing weights if edge already exists)
    # ------------------
    for _, row in hm_pairs_df.iterrows():
        q = row["Query_HmRegionNum"]
        t = row["Target_HmRegionNum"]

        # Safely get events, defaulting to 0
        if has_edge_evts:
            events = row["EventsMapped"]
            if pd.isna(events):
                events = 0
        else:
            events = 0

        # Make sure it's numeric/float
        try:
            events = float(events)
        except (ValueError, TypeError):
            events = 0

        qnode = f"{q} - {RegionNum_To_Genes.get(q, 'None')}"
        tnode = f"{t} - {RegionNum_To_Genes.get(t, 'None')}"

        if qnode not in G.nodes or tnode not in G.nodes:
            print(f"WARNING: missing node(s) for edge {qnode} -- {tnode}")
            # continue  # keep behavior as in your original if you prefer

        # Decide if we should add/update an edge for this row
        if (events > 0) or (add_nonzero_edges and events == 0):
            if G.has_edge(qnode, tnode):
                # Sum with existing weight
                prev_weight = G[qnode][tnode].get("weight", 0)
                new_weight  = prev_weight + events
                G[qnode][tnode]["weight"] = new_weight
                G[qnode][tnode]["label"]  = str(new_weight)
            else:
                # Create new edge
                G.add_edge(
                    qnode,
                    tnode,
                    weight = events,
                    label  = str(events),
                )

    # Remove isolated nodes
    G.remove_nodes_from(list(nx.isolates(G)))

    return G





