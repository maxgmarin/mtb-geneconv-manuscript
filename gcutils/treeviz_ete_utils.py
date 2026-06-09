# gcutils/treeviz_utils.py

import pandas as pd
import numpy as np

from ete3 import Tree, NodeStyle

import ete3 as ETE


### Link to ETE3 documentation for tree layout style
# http://etetoolkit.org/docs/latest/tutorial/tutorial_drawing.html#interactive-visualization-of-trees


def style_ETE3_Tree_ByMtbLineage_V1(input_Tree,
                                    in_LinToColor_Dict,
                                    in_Node_To_PrimaryLin_Dict,
                                    verbose = False):

    if in_LinToColor_Dict == None: in_LinToColor_Dict = {}
        
    if in_Node_To_PrimaryLin_Dict == None: in_Node_To_PrimaryLin_Dict = {}
        
    input_Tree_Labeled = input_Tree.copy()
        
    for n in input_Tree_Labeled.traverse():

        Node_PrimaryLin = in_Node_To_PrimaryLin_Dict.get(n.name, "None")
        Node_LinColor = in_LinToColor_Dict.get(Node_PrimaryLin, "black")
        
        nstyle = ETE.NodeStyle() # http://etetoolkit.org/docs/latest/tutorial/tutorial_drawing.html#node-style

        nstyle["size"] = 0
        nstyle["fgcolor"] = "black"

        nstyle["vt_line_color"] = Node_LinColor
        nstyle["hz_line_color"] = Node_LinColor

        nstyle["vt_line_type"] = 0 # 0 solid, 1 dashed, 2 dotted
        nstyle["hz_line_type"] = 0 # 0 solid, 1 dashed, 2 dotted
#        nstyle["hz_line_width"] = 1
#        nstyle["vt_line_width"] = 1
        
        if n.is_leaf(): nstyle["size"] = 1
        
        n.set_style(nstyle)

        if verbose:
            print(n.name, Node_PrimaryLin, Node_LinColor)

    return input_Tree_Labeled







from ete3 import Tree

def style_single_node_in_tree(
    input_tree: Tree,
    target_node_name: str,
    *,
    node_color: str = "red",
    branch_width: int = 1,
    branch_color: str = "red",
    descendant_leaf_size: int | None = None,
) -> Tree:
    """
    Return a copy of the tree with one node visually emphasized, and
    optionally recolor its descendant leaf nodes.

    Styles applied:
      • Branch leading into the target node: width + color
      • Descendant leaf nodes: fgcolor (and optionally size)

    Parameters
    ----------
    input_tree : ete3.Tree
        Original tree. (Not modified.)
    target_node_name : str
        Name of the node to highlight.
    node_color : str, optional
        Color to apply to descendant leaf markers (fgcolor).
    branch_width : int, optional
        Width of the branch leading to the node.
    branch_color : str, optional
        Color of the branch leading to the node.
    descendant_leaf_size : int or None, optional
        If given, update the 'size' of descendant leaf markers to this value.

    Returns
    -------
    Tree
        A *copy* of the input tree with styling applied to the target node
        and its descendant leaves.
    """

    # 1) Copy the full tree
    new_tree = input_tree.copy()

    # 2) Find the target node inside the copy
    matches = new_tree.search_nodes(name=target_node_name)
    if not matches:
        raise ValueError(f"Node '{target_node_name}' not found in the tree.")
    target_node = matches[0]

    # 3) Style the branch leading into the target node (via img_style)
    target_node.img_style.update({
        #"hz_line_width": branch_width,
        #"vt_line_width": branch_width,
        "hz_line_color": branch_color,
        "vt_line_color": branch_color,
    })

    # 4) Update attributes of descendant leaf nodes
    for desc in target_node.traverse():
        if not desc.is_leaf():
            continue

        # Get current style and modify only what we care about
        style = desc._get_style()
        style["fgcolor"] = node_color
        if descendant_leaf_size is not None:
            style["size"] = descendant_leaf_size

        # Re-attach the modified style
        desc.set_style(style)

    return new_tree















## Define functions for subsetting ETE3 trees

def subset_ETE_Tree(i_Target_Name, input_Tree):

    i_SubTree_AtTargetNode = input_Tree.search_nodes(name = i_Target_Name)[0].copy()

    return i_SubTree_AtTargetNode
    


def subset_tree_by_leaf_ids(tree, leaf_ids, keep_only_targets=False):
    """
    Given an ETE3 Tree and a list of leaf IDs (names), return a new Tree
    corresponding to the smallest subtree containing all those leaves.

    Parameters
    ----------
    tree : ete3.Tree
        The original tree.
    leaf_ids : list of str
        List of leaf names (tree leaf .name attributes).
    keep_only_targets : bool, default False
        - If False: return the subtree rooted at the MRCA of leaf_ids
          (may still contain other leaves).
        - If True: return the MRCA subtree but pruned down so that only
          leaves in `leaf_ids` remain.

    Returns
    -------
    ete3.Tree
        A *new* tree object (the original is not modified).
    """
    # Make sure requested leaves exist in the tree
    tree_leaf_names = {leaf.name for leaf in tree.iter_leaves()}
    missing = set(leaf_ids) - tree_leaf_names
    if missing:
        raise ValueError(f"The following leaf IDs are not in the tree: {missing}")

    # Find the MRCA node (this is a Node object within `tree`)
    mrca_node = tree.get_common_ancestor(leaf_ids)

    # Make a deep copy so we don't alter the original tree
    sub_tree = mrca_node.copy()

    if keep_only_targets:
        # Prune to keep only the requested leaf IDs
        # Note: prune() modifies in place, so we're operating on the copy
        sub_tree.prune(leaf_ids, preserve_branch_length=True)

    return sub_tree









def subset_ETE_Tree_Style_GCEChildNode(
    i_ParentNode_Name: str,
    i_ChildNode_Name: str,
    input_Tree: Tree,
    leaf_fgcolor: str = "red",
):
    """
    Return a copy of the subtree rooted at i_ParentNode_Name.

    Effects applied:
      • Branch leading to i_ChildNode_Name → thick dark red
      • Node marker at i_ChildNode_Name → larger dark red sphere
      • Leaf nodes *descending from the child node* → recolored (default: dark red)

    No styling is applied to upstream nodes or sibling clades.
    """

    # 1. Copy subtree at parent node
    sub_tree = input_Tree.search_nodes(name=i_ParentNode_Name)[0].copy()

    # 2. Find the child node inside the copied subtree
    try:
        child_node = sub_tree.search_nodes(name=i_ChildNode_Name)[0]
    # except IndexError:
    #     raise ValueError(
    #         f"Child node '{i_ChildNode_Name}' not found under parent node '{i_ParentNode_Name}'."
    #     )
    
        # 3. Style the branch leading into the child node
        child_node.img_style.update({
            "hz_line_width": 1,
            "vt_line_width": 1,
            "hz_line_color": leaf_fgcolor,
            "vt_line_color": leaf_fgcolor,
        })
    
        # 5. Highlight *descendant leaf nodes* below the child node
        for desc in child_node.traverse():
            if not desc.is_leaf():
                continue
    
            # Get current style and modify only what we care about
            style = desc._get_style()
            style["fgcolor"] = leaf_fgcolor
    
            # Re-attach the modified style
            desc.set_style(style)
    except: 
        print(f"WARNING: Child node '{i_ChildNode_Name}' not found under parent node '{i_ParentNode_Name}'.")
        
    return sub_tree



















#def plotRecombOnPhylo(i_ParentNode_Name, i_ChildNode_Name, input_Tree):
def subset_ETE_Tree_ForRecombViz_OLD(i_ParentNode_Name, i_ChildNode_Name, input_Tree):

    # Draws nodes as small red spheres of diameter equal to 10 pixels
    i_nstyle = ETE.NodeStyle() # http://etetoolkit.org/docs/latest/tutorial/tutorial_drawing.html#node-style
    i_nstyle["shape"] = "sphere"
    i_nstyle["size"] = 10
    i_nstyle["fgcolor"] = "darkred"

    i_SubTree_AtParentNode = input_Tree.search_nodes(name = i_ParentNode_Name)[0].copy()

    for n in i_SubTree_AtParentNode.traverse():
        if n.name == i_ChildNode_Name:
            n.set_style(i_nstyle)
    
    return i_SubTree_AtParentNode






    