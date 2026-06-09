# gcutils/mpl_treeviz.py

import pandas as pd
import numpy as np


# from .general import check_overlap_with_gene_group

from itertools import chain

from matplotlib.collections import LineCollection
from matplotlib import markers
from matplotlib.path import Path

import numpy as np

from ete3 import Tree, NodeStyle

from math import floor

from itertools import chain
from matplotlib.collections import LineCollection
from math import floor


# Custom `TreeViz` Matplotlib Functions 

def round_sig(x, sig=2):
    return round(x, sig - int(floor(np.log10(abs(x)))) - 1)


def to_coord(x, y, xmin, xmax, ymin, ymax, plt_xmin, plt_ymin, plt_width, plt_height):
    x = (x - xmin) / (xmax - xmin) * plt_width  + plt_xmin
    y = (y - ymin) / (ymax - ymin) * plt_height + plt_ymin
    return x, y

def mpl_plot_tree(
    tree,
    align_names: bool = False,
    name_offset: float | None = None,
    max_dist: float | None = None,
    font_size: int = 9,
    axe=None,
    show_names: bool = True,
    show_scale_bar: bool = True,
    *,
    # scaling
    x_scale: float = 1.0,      # 0.5 = 50% width (compress), 2.0 = expand
    line_scale: float = 1.0,   # multiply branch line widths
    pad_right: float = 0.02,   # fraction of max_x to pad on the right
    # LEAF markers
    draw_leaf_nodes: bool = True,
    leaf_node_size: float = 8.0,
    leaf_node_shape: str = "o",
    leaf_node_color: str = "black",   # default / fallback color
):
    import matplotlib.pyplot as plt
    from itertools import chain
    from matplotlib.collections import LineCollection
    import numpy as np

    if axe is None:
        axe = plt.subplot(111)

    t = tree.copy()

    # Possibly cut long branches visually
    cut_edge = set()
    if max_dist is not None:
        for n in t.iter_descendants():
            if n.dist > max_dist:
                n.dist -= max_dist
                cut_edge.add(n)

    # Compute max_x AFTER cutting, then scale
    raw_max_x = max((n.get_distance(t) for n in t.iter_leaves()), default=0.0)
    max_x = x_scale * raw_max_x

    # Default label offset (scaled to x_scale so it shrinks with the tree)
    if name_offset is None:
        name_offset = (raw_max_x / 100.0) if raw_max_x > 0 else 0.1
    name_offset_scaled = x_scale * name_offset

    vlinec, vlines = [], []
    hlinec, hlines = [], []
    ali_lines = []
    coords, NameToCoords = {}, {}

    leaves_rev = t.get_leaves()[::-1]
    node_pos = {n2: i for i, n2 in enumerate(leaves_rev)}

    leaf_xs, leaf_ys = [], []
    leaf_colors = []   # <-- NEW: per-leaf colors

    def _draw_edge(child, x_left, cstyle):
        h = node_pos[child]
        x_span = x_scale * child.dist
        x_right = x_left + x_span
        if child in cut_edge:
            offset = (max_x / 600.0) if max_x > 0 else 0.001
            mid = x_left + x_span / 2.0
            hlinec.append(((x_left, h), (mid - offset, h))); hlines.append(cstyle)
            hlinec.append(((mid + offset, h), (x_right, h))); hlines.append(cstyle)
            hlinec.append(((mid, h - 0.05), (mid - 2*offset, h + 0.05))); hlines.append(cstyle)
            hlinec.append(((mid + 2*offset, h - 0.05), (mid, h + 0.05))); hlines.append(cstyle)
            axe.text(mid, h - 0.07, f'+{max_dist:g}', va='top', ha='center', size=2*font_size/3)
        else:
            hlinec.append(((x_left, h), (x_right, h))); hlines.append(cstyle)
        return (x_right, h)

    # Traverse & place
    for n in chain(t.iter_descendants(strategy='postorder'), [t]):
        style = n._get_style()   # ETE NodeStyle dict-like
        x = x_scale * (sum(a.dist for a in n.iter_ancestors()) + n.dist)

        if n.is_leaf():
            y = node_pos[n]
            if draw_leaf_nodes:
                leaf_xs.append(x)
                leaf_ys.append(y)
                # NEW: get per-leaf color from style, fallback to leaf_node_color
                fg = style.get("fgcolor", leaf_node_color)
                leaf_colors.append(fg)

            if show_names:
                if align_names:
                    ali_lines.append(((x, y), (max_x + name_offset_scaled, y)))
                    axe.text(max_x + name_offset_scaled, y, n.name, va='center', size=font_size)
                else:
                    axe.text(x + name_offset_scaled, y, n.name, va='center', size=font_size)
        else:
            y = np.mean([node_pos[ch] for ch in n.children])
            node_pos[n] = y
            vlinec.append(((x, node_pos[n.children[0]]), (x, node_pos[n.children[-1]])))
            vlines.append(style)
            for child in n.children:
                cstyle = child._get_style()
                xy = _draw_edge(child, x, cstyle)
                coords[child] = xy
                if child.name:
                    NameToCoords[child.name] = xy

    # Root left segment (optional/harmless)
    _draw_edge(t, 0.0, t._get_style())

    # Make collections
    lstyles = ['-', '--', ':']
    hline_col = LineCollection(
        hlinec,
        colors=[l['hz_line_color'] for l in hlines],
        linestyle=[lstyles[l['hz_line_type']] for l in hlines],
        linewidth=[line_scale * ((l['hz_line_width'] + 1.0) / 2.0) for l in hlines],
    )
    vline_col = LineCollection(
        vlinec,
        colors=[l['vt_line_color'] for l in vlines],
        linestyle=[lstyles[l['vt_line_type']] for l in vlines],
        linewidth=[line_scale * ((l['vt_line_width'] + 1.0) / 2.0) for l in vlines],
    )
    ali_line_col = LineCollection(ali_lines, colors='k')

    axe.add_collection(hline_col)
    axe.add_collection(vline_col)
    axe.add_collection(ali_line_col)

    # Leaf markers
    if draw_leaf_nodes and leaf_xs:
        size_pts2 = (leaf_node_size ** 2) / 2.0
        axe.scatter(
            leaf_xs,
            leaf_ys,
            s=size_pts2,
            marker=leaf_node_shape,
            c=leaf_colors,        # <-- use per-leaf colors here
            zorder=10,
        )

    # Set xlim to compressed width + some right padding for labels
    pad = pad_right * (max_x if max_x > 0 else 1.0)
    right = max_x + (name_offset_scaled if show_names and align_names else 0.0) + pad
    axe.set_xlim(0.0, right)

    # Scale bar (Optional)
    if show_scale_bar == True:
        xmin, xmax = axe.get_xlim()
        ymin, ymax = axe.get_ylim()
        diffy = max(ymax - ymin, 1.0)
        def _round_sig(x, sig=1):
            from math import floor, log10
            return round(x, sig - int(floor(log10(abs(x)))) - 1) if x > 0 else 1.0
        dist = _round_sig((xmax - xmin) / 5.0, sig=1)
        ybar = ymin - diffy / 100.0
        axe.plot([xmin, xmin + dist], [ybar, ybar], color='k')
        axe.plot([xmin, xmin], [ybar - diffy/200.0, ybar + diffy/200.0], color='k')
        axe.plot([xmin + dist, xmin + dist], [ybar - diffy/200.0, ybar + diffy/200.0], color='k')
        axe.text((2*xmin + dist)/2.0, ybar - diffy/200.0, dist, va='top', ha='center', size=font_size)

    axe.set_axis_off()
    return coords, NameToCoords




    