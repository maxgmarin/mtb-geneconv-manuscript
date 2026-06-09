# _____.py


import pandas as pd
import numpy as np

from Bio import SeqIO
import bioframe as bf
import seaborn as sns

import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D

from pathlib import Path

from typing import Tuple, Dict, Optional, List

from dna_features_viewer import GraphicFeature, GraphicRecord, CircularGraphicRecord, BiopythonTranslator

from .treeviz_mpl import mpl_plot_tree


NT_COLORS_DEFAULT = {"A":"green","T":"red","C":"blue","G":"orange"}

Query_CoordCols = ("Query_Name", "Query_Start", "Query_End")
HmReg_CoordCols = ("Chr", "Start", "End")
HmRegion_CoordCols = HmReg_CoordCols
Epitope_CoordCols = ("Chrom", "Rv_Start", "Rv_End")
RE_CoordCols = ("seqname", "start_0based", "end_1based")
GenomeAnno_CoordCols = ("Chrom", "Start", "End")






# Function to plot genome/gene annotations on an axis, with options for labeling and styling.



def plot_genome_annotation_on_ax(
    Viz_Start: int,
    Viz_End: int,
    in_Genome_Graphic_Record,
    i_ax,
    *,
    genome_label_threshold = 5,
    top_title = None,
):
    """
    Plot genome annotations for [Viz_Start, Viz_End] onto an existing Axes.

    Parameters
    ----------
    Viz_Start, Viz_End : int
        Genomic window (1-based).
    in_Genome_Graphic_Record :
        dna_features_viewer GraphicRecord object.
    i_ax : matplotlib Axes
        Axis to plot the genome annotations on.
    genome_label_threshold : int
        Label threshold for dna_features_viewer.
    top_title : str, optional
        If provided, place a title on the axis.
    """

    # Crop record
    graphic_cropped = in_Genome_Graphic_Record.crop((Viz_Start, Viz_End + 1))

    # Plot annotations
    if hasattr(graphic_cropped, "plot"):
        graphic_cropped.plot(
            strand_in_label_threshold=genome_label_threshold,
            ax=i_ax,
        )

    # Axis formatting
    i_ax.set_xlim(Viz_Start, Viz_End)
    i_ax.set_xticks([])
    i_ax.set_xlabel("")

    if top_title:
        i_ax.set_title(top_title, fontsize=10)

    return i_ax


###### Functions for viz of variant calls in set of genomes/assemblies #########################

def prepare_PAF_Variants_InWindowForViz(
    df: pd.DataFrame,
    start1: int,
    end1: int,
    *,
    start0_col: str = "Start_0",     # 0-based start (common in your table)
    end1_col: str = "End",           # 1-based inclusive end
    pos1_col: Optional[str] = None,  # alternative to Start_0/End (already 1-based)
    type_col: Optional[str] = "Variant_Type",
    snp_flag_col: Optional[str] = "SNP",
    ref_col: str = "Ref",
    alt_col_candidates: Tuple[str, ...] = ("Alt","ALT","Child_Call","MutAllele","Alt_Base"),
) -> pd.DataFrame:
    """
    Build a windowed variants DF with:
      - pos1 (1-based), end1 (1-based, inclusive)
      - Type in {'SNP','INS','DEL','INDEL','UNKNOWN'}
      - Mut_NT for SNP color (A/T/C/G when available)
    """
    x = df.copy()

    # Compute pos1 / end1
    if pos1_col and pos1_col in x.columns:
        x["pos1"] = x[pos1_col].astype(int)
        x["end1"] = x.get(end1_col, x["pos1"]).astype(int)
    else:
        if start0_col not in x.columns:
            raise ValueError(f"Expected either {pos1_col=} or a 0-based start column '{start0_col}'.")
        x["pos1"] = x[start0_col].astype(int) + 1
        x["end1"] = x.get(end1_col, x["pos1"]).astype(int)

    # Infer Type
    def _infer_type(row):
        vt = str(row.get(type_col, "")).upper() if (type_col and type_col in row) else ""
        if "DEL" in vt: return "DEL"
        if "INS" in vt: return "INS"
        if "SNP" in vt or "SNV" in vt: return "SNP"
        if snp_flag_col and (snp_flag_col in row) and pd.notnull(row[snp_flag_col]):
            return "SNP" if bool(row[snp_flag_col]) else "INDEL"
        ref = str(row.get(ref_col, ""))
        alt = ""
        for c in alt_col_candidates:
            if c in row and pd.notnull(row[c]) and str(row[c]):
                alt = str(row[c])
                break
        if ref and alt:
            if len(ref) == len(alt) == 1: return "SNP"
            if len(alt) > len(ref): return "INS"
            if len(alt) < len(ref): return "DEL"
            return "INDEL"
        return "UNKNOWN"
    x["Type"] = x.apply(_infer_type, axis=1)

    # SNP mutant base for coloring
    def _mut_nt(row):
        for c in alt_col_candidates:
            if c in row and pd.notnull(row[c]) and str(row[c]):
                return str(row[c])[0].upper()
        return ""
    x["Mut_NT"] = x.apply(_mut_nt, axis=1)

    # Window filter (inclusive overlap)
    m = (x["end1"] >= start1) & (x["pos1"] <= end1)
    return x.loc[m].copy()



# ---------- small utilities ----------

def _dedupe_preserve_order(items: List[str]) -> List[str]:
    seen, out = set(), []
    for x in items:
        if x not in seen:
            out.append(x); seen.add(x)
    return out

def _subset_variants(variants_df: pd.DataFrame, isolates: List[str]) -> pd.DataFrame:
    return variants_df[variants_df["SampleID"].isin(isolates)].copy()

def _row_positions(isolates: List[str]) -> tuple[Dict[str, float], np.ndarray]:
    """Top-to-bottom rows with y=0 at bottom (to match your tree alignment)."""
    n = len(isolates)
    y_levels = np.arange(n, dtype=float)[::-1]
    return dict(zip(isolates, y_levels)), y_levels

def _ensure_axes(ax: Optional[plt.Axes], n_rows: int, show_labels: bool) -> tuple[Optional[plt.Figure], plt.Axes]:
    fig = None
    if ax is None:
        base = 0.22 if not show_labels else 0.30
        fig_height = max(3.0 if show_labels else 2.6, base * n_rows + 0.8)
        fig, ax = plt.subplots(figsize=(10, fig_height), dpi=180)
    return fig, ax

# ---------- drawing primitives ----------

def _draw_row_backgrounds(ax: plt.Axes,
                          isolates: List[str],
                          y_map: Dict[str, float],
                          x0: int, x1: int,
                          row_h: float,
                          bg_color: str, bg_alpha: float) -> None:
    width = x1 - x0
    for sid in isolates:
        y = y_map[sid]
        ax.add_patch(Rectangle((x0, y - row_h/2), width, row_h,
                               facecolor=bg_color, edgecolor="none",
                               alpha=bg_alpha, zorder=0))

def _draw_snps(ax: plt.Axes,
               sdf: pd.DataFrame,
               ymin: float, ymax: float,
               x0: int, x1: int,
               nt_colors: Dict[str, str],
               snp_line_width: float) -> None:
    snps = sdf[sdf["Type"] == "SNP"]
    if snps.empty:
        return
    # Group by first NT of alt for color
    for nt, g in snps.groupby(snps["Mut_NT"].str.upper().str[0].fillna("")):
        color = nt_colors.get(nt, "black")
        for xv in g["pos1"].astype(float).values:
            if x0 <= xv <= x1:
                ax.vlines(x=xv, ymin=ymin, ymax=ymax,
                          lw=snp_line_width, color=color, alpha=0.95, zorder=3)

def _draw_insertions(ax: plt.Axes,
                     sdf: pd.DataFrame,
                     ymin: float, ymax: float,
                     x0: int, x1: int,
                     ins_color: str, ins_line_width: float) -> None:
    ins = sdf[sdf["Type"] == "INS"]
    if ins.empty:
        return
    for xv in ins["pos1"].astype(float).values:
        if x0 <= xv <= x1:
            ax.vlines(x=xv, ymin=ymin, ymax=ymax,
                      lw=ins_line_width, color=ins_color, alpha=0.95, zorder=4)

def _draw_deletions(ax: plt.Axes,
                    sdf: pd.DataFrame,
                    y: float,
                    row_h: float,
                    x0: int, x1: int,
                    del_color: str, del_hline_width: float, del_vcap_width: float, del_alpha: float) -> None:
    dels = sdf[sdf["Type"] == "DEL"]
    if dels.empty:
        return
    ymin = y - row_h/2
    ymax = y + row_h/2
    for _, r in dels.iterrows():
        a = max(x0, float(r["pos1"]))
        b = min(x1, float(r["end1"]))
        if b < a: a, b = b, a
        if np.isclose(a, b):
            ax.vlines(a, ymin=ymin, ymax=ymax, color=del_color, lw=del_vcap_width, alpha=del_alpha, zorder=4)
        else:
            ax.hlines(y=y, xmin=a, xmax=b, color=del_color, lw=del_hline_width, alpha=del_alpha, zorder=4)
            ax.vlines([a, b], ymin=ymin, ymax=ymax, color=del_color, lw=del_vcap_width, alpha=del_alpha, zorder=4)

def _format_x_axis(ax: plt.Axes, start: int, end: int, nbins: int = 6) -> None:
    ax.set_xlim(start, end)
    fmt = mticker.ScalarFormatter(useOffset=False)
    fmt.set_scientific(False)
    ax.xaxis.set_major_formatter(fmt)
    ax.ticklabel_format(style='plain', axis='x', useOffset=False)
    ax.xaxis.set_major_locator(mticker.MaxNLocator(nbins=nbins, integer=True))



# ---------- main wrapper ----------
def plot_variants_for_isolates_simple(
    variants_df: pd.DataFrame,
    isolates: List[str],
    region_start1: int,
    region_end1: int,
    *,
    ax: Optional[plt.Axes] = None,
    show_labels: bool = True,
    title: Optional[str] = None,
    nt_colors: Optional[Dict[str, str]] = None,
    bar_height: float = 0.34,
    snp_tick_height_factor: float = 1.0,
    snp_line_width: float = 0.9,
    bg_color: str = "#CDCDCD",
    bg_alpha: float = 0.8,
    ins_color: str = "#6a0dad",
    ins_line_width: float = 2.0,
    del_color: str = "black",
    del_hline_width: float = 1.2,
    del_vcap_width: float = 1.2,
    del_alpha: float = 1.0,
    return_summary: bool = True,
):
    """
    Plot per-isolate SNPs/INS/DEL across a genomic window.

    Input variants_df is automatically normalized + filtered using
    prepare_PAF_Variants_InWindowForViz().
    """

    nt_colors = nt_colors or NT_COLORS_DEFAULT

    # --- 0) Normalize + window variants ---
    window_df = prepare_PAF_Variants_InWindowForViz(
        variants_df,
        start1=region_start1,
        end1=region_end1,
    )

    # --- 1) Order isolates & subset ---
    ordered = _dedupe_preserve_order(isolates)
    sub = _subset_variants(window_df, ordered)

    # --- 2) Axes setup ---
    fig, ax = _ensure_axes(ax, n_rows=len(ordered), show_labels=show_labels)
    y_map, y_levels = _row_positions(ordered)
    row_h = bar_height * snp_tick_height_factor

    # --- 3) Background rows ---
    _draw_row_backgrounds(ax, ordered, y_map,
                          region_start1, region_end1,
                          row_h, bg_color, bg_alpha)

    # --- 4) Draw variants ---
    had_variants = set()
    for sid in ordered:
        y = y_map[sid]
        sdf = sub[sub["SampleID"] == sid]
        if sdf.empty:
            continue
        had_variants.add(sid)

        ymin, ymax = y - row_h/2, y + row_h/2
        _draw_snps(ax, sdf, ymin, ymax, region_start1, region_end1,
                   nt_colors, snp_line_width)
        _draw_insertions(ax, sdf, ymin, ymax, region_start1, region_end1,
                         ins_color, ins_line_width)
        _draw_deletions(ax, sdf, y, row_h, region_start1, region_end1,
                        del_color, del_hline_width, del_vcap_width, del_alpha)

    # --- 5) Cosmetics ---
    _format_x_axis(ax, region_start1, region_end1, nbins=6)

    ax.set_ylim(-0.8, len(ordered) - 1 + 0.8)
    if show_labels:
        ax.set_yticks(y_levels)
        ax.set_yticklabels(ordered, fontsize=8)
    else:
        ax.set_yticks([])

    ax.set_xlabel("Genomic position (1-based)")
    ax.grid(axis="x", linestyle=":", linewidth=0.5, alpha=0.35)

    if title:
        ax.set_title(title, fontsize=10)

    for spine in ("left", "top", "right"):
        ax.spines[spine].set_visible(False)

    # --- 6) Summary ---
    summary = None
    if return_summary:
        missing = [sid for sid in ordered if sid not in had_variants]
        summary = {
            "window": (region_start1, region_end1),
            "isolates_requested": ordered,
            "isolates_with_variants": sorted(list(had_variants), key=ordered.index),
            "isolates_without_variants": missing,
            "n_rows": len(ordered),
        }

    return (fig, ax, summary) #if fig is not None else (None, ax, summary)





def plot_paralog_variants_on_ax(
    ax,
    in_df_all,
    region_start: int,
    region_end: int,
    *,
    group_col: str = "QueryParalog_RegionID",
    pos1_col: str = "Target_Start",
    end1_col: str = "Target_End",
    type_col: str = "Type",
    alt_col: str = "Alt",
    nt_colors=None,
    bar_height: float = 0.36,
    snp_line_width: float = 0.9,
    ins_color: str = "#6a0dad",
    ins_line_width: float = 1.8,
    del_color: str = "black",
    del_hline_width: float = 1.2,
    del_vcap_width: float = 1.2,
    bg_color: str = "#CDCDCD",
    bg_alpha: float = 0.8,
    show_labels: bool = True,
    x_integer_ticks_nbins: int = 6,
    title: str | None = None,
    hide_xaxis: bool = False,
):
    """
    Minimal visualization:
    - one row per QueryParalog_RegionID
    - SNPs colored by Alt base
    - INS = purple vertical ticks
    - DEL = black horizontal bars with caps
    """
    nt_colors = nt_colors or NT_COLORS_DEFAULT
    
    print(in_df_all.shape, in_df_all[group_col].nunique())

    
    target_region = f"NC_000962.3:{region_start}-{region_end}" 
        
    df = bf.select(in_df_all, target_region,
                        cols=["Target_Name", "Target_Start", "Target_End"])

    print(df.shape, df[group_col].nunique())

    labels = list(dict.fromkeys(df[group_col].astype(str)))  # preserve order
    y_levels = np.arange(len(labels))[::-1]
    y_map = dict(zip(labels, y_levels))

    # background bars
    for lab in labels:
        y = y_map[lab]

        row = df[df[group_col] == lab]
        
        paralog_Query_Aln_Start = row["Aln_Target_Start"].values[0]
        paralog_Query_Aln_End   = row["Aln_Target_End"].values[0]

        win_w = paralog_Query_Aln_End - paralog_Query_Aln_Start

        print(lab, "---", paralog_Query_Aln_Start, paralog_Query_Aln_End, win_w)

        ax.add_patch(Rectangle(
            (paralog_Query_Aln_Start, y - bar_height/2), win_w, bar_height,
            facecolor=bg_color, edgecolor="none", alpha=bg_alpha, zorder=0
        ))

    # draw variants
    for lab in labels:
        row = df[df[group_col] == lab]
        y = y_map[lab]
        ymin, ymax = y - bar_height/2, y + bar_height/2

        for _, r in row.iterrows():
            t = r[type_col]
            a, b = float(r[pos1_col]), float(r[end1_col])
            alt = str(r[alt_col]).upper()[0] if isinstance(r[alt_col], str) else ""

            if t == "SNP":
                if region_start <= a <= region_end:
                    ax.vlines(a, ymin=ymin, ymax=ymax,
                              lw=snp_line_width,
                              color=nt_colors.get(alt, "black"), zorder=3)
            elif t == "INS":
                if region_start <= a <= region_end:
                    ax.vlines(a, ymin=ymin, ymax=ymax,
                              lw=ins_line_width, color=ins_color, zorder=4)
            elif t == "DEL":
                if b < a: a, b = b, a
                a = max(region_start, a)
                b = min(region_end, b)
                if np.isclose(a, b):
                    ax.vlines(a, ymin=ymin, ymax=ymax,
                              color=del_color, lw=del_vcap_width, zorder=4)
                else:
                    ax.hlines(y, xmin=a, xmax=b,
                              color=del_color, lw=del_hline_width, zorder=4)
                    ax.vlines([a, b], ymin=ymin, ymax=ymax,
                              color=del_color, lw=del_vcap_width, zorder=4)

    # cosmetics
    ax.set_xlim(region_start, region_end)
    ax.set_ylim(-0.8, len(labels)-1+0.8)
    if show_labels:
        ax.set_yticks(y_levels)
        ax.set_yticklabels(labels, fontsize=8)
    else:
        ax.set_yticks([])

    ax.xaxis.set_major_formatter(mticker.ScalarFormatter(useOffset=False))
    ax.ticklabel_format(style='plain', axis='x', useOffset=False)
    ax.xaxis.set_major_locator(mticker.MaxNLocator(nbins=x_integer_ticks_nbins, integer=True))
    ax.set_xlabel("Genomic position (Query coordinates)")
    for spine in ("left","top","right"):
        ax.spines[spine].set_visible(False)
    ax.grid(axis="x", linestyle=":", linewidth=0.5, alpha=0.35)

    if hide_xaxis:
        ax.set_xticks([])
        ax.set_xlabel("")
    
    if title:
        ax.set_title(title, fontsize=10)

    return ax























def plot_Anno_with_empty_bottom(
    Viz_Start: int,
    Viz_End: int,
    in_Genome_Graphic_Record,   # dna_features_viewer GraphicRecord (H37Rv)
    *,
    figsize: Tuple[float, float] = (9, 3.6),
    dpi: int = 180,
    height_ratios: Tuple[float, float] = (1.0, 4.0),
    genome_label_threshold: int = 5,
    top_title: Optional[str] = None,
    show_bottom_xlabel: bool = True,
    x_integer_ticks_nbins: int = 6,
    tight_layout: bool = True,
) -> Tuple[plt.Figure, Dict[str, plt.Axes]]:
    """
    Make a two-panel canvas:
      - Top: H37Rv genome annotations over [Viz_Start, Viz_End]
      - Bottom: empty axis (pre-formatted) to overlay mutations/variants later

    Returns
    -------
    fig, {"Genome_Anno_ax": ax0, "Bottom_ax": ax1}
    """
    # Crop to window
    graphic_cropped = in_Genome_Graphic_Record.crop((Viz_Start, Viz_End + 1))

    # Figure & axes
    fig, (ax_top, ax_bottom) = plt.subplots(
        2, 1, figsize=figsize, dpi=dpi,
        gridspec_kw={"height_ratios": height_ratios},
        sharex=False
    )

    # --- Top: genome annotations ---
    if hasattr(graphic_cropped, "plot"):
        graphic_cropped.plot(
            strand_in_label_threshold=genome_label_threshold,
            ax=ax_top,
        )
    ax_top.set_xlim(Viz_Start, Viz_End)
    ax_top.set_xticks([])
    ax_top.set_xlabel("")
    if top_title:
        ax_top.set_title(top_title, fontsize=10)

    # --- Bottom: leave empty but pre-format x-axis for the same window ---
    ax_bottom.set_xlim(Viz_Start, Viz_End)
    if show_bottom_xlabel:
        ax_bottom.set_xlabel("Genomic position (1-based)")
    # integer ticks, no scientific notation
    fmt = mticker.ScalarFormatter(useOffset=False)
    fmt.set_scientific(False)
    ax_bottom.xaxis.set_major_formatter(fmt)
    ax_bottom.ticklabel_format(style="plain", axis="x", useOffset=False)
    ax_bottom.xaxis.set_major_locator(mticker.MaxNLocator(nbins=x_integer_ticks_nbins, integer=True))
    # keep the bottom panel minimal
    for spine in ("left", "top", "right"):
        ax_bottom.spines[spine].set_visible(False)
    ax_bottom.grid(axis="x", linestyle=":", linewidth=0.5, alpha=0.35)

    if tight_layout:
        fig.tight_layout()

    return fig, {"Genome_Anno_ax": ax_top, "Bottom_ax": ax_bottom}




def plot_Anno_WiVariantCalls_2Panel(
    Viz_Start: int,
    Viz_End: int,
    in_Genome_Graphic_Record,
    variants_df: pd.DataFrame,
    isolates: List[str],
    *,
    # canvas options
    figsize: Tuple[float, float] = (7, 9),
    height_ratios: Tuple[float, float] = (1, 25),
    genome_label_threshold: int = 5,
    top_title: Optional[str] = None,
    show_bottom_xlabel: bool = True,
    x_integer_ticks_nbins: int = 6,
    tight_layout: bool = True,
    # variant plotting options (forwarded)
    show_labels: bool = False,
    bar_height: float = 0.34,
    snp_tick_height_factor: float = 2.0,
    snp_line_width: float = 1.0,
    nt_colors: Optional[Dict[str, str]] = None,
    bg_color: str = "#CDCDCD",
    bg_alpha: float = 0.8,
    ins_color: str = "#6a0dad",
    ins_line_width: float = 2.0,
    del_color: str = "black",
    del_hline_width: float = 1.2,
    del_vcap_width: float = 1.2,
    del_alpha: float = 1.0,
):
    """
    Make a 2-panel figure:
      Top  : genome annotations over [Viz_Start, Viz_End]
      Bottom: per-isolate variant track over the same window

    Returns
    -------
    fig, axes_dict, var_in_win, summary
      axes_dict = {"Genome_Anno_ax": ax_top, "Bottom_ax": ax_bottom}
      summary   = return from plot_variants_for_isolates_simple
    """
    # 1) Canvas with top annotations + empty bottom
    fig, axes = plot_Anno_with_empty_bottom(
        Viz_Start=Viz_Start,
        Viz_End=Viz_End,
        in_Genome_Graphic_Record=in_Genome_Graphic_Record,
        figsize=figsize,
        height_ratios=height_ratios,
        genome_label_threshold=genome_label_threshold,
        top_title=top_title,
        show_bottom_xlabel=show_bottom_xlabel,
        x_integer_ticks_nbins=x_integer_ticks_nbins,
        tight_layout=False,  # we'll tight_layout at the very end
    )

    # 2) Window the variants DF to the same region
    var_in_win = prepare_PAF_Variants_InWindowForViz(
        variants_df,
        start1=Viz_Start,
        end1=Viz_End
    )

    # 3) Plot variants on the bottom axis
    _, _, summary = plot_variants_for_isolates_simple(
        variants_df=var_in_win,
        isolates=isolates,
        region_start1=Viz_Start,
        region_end1=Viz_End,
        ax=axes["Bottom_ax"],
        show_labels=show_labels,
        nt_colors=nt_colors,
        bar_height=bar_height,
        snp_tick_height_factor=snp_tick_height_factor,
        snp_line_width=snp_line_width,
        bg_color=bg_color,
        bg_alpha=bg_alpha,
        ins_color=ins_color,
        ins_line_width=ins_line_width,
        del_color=del_color,
        del_hline_width=del_hline_width,
        del_vcap_width=del_vcap_width,
        del_alpha=del_alpha,
        return_summary=True,
    )

    if tight_layout:
        fig.tight_layout()

    return fig, axes, var_in_win, summary








def setup_3panelfig_phylo_variants_anno_axes(
    Viz_Start: int,
    Viz_End: int,
    *,
    figsize: Tuple[float, float] = (16, 4),
    dpi: int = 180,
    # [left, bottom, width, height] in figure coordinates
    tree_rect = (0.1, 0.1, 0.3, 0.8),
    variants_rect = (0.45, 0.1, 0.7, 0.8),
    anno_rect = (0.45, 0.9, 0.7, 0.1),
) -> Tuple[plt.Figure, Dict[str, plt.Axes]]:
    """
    Set up a figure with three axes:

      - ax1_L_Phylo   : left panel for phylogenetic tree
      - ax2_BR_Variants: right panel for per-isolate variants (shares y with tree)
      - ax3_UR_Anno   : upper-right panel for genome annotations

    Parameters
    ----------
    Viz_Start, Viz_End : int
        Genomic window for the right-side panels' x-limits.
    figsize : (float, float)
        Figure size passed to plt.figure.
    dpi : int
        Figure DPI.
    tree_rect, variants_rect, anno_rect : tuple
        [left, bottom, width, height] in figure coordinates for each axis.

    Returns
    -------
    fig : Figure
    axes : dict
        {
          "ax1_L_Phylo": ax1_L_Phylo,
          "ax2_R_Variants": ax2_R_Variants,
          "ax3_UR_Anno": ax3_UR_Anno,
        }
    """
    fig = plt.figure(figsize=figsize, dpi=dpi)

    # Left: phylogeny axis
    ax1_L_Phylo = fig.add_axes(tree_rect)

    # Right: variants axis (shares y with tree)
    ax2_R_Variants = fig.add_axes(variants_rect, sharey=ax1_L_Phylo)
    ax2_R_Variants.set_xlim(Viz_Start, Viz_End)

    # Upper-right: annotation axis (independent x, but matched limits)
    ax3_UR_Anno = fig.add_axes(anno_rect)
    ax3_UR_Anno.set_xlim(Viz_Start, Viz_End)

    axes = {
        "ax1_L_Phylo": ax1_L_Phylo,
        "ax2_BR_Variants": ax2_R_Variants,
        "ax3_UR_Anno": ax3_UR_Anno,
    }

    return fig, axes


def setup_4panelfig_phylo_variants_anno_paralogvar_axes(
    Viz_Start: int,
    Viz_End: int,
    *,
    figsize: Tuple[float, float] = (16, 5),
    dpi: int = 180,
    # [left, bottom, width, height] in figure coordinates
    tree_rect = (0.1, 0.1, 0.3, 0.8),
    variants_rect = (0.45, 0.1, 0.7, 0.8),
    anno_rect = (0.45, 0.9, 0.7, 0.1),
    paralogvariants_rect = (0.45, 1.1, 0.7, 0.3),
) -> Tuple[plt.Figure, Dict[str, plt.Axes]]:
    """
    Set up a figure with three axes:

      - ax1_L_Phylo   : left panel for phylogenetic tree
      - ax2_BR_Variants: right panel for per-isolate variants (shares y with tree)
      - ax3_UR_Anno   : upper-right panel for genome annotations
      - ax3_UR_Anno   : upper-right panel for genome annotations

    Parameters
    ----------
    Viz_Start, Viz_End : int
        Genomic window for the right-side panels' x-limits.
    figsize : (float, float)
        Figure size passed to plt.figure.
    dpi : int
        Figure DPI.
    tree_rect, variants_rect, anno_rect : tuple
        [left, bottom, width, height] in figure coordinates for each axis.

    Returns
    -------
    fig : Figure
    axes : dict
        {
        "ax1_L_Phylo": ax1_L_Phylo,
        "ax2_BR_Variants": ax2_R_Variants,
        "ax3_UR_Anno": ax3_UR_Anno,
        "ax4_UUR_ParalogVariants": ax3_UR_Anno,
        }
    """
    fig = plt.figure(figsize=figsize, dpi=dpi)

    # Left: phylogeny axis
    ax1_L_Phylo = fig.add_axes(tree_rect)

    # Right: variants axis (shares y with tree)
    ax2_R_Variants = fig.add_axes(variants_rect, sharey=ax1_L_Phylo)
    ax2_R_Variants.set_xlim(Viz_Start, Viz_End)

    # Upper-right: annotation axis (independent x, but matched limits)
    ax3_UR_Anno = fig.add_axes(anno_rect)
    ax3_UR_Anno.set_xlim(Viz_Start, Viz_End)


    # Upper-Upper-right: paralog variants axis (independent x, but matched limits)
    ax4_UUR_ParalogVariants = fig.add_axes(paralogvariants_rect)
    ax4_UUR_ParalogVariants.set_xlim(Viz_Start, Viz_End)

    axes = {
        "ax1_L_Phylo": ax1_L_Phylo,
        "ax2_BR_Variants": ax2_R_Variants,
        "ax3_UR_Anno": ax3_UR_Anno,
        "ax4_UUR_ParalogVariants": ax4_UUR_ParalogVariants,
    }

    return fig, axes










def viz_region_Tree_Variants_Anno_ParalogVariants(
    Viz_Start: int,
    Viz_End: int,
    *,
    in_Genome_Graphic_Record,
    inputTree,
    i_AllVar_DF: pd.DataFrame,
    i_ParalogVar_DF: pd.DataFrame,
    i_GCE_DF: bool = None,
    # layout
    figsize: Tuple[float, float] = (16, 7),
    dpi: int = 300,
    tree_rect = (0.28, 0.1, 0.3, 0.8),
    variants_rect = (0.45, 0.1, 0.7, 0.8),
    anno_rect = (0.45, 0.9, 0.7, 0.1),
    paralogvariants_rect = (0.45, 1.0, 0.7, 0.2),
    # options
    show_tree_names: bool = False,
    show_variant_labels: bool = False,
    paralog_title = None,
    snp_line_width: float = 0.7,
    paralog_snp_line_width: float = 0.4,
    tree_pad_right: float = 1.0,
    tree_leaf_node_size: float = 3.0,
    # recombination overlay style
    recomb_facecolor: str = "red",
    recomb_alpha: float = 0.2,
    recomb_linewidth: float = 1.0,
    recomb_row_pad: float = 0.05,
    recomb_show_event_id: bool = True,
    use_one_based_axis: bool = True,
    return_summaries: bool = True,
):

    fig, axes = setup_4panelfig_phylo_variants_anno_paralogvar_axes(
        Viz_Start,
        Viz_End,
        figsize=figsize,
        dpi=dpi,
        tree_rect=tree_rect,
        variants_rect=variants_rect,
        anno_rect=anno_rect,
        paralogvariants_rect=paralogvariants_rect,
    )

    ax_tree = axes["ax1_L_Phylo"]
    ax_vars = axes["ax2_BR_Variants"]
    ax_anno = axes["ax3_UR_Anno"]
    ax_paralog = axes["ax4_UUR_ParalogVariants"]

    # 1) Infer isolate order from tree
    isolates = [lf.name for lf in inputTree.get_leaves()]

    # 2) Plot phylogeny
    coords, NameToCoords = mpl_plot_tree(
        inputTree,
        name_offset=1.2,
        pad_right=tree_pad_right,
        axe=ax_tree,
        show_names=show_tree_names,
        draw_leaf_nodes = True,
        leaf_node_size = tree_leaf_node_size,
        show_scale_bar=False,
    )

    # 3) Plot per-isolate variants
    _, _, variants_summary = plot_variants_for_isolates_simple(
        i_AllVar_DF,
        isolates,
        region_start1=Viz_Start,
        region_end1=Viz_End,
        ax=ax_vars,
        show_labels=show_variant_labels,
        bar_height=0.7,
        snp_tick_height_factor=1.0,
        snp_line_width=snp_line_width,
    )

    # 4) OPTIONAL: overlay recombination tracts
    if i_GCE_DF is not None:
        patches, gce_summary = overlay_gc_event_highlights_on_variants(
            ax=ax_vars,
            events_df=i_GCE_DF,
            input_tree=inputTree,
            isolates_order=isolates,
            region_start1=Viz_Start,
            region_end1=Viz_End,
            facecolor=recomb_facecolor,
            edgecolor="black",
            alpha=recomb_alpha,
            linewidth=recomb_linewidth,
            row_pad = recomb_row_pad,
            show_event_id = recomb_show_event_id,
            use_one_based_axis=use_one_based_axis,
            zorder=2,
            return_summary=True,
        )
    else:
        gce_summary = None

    # 5) Genome annotations
    plot_genome_annotation_on_ax(
        Viz_Start,
        Viz_End,
        in_Genome_Graphic_Record,
        i_ax=ax_anno,
    )

    # 6) Paralog variants
    plot_paralog_variants_on_ax(
        ax_paralog,
        i_ParalogVar_DF,
        region_start=Viz_Start,
        region_end=Viz_End,
        snp_line_width=paralog_snp_line_width,
        title=paralog_title,
        hide_xaxis=True,
    )

    # 7) Summaries (optional)
    summaries = None
    if return_summaries:
        summaries = {
            "variants_summary": variants_summary,
            "gce_summary": gce_summary,
        }

    return fig, axes, summaries




def viz_region_Tree_Variants_Anno(
    Viz_Start: int,
    Viz_End: int,
    *,
    in_Genome_Graphic_Record,
    inputTree,
    i_AllVar_DF: pd.DataFrame,
    i_GCE_DF: pd.DataFrame = None,
    # layout
    figsize: Tuple[float, float] = (16, 5),
    dpi: int = 300,
    tree_rect = (0.28, 0.1, 0.3, 0.8),
    variants_rect = (0.45, 0.1, 0.7, 0.8),
    anno_rect = (0.45, 0.9, 0.7, 0.1),
    # options
    show_tree_names: bool = False,
    show_variant_labels: bool = False,
    snp_line_width: float = 0.7,
    tree_pad_right: float = 1.0,
    tree_leaf_node_size: float = 3.0,
    # recombination overlay style
    recomb_facecolor: str = "red",
    recomb_alpha: float = 0.2,
    recomb_linewidth: float = 1.0,
    recomb_row_pad: float = 0.05,
    recomb_show_event_id: bool = True,
    use_one_based_axis: bool = True,
    return_summaries: bool = True,
):

    # --- 0) Setup figure + axes
    fig, axes = setup_3panelfig_phylo_variants_anno_axes(
        Viz_Start,
        Viz_End,
        figsize=figsize,
        dpi=dpi,
        tree_rect=tree_rect,
        variants_rect=variants_rect,
        anno_rect=anno_rect,
    )

    ax_tree = axes["ax1_L_Phylo"]
    ax_vars = axes["ax2_BR_Variants"]
    ax_anno = axes["ax3_UR_Anno"]

    # --- 1) Infer isolate order from input tree
    isolates = [lf.name for lf in inputTree.get_leaves()]

    # --- 2) Plot phylogeny
    coords, NameToCoords = mpl_plot_tree(
        inputTree,
        name_offset=1.2,
        pad_right=tree_pad_right,
        axe=ax_tree,
        show_names=show_tree_names,
        draw_leaf_nodes=True,
        leaf_node_size = tree_leaf_node_size,
        show_scale_bar=False,
    )

    # --- 3) Plot variants (per isolate)
    _, _, variants_summary = plot_variants_for_isolates_simple(
        i_AllVar_DF,
        isolates,
        region_start1=Viz_Start,
        region_end1=Viz_End,
        ax=ax_vars,
        show_labels=show_variant_labels,
        bar_height=0.7,
        snp_tick_height_factor=1.0,
        snp_line_width=snp_line_width,
    )

    # --- 4) OPTIONAL: overlay recombination tracts
    if i_GCE_DF is not None:
        patches, gce_summary = overlay_gc_event_highlights_on_variants(
            ax=ax_vars,
            events_df=i_GCE_DF,
            input_tree=inputTree,
            isolates_order=isolates,
            region_start1=Viz_Start,
            region_end1=Viz_End,
            facecolor=recomb_facecolor,
            edgecolor="black",
            alpha=recomb_alpha,
            linewidth=recomb_linewidth,
            row_pad=recomb_row_pad,
            show_event_id = recomb_show_event_id,
            use_one_based_axis=use_one_based_axis,
            zorder=2,
            return_summary=True,
        )
    else:
        gce_summary = None

    # --- 5) Genome annotations (upper-right panel)
    plot_genome_annotation_on_ax(
        Viz_Start,
        Viz_End,
        in_Genome_Graphic_Record,
        i_ax=ax_anno,
    )

    # --- 6) Summaries
    summaries = None
    if return_summaries:
        summaries = {
            "variants_summary": variants_summary,
            "gce_summary": gce_summary,
        }

    return fig, axes, summaries







#################################################################################################




######### Functions for overlaying Gubbins recombination event blocks on top of variant visualizations #########

import numpy as np
import pandas as pd
from matplotlib.patches import Rectangle
from typing import Iterable, Dict, List, Tuple, Optional

# ---------- helpers ----------

def _descendant_leaves_names(tree, node_name: str) -> List[str]:
    """Return leaf names under the node named `node_name` (including the node if it's a leaf)."""
    hits = tree.search_nodes(name=node_name)
    if not hits:
        return []
    node = hits[0]
    return [n.name for n in node.iter_leaves()]

def _isolate_y_map(isolates: List[str]) -> Dict[str, float]:
    """
    Map isolate -> y row (top to bottom), matching plot_variants_for_isolates_simple().
    That function uses np.arange(n)[::-1], so we replicate that here.
    """
    n = len(isolates)
    return dict(zip(isolates, np.arange(n, dtype=float)[::-1]))

def _event_x_bounds(start_0based: int, end_1based: int, *, use_one_based: bool = True) -> Tuple[float, float]:
    """Convert event coords to plotting coords. If 1-based axis: start+1, end stays inclusive."""
    x0 = (int(start_0based) + 1) if use_one_based else int(start_0based)
    x1 = int(end_1based)
    if x1 < x0:
        x0, x1 = x1, x0
    return float(x0), float(x1)

def _clip_interval(x0: float, x1: float, win0: float, win1: float) -> Optional[Tuple[float, float]]:
    """Clip [x0,x1] to [win0,win1]; return None if no overlap."""
    a = max(x0, win0)
    b = min(x1, win1)
    if b < a:
        return None
    return (a, b)

def _descendant_row_span(desc_leaves: Iterable[str], y_map: Dict[str, float], row_pad: float = 0.0) -> Optional[Tuple[float, float]]:
    """
    From descendant leaf names, compute vertical extent in row units.
    We expand by ±row_pad to give a little breathing room if desired.
    """
    ys = [y_map[s] for s in desc_leaves if s in y_map]
    if not ys:
        return None
    y_min = min(ys) - 0.5 + row_pad
    y_max = max(ys) + 0.5 - row_pad
    # Ensure non-negative height
    return (y_min, max(y_max, y_min + 1e-9))

# ---------- main overlay ----------

def overlay_gc_event_highlights_on_variants(
    ax,
    events_df: pd.DataFrame,
    input_tree,
    isolates_order: List[str],
    region_start1: int,
    region_end1: int,
    *,
    # column names in events_df
    start_col: str = "start_0based",
    end_col: str   = "end_1based",
    node_col: str  = "Child_Node",
    show_event_id: bool = True,
    # rendering
    facecolor: str = "red",
    edgecolor: str = "black", #Optional[str] = None,
    alpha: float = 0.20,
    linewidth: float = 1.0,
    row_pad: float = 0.05,          # shrink vertical box a touch so it sits inside row bands
    use_one_based_axis: bool = True, # your variants axis is 1-based; keep True
    zorder: int = 2,
    return_summary: bool = True,
):
    """
    Draw a transparent rectangle for each GC event:
      x-span = event genomic range (clipped to [region_start1, region_end1])
      y-span = rows corresponding to descendant leaves of the event's Child_Node

    Parameters
    ----------
    ax : matplotlib Axes where variants are drawn.
    events_df : DataFrame with at least [start_col, end_col, node_col].
    input_tree : ete3.Tree
    isolates_order : list of isolate IDs, in the same order passed to plot_variants_for_isolates_simple.
    region_start1, region_end1 : viewing window (1-based, inclusive).

    Returns
    -------
    patches : list[Rectangle]
    summary_df : pd.DataFrame (optional; if return_summary=True)
    """
    edgecolor = edgecolor if edgecolor is not None else facecolor
    
    y_map = _isolate_y_map(isolates_order)
    patches = []
    summaries = []

    # Step 1: Subset for only events that overlap the region
    region_start0 = region_start1 - 1
    Region_Coords = f"NC_000962.3:{region_start0}-{region_end1}"

    i_Ovrlap_GCE_DF = bf.select(events_df, Region_Coords, cols = ("seqname", "start_1based", "end_1based"))


    for i, row in i_Ovrlap_GCE_DF.iterrows():
        # 1) x-range (axis units)
        try:
            x0_raw, x1_raw = _event_x_bounds(row[start_col], row[end_col], use_one_based=use_one_based_axis)
        except Exception:
            continue

        clipped = _clip_interval(x0_raw, x1_raw, float(region_start1), float(region_end1))
        if clipped is None:
            continue
        x0, x1 = clipped
        width = max(0.0, x1 - x0)

        # 2) descendant rows
        node_name = str(row.get(node_col, ""))
        desc = _descendant_leaves_names(input_tree, node_name)
        span = _descendant_row_span(desc, y_map, row_pad=row_pad)
        if span is None:
            # no descendants in panel
            continue
        y0, y1 = span
        height = y1 - y0

        rect = Rectangle(
            (x0, y0),
            width,
            height,
            facecolor=facecolor,
            edgecolor=edgecolor,
            alpha=alpha,
            linewidth=linewidth,
            zorder=zorder
        )
        ax.add_patch(rect)
        patches.append(rect)


        # draw EventID label to the left of the rectangle
        if show_event_id:
            y_mid = y0 + height / 2.0
            ax.text(
                x0 - 5.0,          # fixed left padding
                y_mid,
                str(row["EventID"]),
                ha="right",
                va="center",
                fontsize=10,
                color= "black", #edgecolor,
                alpha=0.8,
                zorder=zorder + 1,
                clip_on=False,
            )
            
        
        if return_summary:
            summaries.append({
                "Event_index": i,
                "Child_Node": node_name,
                "x0": x0, "x1": x1,
                "y0": y0, "y1": y1,
                "n_descendants_in_panel": len([d for d in desc if d in y_map]),
                "all_descendants": desc,
            })

    if return_summary:
        return patches, pd.DataFrame(summaries).reset_index(drop=True)
    return patches




#################################################################################################









#### Define function for plotting individual Gene Conversion Events in a defined region ####

# ------------------------
# Constants (colors + sizes)
# ------------------------
# NT_COLORS = {"A": "green", "T": "red", "C": "blue", "G": "orange"}
# EVENT_BAR_COLOR = "#C9A0DC"   # light purple

# ------------------------
# Helper functions
# ------------------------



def add_cdseffect_legend(ax, missense_color="red", synonymous_color="blue"):
    """
    Add a CDS Effect legend to the given matplotlib axis.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        The axis to which the legend should be added.
    missense_color : str, optional
        Color for the 'Missense' legend entry (default: 'red').
    synonymous_color : str, optional
        Color for the 'Synonymous' legend entry (default: 'blue').
    """
    legend_elems = [
        Line2D([0], [0], color=missense_color, lw=2, label="Missense"),
        Line2D([0], [0], color=synonymous_color, lw=2, label="Synonymous"),
    ]
    ax.legend(
        handles=legend_elems,
        title="CDS Effect",
        fontsize=7,
        title_fontsize=8,
        loc="upper right",
    )


def _to_int_series(s: pd.Series) -> pd.Series:
    """Convert a Series to Int64 safely (NaNs preserved)."""
    return pd.to_numeric(s, errors="coerce").astype("Int64")


def prepare_snps_positions(snps_df: pd.DataFrame) -> pd.DataFrame:
    """
    Create a '__pos1__' column (1-based position) using Pos_0based or Pos_1based.
    Keeps all other columns intact.
    """
    out = snps_df.copy()
    if "Pos_0based" in out.columns:
        out["__pos1__"] = out["Pos_0based"].astype(int) + 1
    elif "Pos_1based" in out.columns:
        out["__pos1__"] = out["Pos_1based"].astype(int)
    else:
        raise ValueError("SNPs DF must include 'Pos_0based' or 'Pos_1based'.")
    return out

def _format_x_plain(in_ax):
    """Force full integer tick labels (no scientific notation, no offset)."""
    formatter = mticker.ScalarFormatter(useOffset=False)
    formatter.set_scientific(False)
    in_ax.xaxis.set_major_formatter(formatter)
    in_ax.ticklabel_format(style='plain', axis='x', useOffset=False)
    in_ax.xaxis.set_major_locator(mticker.MaxNLocator(nbins=6, integer=True))

def _draw_event_rect(ax, y_center: float, x0: float, x1: float, bar_height: float, color: str):
    """Draw an event as a filled rectangle centered at y_center with given height."""
    width = max(0.0, x1 - x0)
    rect = Rectangle((x0, y_center - bar_height/2), width, bar_height,
                     facecolor=color, edgecolor='none', alpha=0.9)
    ax.add_patch(rect)


# Reusable greedy row packer
def pack_intervals_greedy(starts: np.ndarray, ends: np.ndarray, min_gap: int = 0) -> np.ndarray:
    """
    Assign each interval [starts[i], ends[i]] to the first row where it doesn't overlap
    (with optional min_gap). Intervals should be pre-sorted by start (then end).
    Returns an array of row indices (0,1,2,...).
    """
    rows_last_end: list[float] = []
    rows: list[int] = []
    for s, e in zip(starts, ends):
        placed = False
        for ri in range(len(rows_last_end)):
            if s > (rows_last_end[ri] + min_gap):
                rows.append(ri)
                rows_last_end[ri] = e
                placed = True
                break
        if not placed:
            rows.append(len(rows_last_end))
            rows_last_end.append(e)
    return np.asarray(rows, dtype=int)



# ------------------------
# Plotting functions
# ------------------------

def plot_multi_event_snp_packed_on_ax_V2(
    ax,
    region_start1: int,
    region_end1:   int,
    pGCE_df: pd.DataFrame,
    snps_df: pd.DataFrame,
    *,
    min_rows: int | None = None,
    min_gap: int = 0,                 # bp gap required to share a row
    bar_height: float = 0.5,          # thickness of each event bar
    snp_line_width: float = 1,      # SNP tick thickness
    label_fontsize: int = 8,
    nt_colors: dict | None = None,
    event_bar_color: str = "#FF8488", # light red color, instead of #"#C9A0DC", # light purple
    event_edge_color: str = "black", # black edge
    event_bar_alpha: float = 0.6,
    hide_spines: tuple[str, ...] = ("left", "top", "right"),
    show_xlabel: bool = True,
    title: str | None = None,
    add_legend: bool = False,
    # Effect-coloring toggle
    color_by_CDSeffect: bool = False,    # False -> A/T/C/G colors; True -> missense/synonymous colors
    missense_color: str = "red",
    synonymous_color: str = "blue",  # blue
    add_ParalogMatchNames_ToLabel = False,
    # Bioframe options
    i_chrom: str = "NC_000962.3",
):
    """
    Packed view of GC events + SNP ticks over [region_start1, region_end1] (1-based, inclusive).

    Uses bioframe (bf) to select overlapping events, assuming:
      - pGCE_df has columns: start_0based (0-based), end_1based (1-based, inclusive), EventID
      - snps_df has: Pos_0based or Pos_1based, EventID, Child_Call, and (optionally) MissenseMut
    """
    if nt_colors is None:
        nt_colors = {"A":"green", "T":"red", "C":"blue", "G":"orange"}

    # ---- 1) Select overlapping events with bioframe
    
    region_start0 = region_start1 - 1
    target_region = f"{i_chrom}:{region_start0}-{region_end1}"

    pGCE_df = pGCE_df.copy()
    pGCE_df["chrom"] = pGCE_df["seqname"].copy()
    pGCE_df["start"] = pGCE_df["start_0based"].copy()
    pGCE_df["end"] = pGCE_df["end_1based"].copy()

    # pGCE_df = pGCE_df.copy()
    # pGCE_df.loc[:, "chrom"] = pGCE_df["seqname"].values
    # pGCE_df.loc[:, "start"] = pGCE_df["start_0based"].values
    # pGCE_df.loc[:, "end"] = pGCE_df["end_1based"].values

    
    events = bf.select(pGCE_df, target_region, cols = ("chrom", "start", "end")) 

    
    if events.empty:
        print("No events overlap the requested region.")
        return ax, None

    # Convert to 1-based inclusive for display & clip to window
    events.loc[:, "__start1__"]     = events["start"].astype(float) + 1.0
    events.loc[:, "__end1__"]       = events["end"].astype(float)
    events.loc[:, "__clip_start__"] = np.maximum(events["__start1__"].to_numpy(),
                                                float(region_start1))
    events.loc[:, "__clip_end__"]   = np.minimum(events["__end1__"].to_numpy(),
                                                float(region_end1))

    # Sort by clipped starts/ends
    events = events.sort_values(["__clip_start__", "__clip_end__"]).reset_index(drop=True)

    
    # ---- 2) Greedy row packing (re-usable utility)
    rows = pack_intervals_greedy(
        starts=events["__clip_start__"].to_numpy(dtype=float),
        ends=events["__clip_end__"].to_numpy(dtype=float),
        min_gap=min_gap
    )
    events["__row__"] = rows

    # Compute how many rows to *display*
    nrows_packed = int(events["__row__"].max()) + 1  # what packing actually needs
    nrows_display = max(nrows_packed, int(min_rows)) if (min_rows is not None) else nrows_packed



    # ---- 3) Prep SNPs (1-based positions)
    snps = snps_df.copy()
    if "Pos_0based" in snps.columns:
        snps["__pos1__"] = snps["Pos_0based"].astype(int) + 1
    elif "Pos_1based" in snps.columns:
        snps["__pos1__"] = snps["Pos_1based"].astype(int)
    else:
        raise ValueError("SNPs DF must include 'Pos_0based' or 'Pos_1based'.")
    snps_in_region = snps[(snps["__pos1__"] >= region_start1) & (snps["__pos1__"] <= region_end1)].copy()

    # ---- 4) Draw
    nrows = nrows_display #int(events["__row__"].max()) + 1
    for _, ev in events.iterrows():
        r = int(ev["__row__"])
        y = (nrows - 1 - r)  # row 0 at top
        x0 = float(ev["__clip_start__"])
        x1 = float(ev["__clip_end__"])

        # event rectangle
        rect = Rectangle((x0, y - bar_height/2), max(0.0, x1 - x0),
                         bar_height,
                         facecolor = event_bar_color,
                         edgecolor = event_edge_color,
                         alpha = event_bar_alpha,
                         linewidth=1.0)
        ax.add_patch(rect)

        
        # SNP ticks for this event
        s_ev = snps_in_region[snps_in_region["EventID"] == ev["EventID"]]
        if not s_ev.empty:
            ymin = y - bar_height/2
            ymax = y + bar_height/2
            for _, rS in s_ev.iterrows():
                xv = float(rS["__pos1__"])
                if color_by_CDSeffect:
                    is_missense = bool(rS.get("MissenseMut", False))
                    color = missense_color if is_missense else synonymous_color
                else:
                    nt = str(rS.get("Child_Call", "")).upper()[:1]
                    color = nt_colors.get(nt, "black")
                ax.vlines(x=xv, ymin=ymin, ymax=ymax, lw=snp_line_width, color=color, alpha=0.95)

        # label above bar

        event_label = str(ev["EventID"])

        if (ev["MappedEvent"] == True) & (add_ParalogMatchNames_ToLabel == True):
            Paralogs_WiBestKmerMatches = ev["Top_KmerMatch_HomologGeneIDs"] 

            event_label = str(ev["EventID"]) + " - " + str(Paralogs_WiBestKmerMatches)
        
        x_mid = (x0 + x1) / 2.0
        ax.text(x_mid, y + 0.7 * bar_height, event_label,
                ha="center", va="bottom", fontsize=label_fontsize)

    # ---- 5) Cosmetics
    ax.set_xlim(region_start1, region_end1)
    ax.set_ylim(-0.8, nrows - 1 + 0.8)
    ax.set_yticks([])
    if show_xlabel:
        ax.set_xlabel("Genomic position (1-based)")

    # Force full integer tick labels (no scientific notation, no offset).
    _format_x_plain(ax)


    ax.grid(axis="x", linestyle=":", linewidth=0.5, alpha=0.4)
    for spine in hide_spines:
        ax.spines[spine].set_visible(False)

    if add_legend:
        if color_by_CDSeffect:
            add_cdseffect_legend(ax, missense_color, synonymous_color)

        else:
            legend_elems = [Line2D([0], [0], color=nt_colors[k], lw=2, label=k)
                            for k in ["A","T","C","G"] if k in nt_colors]
            ax.legend(handles=legend_elems, title="Mut allele", fontsize=7, title_fontsize=8,
                      loc="upper right")

    if title is not None:
        ax.set_title(title, fontsize=10)

    summary = {
        "nrows": int(nrows),
        "events_plotted": events["EventID"].tolist(),
        "window": (region_start1, region_end1),
        "color_mode": "effect" if color_by_CDSeffect else "nt",
        "chrom": i_chrom,
    }
    return ax, summary





from typing import Optional, Tuple, Dict
import matplotlib.pyplot as plt

def plot_Anno_with_GCE_V3(
    Viz_Start: int,
    Viz_End: int,
    in_Genome_Graphic_Record,   # dna_features_viewer GraphicRecord (H37Rv)
    pGCE_df,                    # GC events dataframe
    snps_df,                    # SNPs dataframe (Pos_0based or Pos_1based, EventID, Child_Call)
    GCE_color_by_CDSeffect = False,
    ShowVariant_Legend = False,
    *,
    figsize: Tuple[float, float] = (15, 5), 
    dpi: int = 180,                            
    height_ratios: Tuple[float, float] = (1.0, 2.5),
    genome_label_threshold: int = 5,
    # GC-event track style
    
    gce_min_rows = None,
    gce_min_gap: int = 0,
    gce_bar_height: float = 0.5,
    gce_bar_color:  str = "#FF999C",
    gce_edge_color: str = "black",
    gce_snp_line_width: float = 1.2,
    gce_label_fontsize: int = 8,
    show_xlabel_on_gce: bool = True,
    add_ParalogMatchNames_ToGCELabel = False,
    title: Optional[str] = None,
    tight_layout: bool = True,
) -> Tuple[plt.Figure, Dict[str, plt.Axes]]:
    """
    Plot gene annotations (top) and packed GC events with SNPs (bottom).
    Returns (fig, {'Genome_Anno_ax':..., 'GCE_ax':...}).
    """
    # Crop genome record to region
    Graphic_Record_cropped = in_Genome_Graphic_Record.crop((Viz_Start, Viz_End + 1))

    # Figure & axes
    fig, axs = plt.subplots(
        2, 1, figsize=figsize, dpi=dpi,
        gridspec_kw={'height_ratios': height_ratios}
    )
    Genome_Anno_ax, GCE_ax = axs

    # Top: gene annotations
    if hasattr(Graphic_Record_cropped, "plot"):
        Graphic_Record_cropped.plot(
            strand_in_label_threshold=genome_label_threshold, ax=Genome_Anno_ax
        )
    Genome_Anno_ax.set_xlim(Viz_Start, Viz_End)
    Genome_Anno_ax.set_xticks([])
    Genome_Anno_ax.set_xlabel("")

    # Bottom: packed GC events + SNPs (uses your on-axis helper)
    plot_multi_event_snp_packed_on_ax_V2(
        GCE_ax,
        region_start1=Viz_Start, region_end1=Viz_End,
        pGCE_df=pGCE_df, snps_df=snps_df,
        min_gap=gce_min_gap,
        bar_height=gce_bar_height,
        event_bar_color = gce_bar_color,
        event_edge_color = gce_edge_color, 
        event_bar_alpha = 0.6,
        snp_line_width=gce_snp_line_width,
        label_fontsize=gce_label_fontsize,
        show_xlabel=show_xlabel_on_gce,
        add_legend = ShowVariant_Legend,
        title=None,
        color_by_CDSeffect = GCE_color_by_CDSeffect,
        add_ParalogMatchNames_ToLabel = add_ParalogMatchNames_ToGCELabel,
        min_rows=gce_min_rows,
    )

    if title:
        fig.suptitle(title, fontsize=11)

    if tight_layout:
        fig.tight_layout()

    return fig, {"Genome_Anno_ax": Genome_Anno_ax, "GCE_ax": GCE_ax}







from typing import Optional, Tuple, Dict
import matplotlib.pyplot as plt

def plot_ParalogVars_Anno_with_GCE(
    Viz_Start: int,
    Viz_End: int,
    in_Genome_Graphic_Record,   # dna_features_viewer GraphicRecord (H37Rv)
    pGCE_df,                    # GC events dataframe
    snps_df,                    # SNPs dataframe (Pos_0based or Pos_1based, EventID, Child_Call)
    ParalogVar_df,              # Paralog variants dataframe for plot_paralog_variants_on_ax
    GCE_color_by_CDSeffect: bool = False,
    ShowVariant_Legend: bool = False,
    *,
    figsize: Tuple[float, float] = (15, 7),
    dpi: int = 180,
    height_ratios: Tuple[float, float, float] = (1.4, 1.0, 2.5),
    genome_label_threshold: int = 5,
    # Paralog-variant track options
    paralog_bar_height: float = 0.36,
    paralog_snp_line_width: float = 0.9,
    paralog_show_labels: bool = True,
    paralog_title: Optional[str] = "Paralogous Region Alignments",
    # GCE-event track options
    gce_min_rows = None,
    gce_min_gap: int = 0,
    gce_bar_height: float = 0.5,
    gce_bar_color: str = "#FF999C",
    gce_edge_color: str = "black",
    gce_snp_line_width: float = 1.2,
    gce_label_fontsize: int = 8,
    show_xlabel_on_gce: bool = True,
    add_ParalogMatchNames_ToGCELabel: bool = False,
    gce_title: Optional[str] = None,
    anno_title: Optional[str] = None,
    title: Optional[str] = None,
    tight_layout: bool = True,
) -> Tuple[plt.Figure, Dict[str, plt.Axes]]:
    """
    Plot 3 stacked panels over [Viz_Start, Viz_End]:

      1) Paralog variants per paralog-region (top)
      2) Genome / gene annotations (middle)
      3) Packed GC events + SNPs (bottom)

    Returns
    -------
    fig : matplotlib.figure.Figure
    axes_dict : dict
        {
          "ParalogVar_ax": <Axes>,
          "Genome_Anno_ax": <Axes>,
          "GCE_ax": <Axes>,
        }
    """

    # --- Figure & axes: 3 stacked panels ---
    fig, axs = plt.subplots(
        3, 1,
        figsize=figsize,
        dpi=dpi,
        gridspec_kw={'height_ratios': height_ratios}
    )
    ParalogVar_ax, Genome_Anno_ax, GCE_ax = axs

    # --- Top: paralog variants ---
    plot_paralog_variants_on_ax(
        ParalogVar_ax,
        ParalogVar_df,
        region_start=Viz_Start,
        region_end=Viz_End,
        bar_height=paralog_bar_height,
        snp_line_width=paralog_snp_line_width,
        show_labels=paralog_show_labels,
        hide_xaxis=True,   # x-axis only on bottom panel
        title=paralog_title,
    )

    # --- Middle: genome annotations ---
    plot_genome_annotation_on_ax(
        Viz_Start=Viz_Start,
        Viz_End=Viz_End,
        in_Genome_Graphic_Record=in_Genome_Graphic_Record,
        i_ax=Genome_Anno_ax,
        genome_label_threshold=genome_label_threshold,
        top_title=anno_title,
    )

    # --- Bottom: packed GC events + SNPs ---
    plot_multi_event_snp_packed_on_ax_V2(
        GCE_ax,
        region_start1=Viz_Start,
        region_end1=Viz_End,
        pGCE_df=pGCE_df,
        snps_df=snps_df,
        min_gap=gce_min_gap,
        bar_height=gce_bar_height,
        event_bar_color = gce_bar_color,
        event_edge_color = gce_edge_color, # black edge
        event_bar_alpha = 0.6,
        snp_line_width=gce_snp_line_width,
        label_fontsize=gce_label_fontsize,
        show_xlabel=show_xlabel_on_gce,
        add_legend=ShowVariant_Legend,
        title=gce_title,
        color_by_CDSeffect=GCE_color_by_CDSeffect,
        add_ParalogMatchNames_ToLabel=add_ParalogMatchNames_ToGCELabel,
        min_rows=gce_min_rows,
    )

    # --- Overall figure title / layout ---
    if title:
        fig.suptitle(title, fontsize=11)

    if tight_layout:
        fig.tight_layout()

    axes_dict = {
        "ParalogVar_ax": ParalogVar_ax,
        "Genome_Anno_ax": Genome_Anno_ax,
        "GCE_ax": GCE_ax,
    }
    return fig, axes_dict


































############## Define `DNA-Features-Viewer` Preprocessing Functions ##############

#### A) Make graphic features for individual SNPs #####

def generate_SNP_GraphicFeatures(i_Event_SNPs):
    NtColor_Dict = {"A" : "green",
                    "T" : "red",
                    "C" : "blue", 
                    "G" : "orange" }

    L_SNP_GFeats = []
    
    for i, row in i_Event_SNPs.iterrows():
        
        pos_1 = row["Pos_1based"]
        Mut_Allele = row["Child_Call"]
    
        Nt_Color = NtColor_Dict[Mut_Allele]
    
        SNP_Feat = GraphicFeature(start = pos_1 , end = pos_1 + 1 ,
                                  #feature_level_height = 3,
                                  color = Nt_Color,
                                  linecolor = Nt_Color)
        
        L_SNP_GFeats.append(SNP_Feat)

    return L_SNP_GFeats


def generate_SNP_GraphicFeatures_HighlightMissense(i_Event_SNPs,
                                                   show_snp_labels = False):
    """
    Build GraphicFeature objects for SNPs.

    - Missense mutations: red
    - Synonymous (silent): grey
    """
    L_SNP_GFeats = []
    seen_labels = set()
    
    for _, row in i_Event_SNPs.iterrows():
        pos_1 = int(row["Pos_1based"])   # 1-based position
        is_missense = bool(row["MissenseMut"])
        
        # Color and label
        if is_missense:
            color = "red"
            label = None
            if show_snp_labels:
                candidate = f"{row.get('Ref_AA','_')}{int(row.get('Codon','_'))}{row.get('Mut_AA','_')}"  # e.g. A123T
                if candidate not in seen_labels:
                    label = candidate
                    seen_labels.add(candidate)
        else:
            color = "blue"
            label = None

        SNP_Feat = GraphicFeature(
            start=pos_1,
            end=pos_1 + 1,
            color=color,
            linecolor = color, #"grey",
            linewidth = 0.3,
            label=label
        )
        L_SNP_GFeats.append(SNP_Feat)

    return L_SNP_GFeats



def AddSNPs_HighlightMissense_ToGraphicRecord(Graphic_Record,
                                              i_Event_SNPs,
                                              missense_only=False,
                                              show_snp_labels=False):
    """
    Add SNP features to a GraphicRecord.

    - If missense_only=True, only add missense (non-synonymous) SNPs.
    """
    df = i_Event_SNPs
    if missense_only:
        df = df[df["MissenseMut"] == True]

    L_SNP_GFeats = generate_SNP_GraphicFeatures_HighlightMissense(df,
                                                                  show_snp_labels)
    
    Graphic_Record.features = Graphic_Record.features + L_SNP_GFeats

    return Graphic_Record



# def AddSNPs_ToGraphicRecord(Graphic_Record, i_Event_SNPs):

#     L_SNP_GFeats = generate_SNP_GraphicFeatures(i_Event_SNPs)
    
#     Graphic_Record.features = Graphic_Record.features + L_SNP_GFeats

#     return Graphic_Record


##################################################################################################




#### B) Make graphic features for Peptides/Epitopes + GC Events #####

def generate_GCEvent_GraphicFeatures(i_GC_Events_DF):

    L_GFeats = []
    
    for i, row in i_GC_Events_DF.iterrows():
        
        Rv_Start = row["start_0based"]
        Rv_End = row["end_1based"]
        i_Event_ID = row["EventID"]

        Event_Feat = GraphicFeature(start = Rv_Start , end = Rv_End,
                                      label = i_Event_ID[-3:],
                                      linewidth = 0.5,
                                      color = "#CBC3E3", #"purple",
                                      thickness = 5, 
                                      fontdict = {"fontsize": 4},
                                      linecolor = "black")
        
        L_GFeats.append(Event_Feat)

    return L_GFeats


def generate_Epitope_GraphicFeatures(i_Epitopes_DF):

    L_GFeats = []
    
    for i, row in i_Epitopes_DF.iterrows():
        
        Rv_Start = row["Rv_Start"]
        Rv_End = row["Rv_End"]
        i_epitope_seq = row["Epitope_Seq"]
        i_epitope_ID = row["Epitope_ID"]

        Epitope_Feat = GraphicFeature(start = Rv_Start , end = Rv_End ,
                                      #label = i_epitope_ID,
                                      color = "darkred",
                                      linecolor = "black")
        
        L_GFeats.append(Epitope_Feat)

    return L_GFeats


def generate_AssayedPeptide_GraphicFeatures(i_AssayedPeptides_DF):

    L_GFeats = []
    
    for i, row in i_AssayedPeptides_DF.iterrows():
        
        Rv_Start = row["Rv_Start"]
        Rv_End = row["Rv_End"]
        i_epitope_seq = row["Epitope_Seq"]
        i_epitope_ID = row["Epitope_ID"]
        i_PosEpitope = row["PosEpitope_Any"]
        if i_PosEpitope == True:
    
            Peptide_Feat = GraphicFeature(start = Rv_Start , end = Rv_End ,
                                          #label = i_epitope_ID,
                                          linewidth = 0.5,
                                          color = "darkred",
                                          thickness = 5,
                                          linecolor = "black")
            
        else:
            Peptide_Feat = GraphicFeature(start = Rv_Start , end = Rv_End ,
                                          #label = i_epitope_ID,
                                          #linewidth = 0.5,
                                          color = "lightgrey", alpha = 0.3,
                                          thickness = 5,
                                          linewidth = 0.5,
                                          linecolor = "black")

            
        L_GFeats.append(Peptide_Feat)

    return L_GFeats


#####################################################################################################





###### Functions for Composite viz of epitopes + mutations + gene conversio events + gene annotations ######

# `H37Rv Gene Annoations` + `Epitope Mapping`, + `mutation frequency` + `GC Events`, `etc`





def plot_Anno_Epitope_MutFreq(
    Viz_Start: int,
    Viz_End: int,
    in_PerCodon_MutCt_DF,                 # cols: Pos_0based, Mutation_Count
    in_Genome_Graphic_Record,             # dna_features_viewer GraphicRecord (H37Rv)
    in_AllPeptides_GraphicRecord,         # GraphicRecord of all assayed peptides
    in_EpitopeCov_DF,                     # cols: start, end, Cov_PosEpitopes
    *,
    mut_ylim: Optional[Tuple[float, float]] = (0, 15),
    figsize: Tuple[float, float] = (7, 3),
    height_ratios: Tuple[int, int, int, int] = (3, 5, 1, 4),
    cmap_name: str = "Reds",
    point_size: float = 2.0,
    vline_lw: float = 0.9,
    epi_linewidth: float = 2.5,
    genome_label_threshold: int = 5,
    show_bottom_xticks: bool = True,
    tight_layout: bool = True,
    epi_global_max: Optional[float] = None,):   # 👈 new option):
    """
    Plot genomic annotations, per-codon mutation counts, epitope coverage summary,
    and assayed peptides for a window [Viz_Start, Viz_End].

    Epitope coverage colors are normalized to the global max of Cov_PosEpitopes
    across the entire in_EpitopeCov_DF.
    """
    # --- Crop graphic records to region ---
    Graphic_Record_cropped = in_Genome_Graphic_Record.crop((Viz_Start, Viz_End + 1))
    AllPeptides_Records_cropped = in_AllPeptides_GraphicRecord.crop((Viz_Start, Viz_End + 1))

    # --- Create figure & axes ---
    fig, axs = plt.subplots(
        4, 1, figsize=figsize,
        gridspec_kw={'height_ratios': height_ratios}
    )
    Genome_Anno_ax, MF_ax, EpiSumm_ax, AllAssayedPeptides_ax = axs

    # --- Gene Annotations ---
    if hasattr(Graphic_Record_cropped, "plot"):
        Graphic_Record_cropped.plot(strand_in_label_threshold=genome_label_threshold,
                                    ax=Genome_Anno_ax)
    Genome_Anno_ax.set_xlim(Viz_Start, Viz_End)
    Genome_Anno_ax.set_xticks([])

    # --- Mutation Counts ---
    if {"Pos_0based", "Mutation_Count"}.issubset(in_PerCodon_MutCt_DF.columns):
        mut_df = in_PerCodon_MutCt_DF.query(
            f"(Pos_0based + 1) >= {Viz_Start} & (Pos_0based + 1) <= {Viz_End}"
        ).copy()
        if not mut_df.empty:
            xvals = mut_df["Pos_0based"].to_numpy() + 1
            yvals = mut_df["Mutation_Count"].to_numpy()
            MF_ax.vlines(x=xvals, ymin=0, ymax=yvals, color="black", linewidth=vline_lw, alpha=0.6)
            MF_ax.scatter(x=xvals, y=yvals, color="red", s=point_size, alpha=1)

    MF_ax.xaxis.set_major_formatter(mticker.ScalarFormatter())
    MF_ax.ticklabel_format(style='plain', axis='x')
    MF_ax.xaxis.set_major_locator(mticker.MaxNLocator(integer=True))
    MF_ax.set_xlim(Viz_Start, Viz_End)
    if mut_ylim is not None:
        MF_ax.set_ylim(*mut_ylim)
    MF_ax.set_xticks([])
    sns.despine(ax=MF_ax)

    # --- Epitope density summary ---
    cmap = plt.get_cmap(cmap_name)
    if {"start", "end", "Cov_PosEpitopes"}.issubset(in_EpitopeCov_DF.columns):
        Reg_EpiCov_DF = in_EpitopeCov_DF.query(
            f"start >= {Viz_Start} & end <= {Viz_End}"
        ).copy()
        if not Reg_EpiCov_DF.empty:
            if epi_global_max is None:
                global_max = in_EpitopeCov_DF["Cov_PosEpitopes"].max()
            else:
                global_max = epi_global_max
                
            norm = plt.Normalize(vmin=0, vmax=global_max if global_max > 0 else 1)
            for _, row in Reg_EpiCov_DF.iterrows():
                cov = row["Cov_PosEpitopes"]
                if cov > 0:
                    EpiSumm_ax.vlines(x=row["start"], ymin=0, ymax=1,
                                      color=cmap(norm(cov)), linewidth=epi_linewidth, alpha=1)

    EpiSumm_ax.set_xlim(Viz_Start, Viz_End)
    EpiSumm_ax.set_ylim(0, 1)
    EpiSumm_ax.set_xticks([])
    EpiSumm_ax.set_yticks([])
    for spine in ("left", "top", "right"):
        EpiSumm_ax.spines[spine].set_visible(False)
    EpiSumm_ax.set_xlabel("")

    # --- Assayed peptides ---
    #if hasattr(AllPeptides_Records_cropped, "plot"):
    AllPeptides_Records_cropped.plot(strand_in_label_threshold=genome_label_threshold,
                                         ax=AllAssayedPeptides_ax)
    
    AllAssayedPeptides_ax.set_xlim(Viz_Start, Viz_End)
    AllAssayedPeptides_ax.set_ylim(-0.5, 5.5)


    # full integers on x-axis (no sci/offset)
    formatter = mticker.ScalarFormatter(useOffset=False)
    formatter.set_scientific(False)
    AllAssayedPeptides_ax.xaxis.set_major_formatter(formatter)
    AllAssayedPeptides_ax.ticklabel_format(style='plain', axis='x', useOffset=False)
    AllAssayedPeptides_ax.xaxis.set_major_locator(mticker.MaxNLocator(nbins=6, integer=True))
    
    if tight_layout:
        fig.tight_layout()

    axes = {
        "Genome_Anno_ax": Genome_Anno_ax,
        "MF_ax": MF_ax,
        "EpiSumm_ax": EpiSumm_ax,
        "AllAssayedPeptides_ax": AllAssayedPeptides_ax,
    }
    return fig, axes










def plot_mut_epitope_region_with_gce_V2(
    Viz_Start: int,
    Viz_End: int,
    in_PerCodon_MutCt_DF,                 # cols: Pos_0based, Mutation_Count
    in_Genome_Graphic_Record,             # dna_features_viewer GraphicRecord (H37Rv)
    in_AllPeptides_GraphicRecord,         # GraphicRecord of all assayed peptides
    in_EpitopeCov_DF,                     # cols: start, end, Cov_PosEpitopes
    pGCE_df,          
    snps_df,
    *,
    mut_ylim: Optional[Tuple[float, float]] = (0, 15),
    figsize: Tuple[float, float] = (5, 5),
    height_ratios: Tuple[int, int, int, int, int] = (4, 5, 1, 4, 8),
    cmap_name: str = "Reds",
    point_size: float = 2.0,
    vline_lw: float = 0.9,
    epi_linewidth: float = 2.5,
    genome_label_threshold: int = 5,
    # GC-event track style
    gce_min_rows = None,
    GCE_color_by_CDSeffect = True,
    gce_min_gap: int = 0,
    gce_bar_height: float = 0.5,
    gce_snp_line_width: float = 1.2,
    gce_label_fontsize: int = 8,
    show_bottom_xticks: bool = True,
    tight_layout: bool = True,
    epi_global_max: Optional[float] = None,):
    """
    Same as plot_mut_epitope_region, but adds a GC events track at the bottom.

    Axes order (top→bottom):
      0: Genome_Anno_ax (gene annotations)
      1: MF_ax          (mutation counts)
      2: EpiSumm_ax     (epitope density summary bar)
      3: AllAssayedPeptides_ax (all assayed peptides)
      4: GCE_ax         (gene conversion events)
    """
    # --- Crop graphic records to region ---
    Graphic_Record_cropped   = in_Genome_Graphic_Record.crop((Viz_Start, Viz_End + 1))
    AllPeptides_Records_cropped = in_AllPeptides_GraphicRecord.crop((Viz_Start, Viz_End + 1))

    # --- Create figure & axes ---
    fig, axs = plt.subplots(
        5, 1, figsize=figsize,
        gridspec_kw={'height_ratios': height_ratios}
    )
    Genome_Anno_ax, MF_ax, EpiSumm_ax, AllAssayedPeptides_ax, GCE_ax = axs

    # --- Gene Annotations (top track) ---
    if hasattr(Graphic_Record_cropped, "plot"):
        Graphic_Record_cropped.plot(strand_in_label_threshold = genome_label_threshold,
                                    ax=Genome_Anno_ax)
    Genome_Anno_ax.set_xlim(Viz_Start, Viz_End)
    Genome_Anno_ax.set_xticks([])

    # --- Mutation Counts ---
    if {"Pos_0based", "Mutation_Count"}.issubset(in_PerCodon_MutCt_DF.columns):
        mut_df = in_PerCodon_MutCt_DF.query(
            f"(Pos_0based + 1) >= {Viz_Start} & (Pos_0based + 1) <= {Viz_End}"
        ).copy()
        if not mut_df.empty:
            xvals = mut_df["Pos_0based"].to_numpy() + 1
            yvals = mut_df["Mutation_Count"].to_numpy()
            MF_ax.vlines(x=xvals, ymin=0, ymax=yvals, color="black", linewidth=vline_lw, alpha=0.6)
            MF_ax.scatter(x=xvals, y=yvals, color="red", s=point_size, alpha=1)

    MF_ax.xaxis.set_major_formatter(mticker.ScalarFormatter())
    MF_ax.ticklabel_format(style='plain', axis='x')
    MF_ax.xaxis.set_major_locator(mticker.MaxNLocator(integer=True))
    MF_ax.set_xlim(Viz_Start, Viz_End)
    if mut_ylim is not None:
        MF_ax.set_ylim(*mut_ylim)
    MF_ax.set_xticks([])
    sns.despine(ax=MF_ax)

    # --- Epitope density summary (global-normalized) ---
    cmap = plt.get_cmap(cmap_name)
    if {"start", "end", "Cov_PosEpitopes"}.issubset(in_EpitopeCov_DF.columns):
        Reg_EpiCov_DF = in_EpitopeCov_DF.query(
            f"start >= {Viz_Start} & end <= {Viz_End}"
        ).copy()
        if not Reg_EpiCov_DF.empty:
            global_max = (in_EpitopeCov_DF["Cov_PosEpitopes"].max()
                          if epi_global_max is None else float(epi_global_max))
            norm = plt.Normalize(vmin=0, vmax=global_max if global_max > 0 else 1)
            for _, row in Reg_EpiCov_DF.iterrows():
                cov = row["Cov_PosEpitopes"]
                if cov > 0:
                    EpiSumm_ax.vlines(x=row["start"], ymin=0, ymax=1,
                                      color=cmap(norm(cov)), linewidth=epi_linewidth, alpha=1)

    EpiSumm_ax.set_xlim(Viz_Start, Viz_End)
    EpiSumm_ax.set_ylim(0, 1)
    EpiSumm_ax.set_xticks([])
    EpiSumm_ax.set_yticks([])
    for spine in ("left", "top", "right"):
        EpiSumm_ax.spines[spine].set_visible(False)
    EpiSumm_ax.set_xlabel("")


    # --- Assayed peptides ---
    if hasattr(AllPeptides_Records_cropped, "plot"):
        AllPeptides_Records_cropped.plot(strand_in_label_threshold=genome_label_threshold,
                                         ax=AllAssayedPeptides_ax)
    AllAssayedPeptides_ax.set_xlim(Viz_Start, Viz_End)
    AllAssayedPeptides_ax.set_ylim(-0.5, 5.5)
    AllAssayedPeptides_ax.set_xticks([])

    
    # --- GC events (bottom track) ---
    
    # plot_multi_event_snp_packed_on_ax(GCE_ax,
    #                                   region_start1 = Viz_Start, region_end1 = Viz_End,
    #                                   pGCE_df = pGCE_DF, snps_df = snps_df, label_fontsize = 4,
    #                                   add_legend=False, )

    plot_multi_event_snp_packed_on_ax_V2(GCE_ax,
                                      region_start1=Viz_Start, region_end1=Viz_End,
                                      pGCE_df = pGCE_df,
                                      snps_df = snps_df,
                                      min_gap=gce_min_gap,
                                      bar_height=gce_bar_height,
                                      snp_line_width=gce_snp_line_width,
                                      label_fontsize=gce_label_fontsize,
                                      color_by_CDSeffect = GCE_color_by_CDSeffect,
                                      #show_xlabel=show_xlabel_on_gce,
                                      add_legend=False,
                                      title=None,
                                      min_rows = gce_min_rows)

    # GCE_ax.set_xlim(Viz_Start, Viz_End)

    if show_bottom_xticks:
        GCE_ax.xaxis.set_major_formatter(mticker.ScalarFormatter())
        
        GCE_ax.xaxis.set_major_locator(mticker.MaxNLocator(integer=True))
        #GCE_ax.ticklabel_format(style='plain', axis='x')
        
        GCE_ax.ticklabel_format(
            axis="x",
            style="plain",
            useOffset=False
        )

    
    if tight_layout:
        fig.tight_layout()

    axes = {
        "Genome_Anno_ax": Genome_Anno_ax,
        "MF_ax": MF_ax,
        "EpiSumm_ax": EpiSumm_ax,
        "AllAssayedPeptides_ax": AllAssayedPeptides_ax,
        "GCE_ax": GCE_ax,
    }
    return fig, axes







































