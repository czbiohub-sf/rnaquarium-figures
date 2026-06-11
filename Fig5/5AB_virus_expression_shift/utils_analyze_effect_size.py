from __future__ import annotations 
import os, time, math
import numpy as np
import pandas as pd
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor, as_completed
import warnings, time, os, json, time
from typing import Dict, Any, Optional, Tuple, Iterable, Mapping, Callable, List, Literal
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import gaussian_kde
from scipy.stats import ks_2samp
from scipy import stats
from scipy.stats import anderson_ksamp, PermutationMethod
from matplotlib.patches import Patch
import re
from typing import Optional, Tuple, List


def read_parquet(parquet_path: str, index_col: Optional[str]) -> pd.DataFrame:
    df = pd.read_parquet(parquet_path)
    if index_col is not None:
        if index_col in df.columns:
            df = df.set_index(index_col)
        else:
            warnings.warn(
                f"index_col='{index_col}' not found in columns; using existing DataFrame index."
            )
    return df

def compute_pairwise_stat_ci_grid(
    list_a: Iterable,  # rows (e.g., genes)
    list_b: Iterable,  # cols (e.g., viruses)
    fetch_vals_fn: Callable[..., Tuple[np.ndarray, np.ndarray]],
    compute_fn: Callable[..., Mapping[str, object]],
    *,
    fetch_kwargs: dict | None = None,
    compute_kwargs: dict | None = None,
    ci_fmt: str = "[{:.2f}, {:.2f}]",
    progress: bool = True,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Build (stat_df, ci_text_df) over all combinations of list_a x list_b.
    - fetch_vals_fn(a, b, **fetch_kwargs) -> (pos_vals, neg_vals)
      (flip order here if your function returns differently)
    - compute_fn(neg_vals, pos_vals, **compute_kwargs) -> dict with:
        {"point": float, "ci": (low, high)}

    Returns
    -------
    stat_df : DataFrame of test statistics
    ci_df   : DataFrame of CI strings formatted by `ci_fmt`
    """
    fetch_kwargs = fetch_kwargs or {}
    compute_kwargs = compute_kwargs or {}

    list_a = list(list_a)
    list_b = list(list_b)

    stat = np.zeros((len(list_a), len(list_b)), dtype=float)
    ci_text = np.empty((len(list_a), len(list_b)), dtype=object)

    for i, a in enumerate(list_a):
        if progress:
            print(f"computing {a}...", end="")
        t0 = time.time()

        for j, b in enumerate(list_b):
            # Expecting (pos_vals, neg_vals) from fetch; swap if yours differs
            pos_vals, neg_vals = fetch_vals_fn(
                gene_name=a, virus_name=b, **fetch_kwargs
            )
            res = compute_fn(neg_vals, pos_vals, **compute_kwargs)
            s = float(res["point"])
            lo, hi = res["ci"]

            stat[i, j] = s
            ci_text[i, j] = ci_fmt.format(lo, hi)

        if progress:
            print(f"{time.time() - t0:.2f} seconds")

    stat_df = pd.DataFrame(stat, index=list_a, columns=list_b)
    ci_df = pd.DataFrame(ci_text, index=list_a, columns=list_b)
    return stat_df, ci_df



def plot_clustered_heatmap_with_ci(
    stat_df: pd.DataFrame,
    ci_df: pd.DataFrame,
    *,
    categories: Mapping[object, str] | None = None,   # categories for rows or cols
    categories_on: str = "rows",                      # "rows" or "cols"
    cmap: str = "vlag",

    # ---- text sizing ----
    annot_fontsize: float = 7,                        # CI text inside each cell
    row_tick_fontsize: float = 12,
    col_tick_fontsize: float = 12,                     # x labels (rotated 45°)
    axis_label_fontsize: float = 14,                  # "List A/B"
    legend_fontsize: float = 12,
    legend_title_fontsize: float = 12,
    cbar_label_fontsize: float = 12,
    cbar_tick_fontsize: float = 12,

    # ---- sizing/layout ----
    figsize: Tuple[float, float] | None = None,
    cell_size: Tuple[float, float] = (0.45, 0.38),    # width, height (inches) per cell
    dendrogram_ratio: Tuple[float, float] = (0.12, 0.12),
    cbar_label: str = "Test statistic",

    x_label: str = "Virus",
    y_label: str = "Gene",
    # ---- saving options ----
    savepath: str | Path | None = None,               # e.g., "out.png", "out.pdf", "out.svg"
    dpi: int = 300,
    bbox_inches: str = "tight",
    pad_inches: float = 0.02,
    transparent: bool = True,
    rasterize_heatmap: bool = False,                  # useful for compact vector PDFs
    close: bool = False,                              # close the figure after saving (no show)
):
    """
    Make a seaborn clustermap of `stat_df` and overlay CI text from `ci_df`.

    Adjustable sizing:
    - If `figsize` is given, it is used as-is.
    - Otherwise, the size is computed from `cell_size` and the grid shape.

    All text is adjustable via the *_fontsize arguments.
    X-axis tick labels are rotated 45 degrees and right-aligned.

    Returns
    -------
    cg : seaborn.matrix.ClusterGrid
        The clustermap object for further customization if needed.
    """
    # ---- compute figure size if not provided ----
    if figsize is None:
        n_rows, n_cols = stat_df.shape
        # Padding for dendrograms/labels/legend/colorbar
        pad_w, pad_h = 3.5, 2.8
        width = max(6.0, n_cols * cell_size[0] + pad_w)
        height = max(5.0, n_rows * cell_size[1] + pad_h)
        figsize = (width, height)

    # ---- category color bars ----
    row_colors = col_colors = None
    legend_handles = []
    if categories is not None:
        if categories_on == "rows":
            labels = [categories.get(idx, "Unlabeled") for idx in stat_df.index]
            cats = pd.Categorical(labels)
            pal = sns.color_palette("Set2", n_colors=len(cats.categories))
            lut: Dict[str, tuple] = dict(zip(cats.categories, pal))
            row_colors = pd.Series(labels, index=stat_df.index).map(lut)
            legend_handles = [Patch(facecolor=lut[c], label=c) for c in cats.categories]
        else:
            labels = [categories.get(col, "Unlabeled") for col in stat_df.columns]
            cats = pd.Categorical(labels)
            pal = sns.color_palette("Set2", n_colors=len(cats.categories))
            lut: Dict[str, tuple] = dict(zip(cats.categories, pal))
            col_colors = pd.Series(labels, index=stat_df.columns).map(lut)
            legend_handles = [Patch(facecolor=lut[c], label=c) for c in cats.categories]

    # ---- robust color scaling centered at 0 when appropriate ----
    vmin = np.nanpercentile(stat_df.to_numpy(), 2.5)
    vmax = np.nanpercentile(stat_df.to_numpy(), 97.5)
    center = 0.0 if (vmin < 0 < vmax) else None

    cg = sns.clustermap(
        stat_df,
        cmap=cmap,
        center=center,
        vmin=None if center is None else vmin,
        vmax=None if center is None else vmax,
        linewidths=0.5,
        linecolor="white",
        row_colors=row_colors,
        col_colors=col_colors,
        cbar_kws={"label": cbar_label},
        dendrogram_ratio=dendrogram_ratio,
        figsize=figsize,
    )

    # ---- legend for categories ----
    ax_for_legend = cg.ax_row_dendrogram if categories_on == "rows" else cg.ax_col_dendrogram
    if legend_handles:
        leg = ax_for_legend.legend(
            handles=legend_handles,
            title="Category",
            loc="center",
            bbox_to_anchor=(0.5, 0.5),
            frameon=False,
            prop={"size": legend_fontsize},
        )
        if leg.get_title() is not None:
            leg.get_title().set_fontsize(legend_title_fontsize)

    # ---- reorder CI to clustered order & annotate ----
    row_order = [stat_df.index[i] for i in cg.dendrogram_row.reordered_ind]
    col_order = [stat_df.columns[j] for j in cg.dendrogram_col.reordered_ind]
    ci_reordered = ci_df.loc[row_order, col_order].to_numpy()

    ax = cg.ax_heatmap
    n_rows, n_cols = ci_reordered.shape
    for i in range(n_rows):
        for j in range(n_cols):
            ax.text(
                j + 0.5,
                i + 0.5,
                str(ci_reordered[i, j]),
                ha="center",
                va="center",
                fontsize=annot_fontsize,
            )

    # ---- axis labels ----
    ax.set_xlabel(x_label, fontsize=axis_label_fontsize)
    ax.set_ylabel(y_label, fontsize=axis_label_fontsize)

    # ---- tick label fonts & rotation ----
    # x: rotate 45°, right-align for readability
    ax.set_xticklabels(
        ax.get_xticklabels(),
        rotation=45,
        ha="right",
        rotation_mode="anchor",
        fontsize=col_tick_fontsize,
    )
    # y: keep horizontal
    ax.set_yticklabels(
        ax.get_yticklabels(),
        rotation=0,
        fontsize=row_tick_fontsize,
    )

    # ---- colorbar fonts ----
    if cg.cax is not None:
        cg.cax.set_ylabel(cbar_label, fontsize=cbar_label_fontsize)
        cg.cax.tick_params(labelsize=cbar_tick_fontsize)

    # ---- saving (optional) ----
    if savepath is not None:
        # Optionally rasterize the QuadMesh to keep vector files small
        if rasterize_heatmap:
            for coll in getattr(cg.ax_heatmap, "collections", []):
                try:
                    coll.set_rasterized(True)
                except Exception:
                    pass

        fig = getattr(cg, "figure", None) or getattr(cg, "fig", None) or plt.gcf()
        fig.savefig(
            str(savepath),
            dpi=dpi,
            bbox_inches=bbox_inches,
            pad_inches=pad_inches,
            transparent=transparent,
        )
        if close:
            plt.close(fig)
            return cg  # don't call plt.show() if we just closed

    plt.show()
    return cg


def plot_clustered_heatmap_with_ci_and_violin(
    stat_df: pd.DataFrame,
    ci_df: pd.DataFrame,
    *,
    categories: Mapping[object, str] | None = None,   # categories for rows or cols
    categories_on: str = "rows",                      # "rows" or "cols"
    cmap: str = "vlag",

    # ---- text sizing ----
    annot_fontsize: float = 7,                        # CI text inside each cell
    row_tick_fontsize: float = 12,
    col_tick_fontsize: float = 12,                    # x labels (rotated 45°)
    axis_label_fontsize: float = 14,                  # "List A/B"
    legend_fontsize: float = 12,
    legend_title_fontsize: float = 12,
    cbar_label_fontsize: float = 12,
    cbar_tick_fontsize: float = 12,

    # ---- sizing/layout ----
    figsize: Tuple[float, float] | None = None,
    cell_size: Tuple[float, float] = (0.45, 0.38),    # width, height (inches) per cell
    dendrogram_ratio: Tuple[float, float] = (0.12, 0.12),
    cbar_label: str = "Test statistic",

    x_label: str = "Virus",
    y_label: str = "Gene",

    # ---- mini-violin controls ----
    # Expect: {col_name: pandas.DataFrame with a numeric "point" column}
    col_nulls: Mapping[object, pd.DataFrame] | None = None,
    violin_height_ratio: float = 0.6,   # fraction of the col-dendrogram height to devote to violins
    violin_width: float = 0.9,          # relative width per column cell (0..1-ish)
    violin_alpha: float = 0.65,
    violin_facecolor: str = "0.4",      # grayscale ok
    violin_edgecolor: str = "none",
    violin_zero: float | None = 0.0,    # draw dotted zero line if not None

    # ---- saving options ----
    savepath: str | Path | None = None,               # e.g., "out.png", "out.pdf", "out.svg"
    dpi: int = 300,
    bbox_inches: str = "tight",
    pad_inches: float = 0.02,
    transparent: bool = True,
    rasterize_heatmap: bool = False,                  # useful for compact vector PDFs
    close: bool = False,                              # close the figure after saving (no show)
):
    """
    Make a seaborn clustermap of `stat_df` and overlay CI text from `ci_df`.
    Add mini-violins of each column's null distribution above the column dendrogram.

    Color scaling:
      - vmax comes from the actual heatmap matrix (stat_df).
      - vmin comes from the global minimum across all col_nulls[col]['point'] values.
      - Both heatmap and mini-violins share [vmin, vmax].

    Violin clipping:
      - Null samples used for violins are clipped by excluding points > vmax
        (so KDE/shape never extends beyond the heatmap's max).
      - The violin PolyCollections are also hard-clipped to the axes patch.
    """

    # ---- figure size ----
    if figsize is None:
        n_rows, n_cols = stat_df.shape
        pad_w, pad_h = 3.5, 2.8
        width = max(6.0, n_cols * cell_size[0] + pad_w)
        height = max(5.0, n_rows * cell_size[1] + pad_h)
        figsize = (width, height)

    # ---- category color bars ----
    row_colors = col_colors = None
    legend_handles = []
    if categories is not None:
        if categories_on == "rows":
            labels = [categories.get(idx, "Unlabeled") for idx in stat_df.index]
            cats = pd.Categorical(labels)
            pal = sns.color_palette("Set2", n_colors=len(cats.categories))
            lut: Dict[str, tuple] = dict(zip(cats.categories, pal))
            row_colors = pd.Series(labels, index=stat_df.index).map(lut)
            legend_handles = [Patch(facecolor=lut[c], label=c) for c in cats.categories]
        else:
            labels = [categories.get(col, "Unlabeled") for col in stat_df.columns]
            cats = pd.Categorical(labels)
            pal = sns.color_palette("Set2", n_colors=len(cats.categories))
            lut: Dict[str, tuple] = dict(zip(cats.categories, pal))
            col_colors = pd.Series(labels, index=stat_df.columns).map(lut)
            legend_handles = [Patch(facecolor=lut[c], label=c) for c in cats.categories]

    # ---- vmax from stat_df; vmin from col_nulls ----
    vmax = float(np.nanmax(stat_df.to_numpy()))
    if not np.isfinite(vmax):
        vmax = 1.0

    vmin = None
    if col_nulls:
        mins = []
        for k, df in col_nulls.items():
            if isinstance(df, pd.DataFrame) and ("point" in df.columns):
                vals = pd.to_numeric(df["point"], errors="coerce").to_numpy()
                vals = vals[np.isfinite(vals)]
                if vals.size:
                    mins.append(np.min(vals))
        if mins:
            vmin = float(np.min(mins))
    if vmin is None or not np.isfinite(vmin):
        vmin = float(np.nanpercentile(stat_df.to_numpy(), 2.5))
        if not np.isfinite(vmin):
            vmin = -1.0

    if vmin >= vmax:
        pad = 1.0 if vmax == 0 else abs(vmax) * 0.05
        vmin = vmax - pad
    center = 0.0 if (vmin < 0 < vmax) else None

    # ---- clustermap with fixed [vmin, vmax] ----
    cg = sns.clustermap(
        stat_df,
        cmap=cmap,
        center=center,
        vmin=vmin,
        vmax=vmax,
        linewidths=0.5,
        linecolor="white",
        row_colors=row_colors,
        col_colors=col_colors,
        cbar_kws={"label": cbar_label},
        dendrogram_ratio=dendrogram_ratio,
        figsize=figsize,
    )

    # ---- legend ----
    ax_for_legend = cg.ax_row_dendrogram if categories_on == "rows" else cg.ax_col_dendrogram
    if legend_handles:
        leg = ax_for_legend.legend(
            handles=legend_handles,
            title="Category",
            loc="center",
            bbox_to_anchor=(0.5, 0.5),
            frameon=False,
            prop={"size": legend_fontsize},
        )
        if leg.get_title() is not None:
            leg.get_title().set_fontsize(legend_title_fontsize)

    # ---- annotate CI in clustered order ----
    row_order = [stat_df.index[i] for i in cg.dendrogram_row.reordered_ind]
    col_order = [plot_clustered_heatmap_with_ci_and_violin.columns[j] for j in cg.dendrogram_col.reordered_ind]
    ci_reordered = ci_df.loc[row_order, col_order].to_numpy()

    ax = cg.ax_heatmap
    n_rows, n_cols = ci_reordered.shape
    for i in range(n_rows):
        for j in range(n_cols):
            ax.text(
                j + 0.5,
                i + 0.5,
                str(ci_reordered[i, j]),
                ha="center",
                va="center",
                fontsize=annot_fontsize,
            )

    # =========================
    # Mini-violins above columns (clipped at vmax)
    # =========================
    if col_nulls is not None and len(col_nulls) > 0:
        datasets = []
        positions = []
        eps = max(1e-9, (vmax - vmin) * 1e-9)  # tiny epsilon for jitter fallback

        for j, c in enumerate(col_order):
            df = col_nulls.get(c, None)
            if isinstance(df, pd.DataFrame) and ("point" in df.columns):
                arr0 = pd.to_numeric(df["point"], errors="coerce").to_numpy()
                arr0 = arr0[np.isfinite(arr0)]
                if arr0.size == 0:
                    continue
                mu0 = float(np.nanmean(arr0))

                # ---- clip dataset by excluding points above vmax ----
                arr = arr0[arr0 <= vmax]

                if arr.size >= 3 and np.nanstd(arr) > 0:
                    datasets.append(arr)
                    positions.append(j + 0.5)
                else:
                    # If everything was > vmax (or too few points), place a tiny jitter at (vmax - eps)
                    y = min(max(vmin, min(mu0, vmax - eps)), vmax - eps)
                    jitter = np.array([y - eps, y, y + eps])
                    datasets.append(jitter)
                    positions.append(j + 0.5)

        # Split the existing column-dendrogram axis into [dendrogram | violins]
        fig = cg.fig if hasattr(cg, "fig") else cg.figure
        den_bbox = cg.ax_col_dendrogram.get_position()
        h_viol = den_bbox.height * float(violin_height_ratio)
        h_dend = den_bbox.height - h_viol
        cg.ax_col_dendrogram.set_position([den_bbox.x0, den_bbox.y0, den_bbox.width, h_dend])
        ax_violin = fig.add_axes([den_bbox.x0, den_bbox.y0 + h_dend, den_bbox.width, h_viol])

        if datasets:
            vp = ax_violin.violinplot(
                datasets,
                positions=positions,
                widths=violin_width,
                showmeans=False, showmedians=False, showextrema=False,
            )
            # style + hard clipping to axis patch
            for body in vp.get("bodies", []):
                body.set_facecolor(violin_facecolor)
                body.set_edgecolor(violin_edgecolor)
                body.set_alpha(violin_alpha)
                body.set_clip_on(True)
                body.set_clip_path(ax_violin.patch)

        # Zero/reference line only if within [vmin, vmax]
        if violin_zero is not None and np.isfinite(violin_zero) and (vmin <= violin_zero <= vmax):
            ax_violin.axhline(float(violin_zero), lw=0.6, ls=":", color="k", alpha=0.6, zorder=1)

        # Align with heatmap coordinates and hide ticks/frames
        ax_violin.set_xlim(0, n_cols)
        ax_violin.set_ylim(vmin, vmax)
        ax_violin.set_xticks([])
        ax_violin.set_yticks([])
        ax_violin.set_frame_on(False)

    # ---- labels, ticks, colorbar ----
    ax.set_xlabel(x_label, fontsize=axis_label_fontsize)
    ax.set_ylabel(y_label, fontsize=axis_label_fontsize)
    ax.set_xticklabels(
        ax.get_xticklabels(),
        rotation=45, ha="right", rotation_mode="anchor",
        fontsize=col_tick_fontsize,
    )
    ax.set_yticklabels(
        ax.get_yticklabels(),
        rotation=0,
        fontsize=row_tick_fontsize,
    )
    if cg.cax is not None:
        cg.cax.set_ylabel(cbar_label, fontsize=cbar_label_fontsize)
        cg.cax.tick_params(labelsize=cbar_tick_fontsize)

    # ---- saving ----
    if savepath is not None:
        if rasterize_heatmap:
            for coll in getattr(cg.ax_heatmap, "collections", []):
                try:
                    coll.set_rasterized(True)
                except Exception:
                    pass
        fig = getattr(cg, "figure", None) or getattr(cg, "fig", None) or plt.gcf()
        fig.savefig(
            str(savepath),
            dpi=dpi,
            bbox_inches=bbox_inches,
            pad_inches=pad_inches,
            transparent=transparent,
        )
        if close:
            plt.close(fig)
            return cg

    plt.show()
    return cg



def plot_clustered_heatmap_with_ci_and_violin_and_symbols(
    stat_df: pd.DataFrame,
    ci_df: pd.DataFrame,
    *,
    # ---- optional category color bars ----
    # Backward-compatible: `categories` + `categories_on` will populate row OR col categories
    categories: Mapping[object, str] | None = None,
    categories_on: str = "rows",  # "rows" or "cols"
    # New: you can provide BOTH at once
    row_categories: Mapping[object, str] | None = None,
    col_categories: Mapping[object, str] | None = None,
    row_category_palette: str = "Set2",
    col_category_palette: str = "Set2",
    row_category_legend_title: str = "Row category",
    col_category_legend_title: str = "Column category",

    cmap: str = "vlag",

    # ---- text sizing ----
    annot_fontsize: float = 7,
    row_tick_fontsize: float = 12,
    col_tick_fontsize: float = 12,
    axis_label_fontsize: float = 14,
    legend_fontsize: float = 12,
    legend_title_fontsize: float = 12,
    cbar_label_fontsize: float = 12,
    cbar_tick_fontsize: float = 12,

    # ---- sizing/layout ----
    figsize: Tuple[float, float] | None = None,
    cell_size: Tuple[float, float] = (0.45, 0.38),
    dendrogram_ratio: Tuple[float, float] = (0.12, 0.12),
    cbar_label: str = "Test statistic",

    x_label: str = "Virus",
    y_label: str = "Gene",

    # ---- mini-violin controls ----
    col_nulls: Mapping[object, pd.DataFrame] | None = None,
    violin_height_ratio: float = 0.6,
    violin_width: float = 0.9,
    violin_alpha: float = 0.65,
    violin_facecolor: str = "0.4",
    violin_edgecolor: str = "none",
    violin_zero: float | None = 0.0,

    # ---- triangle annotations (per-column thresholds) ----
    upper_thresholds: Mapping[object, float] | None = None,
    upper_symbol: str = "▲",
    upper_fontsize: float | None = None,
    upper_color: str = "k",
    upper_pos: Tuple[float, float] = (0.82, 0.22),

    lower_thresholds: Mapping[object, float] | None = None,
    lower_symbol: str = "▼",
    lower_fontsize: float | None = None,
    lower_color: str = "k",
    lower_pos: Tuple[float, float] = (0.18, 0.78),

    # ---- saving options ----
    savepath: str | Path | None = None,
    dpi: int = 300,
    bbox_inches: str = "tight",
    pad_inches: float = 0.02,
    transparent: bool = True,
    rasterize_heatmap: bool = False,
    close: bool = False,

    # ---- extra x-axis increment labels ----
    show_col_increments: bool = True,
    col_increment_start: int = 1,
    increment_tick_fontsize: float | None = None,
    increment_tick_pad: float = 2.0,
    col_tick_pad: float = 18.0,

    # ---- CI text ----
    show_ci_text: bool = True,

    # ---- sizing/layout ----
    colors_ratio: tuple[float, float] | float = (0.03, 0.10),  # (row_colors, col_colors)
):
    """
    Make a seaborn clustermap of `stat_df` and overlay CI text from `ci_df`.
    Add mini-violins of each column's null distribution above the column dendrogram.

    Optional category color bars:
      - Backward compatible: pass `categories` and `categories_on="rows"/"cols"`.
      - New: pass `row_categories` and/or `col_categories` simultaneously.

    Triangle annotations:
      - If `upper_thresholds` is provided, draw an upward triangle where stat_df[row, col] > upper_thresholds[col].
      - If `lower_thresholds` is provided, draw a downward triangle where stat_df[row, col] < lower_thresholds[col].
      - Marker positions use `upper_pos` and `lower_pos` within each cell (0..1 in x and y).

    Color scaling:
      - vmax comes from stat_df.
      - vmin comes from global minimum across col_nulls[col]['point'] if provided; else from stat_df percentile.
      - Both heatmap and mini-violins share [vmin, vmax].

    Violin clipping:
      - Null samples used for violins are clipped by excluding points > vmax.
      - The violin PolyCollections are hard-clipped to the axes patch.
    """

    def _is_nonempty_mapping(m) -> bool:
        if m is None:
            return False
        try:
            return len(m) > 0
        except Exception:
            return True  # treat non-sized mapping-like objects as enabled

    def _coerce_label(x) -> str:
        # Normalize None/NaN into a stable label and cast everything to str for consistent legend labels
        if x is None:
            return "Unlabeled"
        try:
            if isinstance(x, float) and np.isnan(x):
                return "Unlabeled"
        except Exception:
            pass
        s = str(x)
        return "Unlabeled" if s.strip() == "" else s

    def _make_colorbar(keys, cat_map, palette_name: str):
        labels = [_coerce_label(cat_map.get(k, "Unlabeled")) for k in keys]
        cats = pd.Categorical(labels)
        pal = sns.color_palette(palette_name, n_colors=len(cats.categories)) if len(cats.categories) else []
        lut: Dict[str, tuple] = dict(zip(cats.categories, pal))
        colors = pd.Series(labels, index=keys).map(lut)
        handles = [Patch(facecolor=lut[c], label=str(c)) for c in cats.categories]
        return colors, handles

    # ---- figure size ----
    if figsize is None:
        n_rows0, n_cols0 = stat_df.shape
        pad_w, pad_h = 3.5, 2.8
        width = max(6.0, n_cols0 * cell_size[0] + pad_w)
        height = max(5.0, n_rows0 * cell_size[1] + pad_h)
        figsize = (width, height)

    # ---- category color bars (optional; row and/or col) ----
    # Map backward-compatible `categories` into row/col if user didn't provide row_categories/col_categories explicitly.
    if _is_nonempty_mapping(categories):
        cat_on = str(categories_on).lower()
        if cat_on in {"rows", "row"} and row_categories is None:
            row_categories = categories
        elif cat_on in {"cols", "col", "columns", "column"} and col_categories is None:
            col_categories = categories

    row_colors = col_colors = None
    row_legend_handles: list = []
    col_legend_handles: list = []

    if _is_nonempty_mapping(row_categories):
        row_colors, row_legend_handles = _make_colorbar(stat_df.index, row_categories, row_category_palette)

    if _is_nonempty_mapping(col_categories):
        col_colors, col_legend_handles = _make_colorbar(stat_df.columns, col_categories, col_category_palette)

    # ---- vmax from stat_df; vmin from col_nulls ----
    vmax = float(np.nanmax(stat_df.to_numpy()))
    if not np.isfinite(vmax):
        vmax = 1.0

    vmin = None
    if col_nulls:
        mins = []
        for _, df in col_nulls.items():
            if isinstance(df, pd.DataFrame) and ("point" in df.columns):
                vals = pd.to_numeric(df["point"], errors="coerce").to_numpy()
                vals = vals[np.isfinite(vals)]
                if vals.size:
                    mins.append(np.min(vals))
        if mins:
            vmin = float(np.min(mins))

    if vmin is None or not np.isfinite(vmin):
        vmin = float(np.nanpercentile(stat_df.to_numpy(), 2.5))
        if not np.isfinite(vmin):
            vmin = -1.0

    if vmin >= vmax:
        pad = 1.0 if vmax == 0 else abs(vmax) * 0.05
        vmin = vmax - pad

    center = 0.0 if (vmin < 0 < vmax) else None

    # ---- clustermap with fixed [vmin, vmax] ----
    clustermap_kwargs = dict(
        cmap=cmap,
        center=center,
        vmin=vmin,
        vmax=vmax,
        linewidths=0.5,
        linecolor="white",
        cbar_kws={"label": cbar_label},
        dendrogram_ratio=dendrogram_ratio,
        colors_ratio=colors_ratio,
        figsize=figsize,
    )
    if row_colors is not None:
        clustermap_kwargs["row_colors"] = row_colors
    if col_colors is not None:
        clustermap_kwargs["col_colors"] = col_colors

    cg = sns.clustermap(stat_df, **clustermap_kwargs)

    # ---- clustered order ----
    row_order = [stat_df.index[i] for i in cg.dendrogram_row.reordered_ind]
    col_order = [stat_df.columns[j] for j in cg.dendrogram_col.reordered_ind]
    ci_reordered = ci_df.loc[row_order, col_order].to_numpy()
    stat_reordered = stat_df.loc[row_order, col_order].to_numpy()

    ax = cg.ax_heatmap
    n_rows, n_cols = ci_reordered.shape

    # ---- annotate CI + triangles ----
    if show_ci_text:
        for i in range(n_rows):
            for j in range(n_cols):
                ax.text(
                    j + 0.5,
                    i + 0.5,
                    str(ci_reordered[i, j]),
                    ha="center",
                    va="center",
                    fontsize=annot_fontsize,
                    clip_on=True,
                    zorder=4,
                )

    if upper_thresholds is not None:
        sfs_up = (annot_fontsize + 2) if (upper_fontsize is None) else upper_fontsize
        sx, sy = upper_pos
        for i in range(n_rows):
            for j in range(n_cols):
                thr_up = upper_thresholds.get(col_order[j], None)
                if thr_up is None:
                    continue
                val = stat_reordered[i, j]
                if np.isfinite(val) and (val > float(thr_up)):
                    ax.text(
                        j + sx,
                        i + sy,
                        upper_symbol,
                        ha="center",
                        va="center",
                        fontsize=sfs_up,
                        color=upper_color,
                        clip_on=True,
                        zorder=5,
                    )

    if lower_thresholds is not None:
        sfs_dn = (annot_fontsize + 2) if (lower_fontsize is None) else lower_fontsize
        dx, dy = lower_pos
        for i in range(n_rows):
            for j in range(n_cols):
                thr_dn = lower_thresholds.get(col_order[j], None)
                if thr_dn is None:
                    continue
                val = stat_reordered[i, j]
                if np.isfinite(val) and (val < float(thr_dn)):
                    ax.text(
                        j + dx,
                        i + dy,
                        lower_symbol,
                        ha="center",
                        va="center",
                        fontsize=sfs_dn,
                        color=lower_color,
                        clip_on=True,
                        zorder=5,
                    )

    # =========================
    # Mini-violins above columns (clipped at vmax)
    # =========================
    ax_violin = None
    if col_nulls is not None and len(col_nulls) > 0:
        datasets = []
        positions = []
        eps = max(1e-9, (vmax - vmin) * 1e-9)

        for j, c in enumerate(col_order):
            df = col_nulls.get(c, None)
            if isinstance(df, pd.DataFrame) and ("point" in df.columns):
                arr0 = pd.to_numeric(df["point"], errors="coerce").to_numpy()
                arr0 = arr0[np.isfinite(arr0)]
                if arr0.size == 0:
                    continue
                mu0 = float(np.nanmean(arr0))

                arr = arr0[arr0 <= vmax]
                if arr.size >= 3 and np.nanstd(arr) > 0:
                    datasets.append(arr)
                    positions.append(j + 0.5)
                else:
                    y = min(max(vmin, min(mu0, vmax - eps)), vmax - eps)
                    jitter = np.array([y - eps, y, y + eps])
                    datasets.append(jitter)
                    positions.append(j + 0.5)

        fig = cg.fig if hasattr(cg, "fig") else cg.figure
        den_bbox = cg.ax_col_dendrogram.get_position()
        h_viol = den_bbox.height * float(violin_height_ratio)
        h_dend = den_bbox.height - h_viol
        cg.ax_col_dendrogram.set_position([den_bbox.x0, den_bbox.y0, den_bbox.width, h_dend])
        ax_violin = fig.add_axes([den_bbox.x0, den_bbox.y0 + h_dend, den_bbox.width, h_viol])

        if datasets:
            vp = ax_violin.violinplot(
                datasets,
                positions=positions,
                widths=violin_width,
                showmeans=False,
                showmedians=False,
                showextrema=False,
            )
            for body in vp.get("bodies", []):
                body.set_facecolor(violin_facecolor)
                body.set_edgecolor(violin_edgecolor)
                body.set_alpha(violin_alpha)
                body.set_clip_on(True)
                body.set_clip_path(ax_violin.patch)

        if violin_zero is not None and np.isfinite(violin_zero) and (vmin <= violin_zero <= vmax):
            ax_violin.axhline(float(violin_zero), lw=0.6, ls=":", color="k", alpha=0.6, zorder=1)

        ax_violin.set_xlim(0, n_cols)
        ax_violin.set_ylim(vmin, vmax)
        ax_violin.set_xticks([])
        ax_violin.set_yticks([])
        ax_violin.set_frame_on(False)

    # ---- legends (row + column, independently) ----
    if row_legend_handles:
        leg = cg.ax_row_dendrogram.legend(
            handles=row_legend_handles,
            title=row_category_legend_title,
            loc="center",
            bbox_to_anchor=(0.5, 0.5),
            frameon=False,
            prop={"size": legend_fontsize},
        )
        if leg.get_title() is not None:
            leg.get_title().set_fontsize(legend_title_fontsize)

    if col_legend_handles:
        # Put the column-category legend on the violin axis if it exists (usually cleaner),
        # otherwise on the column dendrogram axis.
        ax_for_col_legend = ax_violin if ax_violin is not None else cg.ax_col_dendrogram
        ncol = min(len(col_legend_handles), 6)
        leg = ax_for_col_legend.legend(
            handles=col_legend_handles,
            title=col_category_legend_title,
            loc="upper center",
            bbox_to_anchor=(0.5, 1.0),
            frameon=False,
            ncol=ncol,
            prop={"size": legend_fontsize},
        )
        if leg.get_title() is not None:
            leg.get_title().set_fontsize(legend_title_fontsize)

    # ---- labels, ticks, colorbar ----
    ax.set_xlabel(x_label, fontsize=axis_label_fontsize)
    ax.set_ylabel(y_label, fontsize=axis_label_fontsize)

    # main column labels (rotated 45°)
    ax.tick_params(axis="x", which="both", bottom=True, top=False, labelbottom=True, pad=col_tick_pad)
    ax.set_xticklabels(
        ax.get_xticklabels(),
        rotation=45,
        ha="right",
        rotation_mode="anchor",
        fontsize=col_tick_fontsize,
    )

    # increment labels (not rotated) ABOVE the rotated column labels
    if show_col_increments:
        ax_inc = ax.twiny()
        ax_inc.set_xlim(ax.get_xlim())
        ax_inc.set_xticks(ax.get_xticks())

        inc_labels = [str(int(col_increment_start) + i) for i in range(n_cols)]
        ax_inc.set_xticklabels(
            inc_labels,
            rotation=0,
            ha="center",
            fontsize=(col_tick_fontsize if increment_tick_fontsize is None else increment_tick_fontsize),
        )

        ax_inc.xaxis.set_ticks_position("bottom")
        ax_inc.xaxis.set_label_position("bottom")
        ax_inc.spines["top"].set_visible(False)
        ax_inc.spines["bottom"].set_visible(False)
        ax_inc.tick_params(axis="x", which="both", length=0, pad=increment_tick_pad)
        ax_inc.set_xlabel("")
        ax_inc.set_frame_on(False)

    ax.set_yticklabels(
        ax.get_yticklabels(),
        rotation=0,
        fontsize=row_tick_fontsize,
    )

    if cg.cax is not None:
        cg.cax.set_ylabel(cbar_label, fontsize=cbar_label_fontsize)
        cg.cax.tick_params(labelsize=cbar_tick_fontsize)

    # ---- saving ----
    if savepath is not None:
        if rasterize_heatmap:
            for coll in getattr(cg.ax_heatmap, "collections", []):
                try:
                    coll.set_rasterized(True)
                except Exception:
                    pass

        fig = getattr(cg, "figure", None) or getattr(cg, "fig", None) or plt.gcf()
        fig.savefig(
            str(savepath),
            dpi=dpi,
            bbox_inches=bbox_inches,
            pad_inches=pad_inches,
            transparent=transparent,
        )
        if close:
            plt.close(fig)
            return cg

    plt.show()
    return cg






def _split_into_n_chunks(n: int, num_chunks: int) -> List[List[int]]:
    """Return `num_chunks` ~equal index chunks covering range(n)."""
    num_chunks = max(1, min(num_chunks, n))
    base, rem = divmod(n, num_chunks)
    sizes = [base + (1 if i < rem else 0) for i in range(num_chunks)]
    chunks, start = [], 0
    for sz in sizes:
        if sz > 0:
            chunks.append(list(range(start, start + sz)))
        start += sz
    return chunks

def _worker_chunk_process(
    chunk: List[Tuple[int, str]],
    virus_name: str,
    fetch_vals_fn: Callable[..., Tuple[np.ndarray, np.ndarray]],
    compute_fn: Callable[..., Mapping[str, object]],
    fetch_kwargs: dict | None,
    compute_kwargs: dict | None,
) -> List[Tuple[int, str, float]]:
    out: List[Tuple[int, str, float]] = []
    for idx, gene in chunk:
        pos_vals, neg_vals = fetch_vals_fn(gene_name=gene, virus_name=virus_name, **(fetch_kwargs or {}))
        res = compute_fn(neg_vals, pos_vals, **(compute_kwargs or {}))
        out.append((idx, gene, float(res["point"])))
    return out

def compute_stat_all_genes(
    list_a: Iterable,  # genes
    virus_name: str,
    fetch_vals_fn: Callable[..., Tuple[np.ndarray, np.ndarray]],
    compute_fn: Callable[..., Mapping[str, object]],
    *,
    fetch_kwargs: dict | None = None,
    compute_kwargs: dict | None = None,
    progress: bool = True,
    n_workers: int | None = None,       # max concurrent chunks
    chunk_size: int | None = None,      # desired genes per chunk (overrides equal-split)
    backend: str = "threads",           # 'threads' or 'processes'
) -> pd.DataFrame:
    """
    Parallel compute of point estimates for each gene vs `virus_name`.

    Chunking:
      - If `chunk_size` is provided: number of chunks = ceil(N / chunk_size).
      - Else: one chunk per worker (≈ equal parts).
    Concurrency is limited by `n_workers` (defaults to CPU count).
    Reports once per finished chunk; no tqdm.
    """
    fetch_kwargs = fetch_kwargs or {}
    compute_kwargs = compute_kwargs or {}

    genes = list(list_a)
    n = len(genes)
    if n == 0:
        return pd.DataFrame(columns=["gene", "point"])

    if backend not in {"threads", "processes"}:
        raise ValueError("backend must be 'threads' or 'processes'")

    desired_workers = n_workers or (os.cpu_count() or 1)

    # Decide number of chunks
    if chunk_size is not None:
        if chunk_size <= 0:
            raise ValueError("chunk_size must be a positive integer")
        num_chunks = max(1, min(math.ceil(n / chunk_size), n))
    else:
        num_chunks = max(1, min(desired_workers, n))  # equal-ish split

    idx_chunks = _split_into_n_chunks(n, num_chunks)
    chunks: List[List[Tuple[int, str]]] = [[(i, genes[i]) for i in idxs] for idxs in idx_chunks]

    # Limit concurrency to available workers but allow more chunks than workers
    max_workers = max(1, min(desired_workers, num_chunks))
    print(f"max_workers: {max_workers}")

    results: List[Tuple[int, str, float]] = []
    total_chunks = len(chunks)
    processed_chunks = 0
    processed_genes = 0
    t0 = time.time()

    if backend == "threads":
        def _worker_chunk_thread(chunk_local: List[Tuple[int, str]]) -> List[Tuple[int, str, float]]:
            out: List[Tuple[int, str, float]] = []
            for idx, gene in chunk_local:
                pos_vals, neg_vals = fetch_vals_fn(gene_name=gene, virus_name=virus_name, **fetch_kwargs)
                res = compute_fn(neg_vals, pos_vals, **compute_kwargs)
                out.append((idx, gene, float(res["point"])))
            return out

        with ThreadPoolExecutor(max_workers=max_workers) as ex:
            futs = [ex.submit(_worker_chunk_thread, ch) for ch in chunks]
            for fut in as_completed(futs):
                chunk_res = fut.result()
                results.extend(chunk_res)
                processed_chunks += 1
                processed_genes += len(chunk_res)
                if progress:
                    elapsed = time.time() - t0
                    pct = 100.0 * processed_genes / n
                    print(f"[{processed_chunks}/{total_chunks}] chunk done "
                          f"({len(chunk_res)} genes). {processed_genes}/{n} genes "
                          f"({pct:.1f}%) — {elapsed:.1f}s elapsed.", flush=True)
    else:
        with ProcessPoolExecutor(max_workers=max_workers) as ex:
            futs = [
                ex.submit(_worker_chunk_process, ch, virus_name, fetch_vals_fn, compute_fn, fetch_kwargs, compute_kwargs)
                for ch in chunks
            ]
            for fut in as_completed(futs):
                chunk_res = fut.result()
                results.extend(chunk_res)
                processed_chunks += 1
                processed_genes += len(chunk_res)
                if progress:
                    elapsed = time.time() - t0
                    pct = 100.0 * processed_genes / n
                    print(f"[{processed_chunks}/{total_chunks}] chunk done "
                          f"({len(chunk_res)} genes). {processed_genes}/{n} genes "
                          f"({pct:.1f}%) — {elapsed:.1f}s elapsed.", flush=True)

    # Restore original order
    results.sort(key=lambda t: t[0])
    _, ordered_genes, ordered_points = zip(*results)
    return pd.DataFrame({"gene": ordered_genes, "point": ordered_points})



def get_expr_vals(gene_name, virus_name, counts_table, species2runs):
    """
    Get the expr values of a specific gene across virus-positive and virus-negative runs.
    counts_table: pandas DataFrame, columns=genes, rows=samples.
    gene_name: str, id of the gene to plot.
    virus_name: str, virus species name
    species2runs: dict, species to a dict of positive runs
    
    """
    gene_series = counts_table[gene_name]
    
    pos_runs = [i.removesuffix("_salmon") for i in species2runs[virus_name]]
    pos_runs_series= gene_series.reindex(pos_runs)
    pos_runs_val = pos_runs_series.dropna().to_list()
    
    neg_runs_val = gene_series.drop(index = pos_runs, errors = "ignore").tolist()

    return pos_runs_val, neg_runs_val



def compare_distributions(a, b, test='ks'):
    """
    Compare two 1D samples with a selected statistical two‐sample test.

    Parameters
    ----------
    a : array‐like of float
        First sample.
    b : array‐like of float
        Second sample.
    test : {'ks', 'ad', 'cvm', 'ttest', 'mannwhitney', 'all'}
        Which test to run:
        - 'ks'          : Kolmogorov–Smirnov two‐sample test [sensitive in both location and shape]
        - 'ad'          : Anderson–Darling k‐sample test [more sensitive to deviations in the tails]
        - 'cvm'         : Cramér–von Mises two‐sample test [intermediate sensitivity]
        - 'ttest'       : Student's T-test [parametric test sensitive to mean differences]
        - 'mannwhitney' : Mann-Whitney U test [nonparametric test sensitive to median differences]
        - 'all'         : run all tests and return a dict

    Returns
    -------
    result : dict
        Returns a dict containing statistics and p‐values for each test.
    """

    a = np.asarray(a)
    b = np.asarray(b)

    def run_ks():
        d, p = stats.ks_2samp(a, b)
        return {'statistic': d, 'pvalue': p}

    def run_ad():
        ad = stats.anderson_ksamp([a, b])
        p = ad.pvalue if hasattr(ad, 'pvalue') else ad.significance_level / 100.0
        return {'statistic': ad.statistic, 'pvalue': p}

    def run_cvm():
        cvm = stats.cramervonmises_2samp(a, b)
        return {'statistic': cvm.statistic, 'pvalue': cvm.pvalue}

    def run_ttest():
        t_stat, p = stats.ttest_ind(a, b, equal_var=False)
        return {'statistic': t_stat, 'pvalue': p}

    def run_mannwhitney():
        u_stat, p = stats.mannwhitneyu(a, b, alternative='two-sided', method= 'auto')
        f = u_stat / (len(a) * len(b))
        return {'statistic': f, 'pvalue': p}

    test = test.lower()
    tests = {
        'ks': run_ks,
        'ad': run_ad,
        'cvm': run_cvm,
        'ttest': run_ttest,
        'mannwhitney': run_mannwhitney
    }

    if test == 'all':
        return {name: func() for name, func in tests.items()}
    elif test in tests:
        return tests[test]()
    else:
        raise ValueError(f"Unknown test '{test}'. Choose from 'ks', 'ad', 'cvm', 'ttest', 'mannwhitney', or 'all'.")


def plot_distribution(
    data,
    bins=50,
    title="Distribution with Markers",
    hist_kwargs=None,
    line_kwargs=None,
    arrow_dx=15,
    arrow_dy=10,
):
    data = np.asarray(data).ravel()

    hist_kwargs = {} if hist_kwargs is None else dict(hist_kwargs)
    line_kwargs = {"linestyle": "--", "linewidth": 1.5} if line_kwargs is None else dict(line_kwargs)
    hist_kwargs.setdefault("alpha", 0.7)  # default alpha only if not provided

    fig, ax = plt.subplots(figsize=(12, 6))
    ax.hist(data, bins=bins, density=True, **hist_kwargs)

    ymin, ymax = ax.get_ylim()

    ax.set_title(title)
    ax.set_xlabel("Value")
    ax.set_ylabel("Density")
    ax.margins(x=0.02)
    plt.show()
    return fig, ax



## outlier detection
# ---------- core helpers ----------

def _coerce_numeric(series: pd.Series) -> pd.Series:
    vals = pd.to_numeric(series, errors="coerce").dropna()
    if len(vals) == 0:
        raise ValueError("No numeric data after coercion/dropna.")
    return vals

def _mad_star(values: pd.Series) -> float:
    """
    Scaled MAD (consistent for sigma at normal): 1.4826 * median(|x - median(x)|)
    """
    med = values.median()
    mad = np.median(np.abs(values - med))
    return 1.4826 * mad

def _medcouple(values: pd.Series, require_mc: bool = False) -> float:
    """
    Robust skewness (Hubert-Vandervieren adjusted boxplot).
    We try to import a fast implementation; otherwise fall back to 0 (no skew),
    which reduces to standard Tukey fences.
    """
    try:
        # pip install medcouple (https://pypi.org/project/medcouple/)
        from medcouple import medcouple as mc  # type: ignore
        return float(mc(np.asarray(values, dtype=float)))
    except Exception as e:
        if require_mc:
            raise RuntimeError(
                "Medcouple ('medcouple' package) is required for method='adjusted_boxplot'."
            ) from e
        warnings.warn("medcouple package not found; using MC=0 (standard Tukey).")
        return 0.0

def _gpd_upper_cutoff(values: pd.Series, base_q: float = 0.95, tail_alpha: float = 0.005) -> float:
    """
    EVT/POT upper cutoff via GPD on exceedances above u=quantile(base_q).
    Returns the cutoff with tail probability tail_alpha relative to the whole sample.
    Requires SciPy.
    """
    try:
        from scipy.stats import genpareto
    except Exception as e:
        raise RuntimeError("SciPy is required for method='evd' (GPD).") from e

    x = np.asarray(values, dtype=float)
    u = np.quantile(x, base_q)
    exc = x[x > u] - u
    if len(exc) < 20:
        warnings.warn(
            f"Only {len(exc)} exceedances above u={base_q:.3f} quantile; EVT fit may be unstable."
        )
    # Fit GPD with location fixed at 0
    c, loc, scale = genpareto.fit(exc, floc=0.0)
    # Desired unconditional tail probability = tail_alpha.
    # Let p_u = P(X>u) = 1 - base_q. We need conditional exceedance prob p* s.t.
    # p* * p_u = tail_alpha  =>  p* = tail_alpha / p_u
    p_u = 1.0 - base_q
    p_star = min(max(tail_alpha / max(p_u, 1e-12), 1e-9), 1 - 1e-9)  # clamp
    # Inverse CDF of GPD for exceedances:
    from scipy.stats import genpareto
    q_exc = genpareto.ppf(1 - p_star, c, loc=0.0, scale=scale)
    return float(u + q_exc)

# ---------- public API ----------

Method = Literal["tukey", "mad", "quantile", "ipr", "adjusted_boxplot", "evd"]
Side = Literal["both", "upper", "lower"]

def compute_fences(
    series: pd.Series,
    method: Method = "tukey",
    *,
    side: Side = "both",
    # common knobs
    factor: float = 1.5,                 # Tukey multiplier
    k: float = 3.5,                      # MAD/IPR multiplier
    alpha: float = 0.01,                 # quantile tail mass for 'quantile'
    ipr_low: float = 0.10,               # lower quantile for IPR
    ipr_high: float = 0.90,              # upper quantile for IPR
    require_mc: bool = False,            # for adjusted_boxplot
    # EVT knobs
    base_q: float = 0.95,                # POT threshold quantile
    tail_alpha: float = 0.005            # desired unconditional tail prob for EVT cutoff
) -> Tuple[Optional[float], Optional[float]]:
    """
    Compute outlier fences without assuming normality.

    Returns (lower, upper). If a side does not apply in a method (e.g., EVT is upper-only),
    the non-applicable side will be None.

    Methods:
    - 'tukey'             : Q1 -/+ factor*IQR (factor=1.5 default)
    - 'mad'               : median -/+ k*MAD* (k=3.5; MAD* = 1.4826*MAD)
    - 'quantile'          : (Q_alpha, Q_{1-alpha})
    - 'ipr'               : median -/+ k*(Q_{ipr_high}-Q_{ipr_low})
    - 'adjusted_boxplot'  : Hubert–Vandervieren skew-adjusted Tukey using medcouple
    - 'evd'               : EVT/POT (GPD) for the upper fence only; lower=None

    Parameters can be tuned via kwargs (see signature).
    """
    values = _coerce_numeric(series)

    if method == "tukey":
        q1, q3 = values.quantile([0.25, 0.75])
        iqr = q3 - q1
        L = float(q1 - factor * iqr)
        U = float(q3 + factor * iqr)

    elif method == "mad":
        med = float(values.median())
        mad_s = _mad_star(values)
        L = med - k * mad_s
        U = med + k * mad_s

    elif method == "quantile":
        if not (0 < alpha < 0.5):
            raise ValueError("alpha must be in (0, 0.5) for method='quantile'.")
        L = float(values.quantile(alpha))
        U = float(values.quantile(1 - alpha))

    elif method == "ipr":
        if not (0 < ipr_low < ipr_high < 1):
            raise ValueError("Require 0 < ipr_low < ipr_high < 1 for method='ipr'.")
        ql = float(values.quantile(ipr_low))
        qh = float(values.quantile(ipr_high))
        ipr = qh - ql
        med = float(values.median())
        L = med - k * ipr
        U = med + k * ipr

    elif method == "adjusted_boxplot":
        q1, q3 = values.quantile([0.25, 0.75])
        iqr = q3 - q1
        MC = _medcouple(values, require_mc=require_mc)
        if MC >= 0:
            L = float(q1 - factor * np.exp(-4 * MC) * iqr)
            U = float(q3 + factor * np.exp( 3 * MC) * iqr)
        else:
            L = float(q1 - factor * np.exp(-3 * MC) * iqr)
            U = float(q3 + factor * np.exp( 4 * MC) * iqr)

    elif method == "evd":
        # Upper-only threshold using POT/GPD
        U = _gpd_upper_cutoff(values, base_q=base_q, tail_alpha=tail_alpha)
        L = None

    else:
        raise ValueError(f"Unknown method: {method}")

    # Respect side argument by nulling out the opposite fence(s)
    if side == "upper":
        return None, U
    elif side == "lower":
        return L, None
    else:
        return L, U

def upper_whisker_threshold(series: pd.Series) -> float:
    """Original behavior: Tukey upper fence (Q3 + 1.5*IQR)."""
    _, U = compute_fences(series, method="tukey", side="upper", factor=1.5)
    assert U is not None
    return float(U)

def upper_threshold(
    series: pd.Series,
    method: Method = "tukey",
    **kwargs
) -> float:
    """
    Generalized upper threshold. Examples:
      upper_threshold(s, method='mad', k=3.5)
      upper_threshold(s, method='quantile', alpha=0.01)
      upper_threshold(s, method='ipr', ipr_low=0.10, ipr_high=0.90, k=0.74)
      upper_threshold(s, method='adjusted_boxplot')  # requires 'medcouple' pkg for skew adjust
      upper_threshold(s, method='evd', base_q=0.95, tail_alpha=0.005)  # EVT POT
    """
    _, U = compute_fences(series, method=method, side="upper", **kwargs)
    if U is None:
        raise ValueError(f"Method '{method}' does not define an upper threshold.")
    return float(U)

def lower_threshold(
    series: pd.Series,
    method: Method = "tukey",
    **kwargs
) -> float:
    """
    Generalized lower threshold. Examples:
    """
    L, _ = compute_fences(series, method=method, side="lower", **kwargs)
    if L is None:
        raise ValueError(f"Method '{method}' does not define a lower threshold.")
    return float(L)

def thresholds_for_column(
    col_name: str,
    df: pd.DataFrame,
    *,
    method: Method = "tukey",
    **kwargs
) -> Tuple[float, float]:
    """
    Backward-compatible column lookup with flexible method selection.
    Expects `all_shift[col_name]` to be a DataFrame with a 'point' column.
    """
    df_for_col = df[col_name]  # KeyError if missing (intentional)
    if "point" not in df_for_col.columns:
        raise KeyError(f"'{col_name}' lookup DataFrame is missing a 'point' column.")
    return upper_threshold(df_for_col["point"], method=method, **kwargs), lower_threshold(df_for_col["point"], method=method, **kwargs)


def parse_virus_categories(
    virus_list: list,
    category_keywords: list = ["endogenous_or_nonfish", "insufficient_evidence"],
    default_label: Optional[str] = None,
) -> Tuple[List[str], List[Optional[str]]]:
    """
    Parse a virus list and create a label per entry based on the presence of
    any of `category_keywords`.

    Parameters
    ----------
    virus_list : list
        List of virus strings (or values castable to str).
    category_keywords : list[str]
        Keywords to detect and strip.
    default_label : Optional[str]
        Label to use when no keyword is found for an entry.
        If None, the label will be None in that case.

    Returns
    -------
    cleaned_viruses : list[str]
        Virus strings with any matched keywords stripped out.
    labels : list[Optional[str]]
        Matched keywords per entry, joined by "+".
        If none matched, `default_label` is used.
    """
    if virus_list is None:
        return [], []

    delim = r"[\s\|,;/\[\]\(\)\{\}_-]"

    kw_patterns = {
        kw: re.compile(rf"(?i)(^|{delim}){re.escape(kw)}(?=$|{delim})")
        for kw in category_keywords
    }

    def _cleanup(s: str) -> str:
        s = s.strip()
        s = re.sub(r"[\[\]\(\)\{\}]", "", s)
        s = re.sub(r"\s+", " ", s)
        s = re.sub(r"_+", "_", s)
        s = re.sub(r"-+", "-", s)
        s = re.sub(r"\|+", "|", s)
        s = re.sub(r",+", ",", s)
        s = re.sub(r";+", ";", s)
        s = re.sub(r"/+", "/", s)
        s = s.strip(" _-|,;/")
        s = re.sub(r"\s*([_-])\s*", r"\1", s)
        return re.sub(r"\s+", " ", s).strip()

    cleaned_viruses: List[str] = []
    labels: List[Optional[str]] = []

    for item in virus_list:
        s = "" if item is None else str(item)
        found: List[str] = []

        for kw in category_keywords:
            pat = kw_patterns[kw]
            if pat.search(s):
                found.append(kw)
                s = pat.sub("", s)

        cleaned_viruses.append(_cleanup(s))
        labels.append("+".join(found) if found else default_label)

    return cleaned_viruses, labels