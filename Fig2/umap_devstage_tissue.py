#!/usr/bin/env python3
"""
Figure 1 Panel C: Transcriptome clustering UMAPs

Creates UMAP plots of the 75k zebrafish transcriptome (61k samples with
expression data), colored by:
  1. devstage_curation (6 categories: Embryo, Larval, Juvenile, Adult,
     Multi-stage, Undetermined)
  2. tissue_curation_coarse (21 anatomical categories)

Legends (with proportional stacked bars) are drawn as separate figures.

Data source: pre-processed anndata with log2 TMM-CPM normalization.
Pre-computed PCA/UMAP embeddings are used if present; otherwise computed
following the reference UMAP-embedding notebook (internal).
"""

import argparse
import numpy as np
import scipy.sparse
import polars as pl
import anndata as ad
import scanpy as sc
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.patches as mpatches
import re
from pathlib import Path

# =============================================================================
# Configuration
# =============================================================================

_parser = argparse.ArgumentParser(description=__doc__)
_parser.add_argument(
    "--anndata-h5ad", type=Path,
    default=Path("data/75k_unstable/75k_unstable_anndata_zfin_aliases_metadata.log2tmmcpm.h5ad"),
    help="Transcriptome anndata (.h5ad) with pre-computed PCA/UMAP embeddings.",
)
_parser.add_argument(
    "--metadata-file", type=Path,
    default=Path("data/metadata/all_zf_dates_devstage_tissue_tech_curated.tsv"),
    help="Curated SRA/GEO metadata TSV (devstage, tissue, tech).",
)
_args, _ = _parser.parse_known_args()

ANNDATA_H5AD = _args.anndata_h5ad
METADATA_FILE = _args.metadata_file
BATLOW_PALETTE = Path("palette/batlow/DiscretePalettes/batlow100.txt")
OUTPUT_DIR = Path("figures")
OUTPUT_DIR.mkdir(exist_ok=True)

CM = 1 / 2.54
FIG_W = 8.24 * CM   # umap panel: 8.24 × 6.25 cm
FIG_H = 6.25 * CM
LEGEND_W = 8.24 * CM  # legend: 8.24 × 0.89 cm (but allow overflow vertically)
LEGEND_H = 0.89 * CM * (7.0 / 3.0)

# "dot"  → circle + label text flowing horizontal-first (new style)
# "bar"  → proportional bar per category (old style)
LEGEND_STYLE = "dot"

# UMAP parameters (used only if UMAP not pre-computed)
N_PCS = 50
N_NEIGHBORS = 100
METRIC = "euclidean"
MIN_DIST = 0.5
SPREAD = 1.0
HVG_N_TOP = 2000

# Font settings (Arial, 5–7pt per layout spec)
mpl.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['font.sans-serif'] = ['Arial', 'Helvetica Neue', 'DejaVu Sans']
mpl.rcParams['svg.fonttype'] = 'none'

FONT_TITLE = 6
FONT_LABEL = 6
FONT_LEGEND = 5
FONT_LEGEND_DOT = 4  # dot-style legend labels
FONT_TICK = 6

# Devstage chronological order and colors (RdYlBu-inspired gradient)
DEVSTAGE_ORDER = [
    "Zygote", "Cleavage", "Blastula", "Gastrula", "Segmentation",
    "Pharyngula", "Hatching", "Larval", "Juvenile", "Adult",
    "Multi-stage", "Undetermined",
]
DEVSTAGE_COLORS = [
    "#a50026",  # Zygote        - dark red
    "#d73027",  # Cleavage      - red
    "#f46d43",  # Blastula      - orange-red
    "#fdae61",  # Gastrula      - orange
    "#fee090",  # Segmentation  - yellow-orange
    "#ffffbf",  # Pharyngula    - light yellow
    "#e0f3f8",  # Hatching      - light cyan
    "#abd9e9",  # Larval        - light blue
    "#74add1",  # Juvenile      - medium blue
    "#4575b4",  # Adult         - blue
    "#313695",  # Multi-stage   - dark blue
    "#808080",  # Undetermined  - gray
]

# Gray overrides for uninformative categories
GRAY_CATS = {
    "Undetermined": "#9E9E9E",
    "Multi-tissue/Undetermined": "#BDBDBD",
}

# Hand-picked tissue palette (Zebrahub-inspired hue families)
TISSUE_COLORS = {
    "Nervous System":            "darkcyan",
    "Sensory System":            "lightseagreen",
    "Muscular System":           "mediumpurple",
    "Skeletal Element":          "rosybrown",
    "Cardiovascular System":     "coral",
    "Hematopoietic System":      "crimson",
    "Digestive System":          "olivedrab",
    "Liver and Biliary System":  "goldenrod",
    "Endocrine System":          "orchid",
    "Renal System":              "steelblue",
    "Respiratory System":        "cornflowerblue",
    "Reproductive System":       "palevioletred",
    "Surface Structure":         "burlywood",
    "Adipose Tissue":            "palegoldenrod",
    "Swim Bladder":              "mediumaquamarine",
    "Cell Line":                 "darkgray",
    "Cancer or Tumor":           "saddlebrown",
    "Multi-tissue/Undetermined": "#BDBDBD",
}

# =============================================================================
# Palette loading
# =============================================================================

def load_batlow_palette():
    """Load all colors from batlow100 discrete palette, return as list of hex."""
    colors = []
    if BATLOW_PALETTE.exists():
        for line in BATLOW_PALETTE.read_text().strip().split('\n'):
            if not line.strip() or line.startswith('S'):
                continue
            parts = line.split()
            hex_parts = [p for p in parts if p.startswith('#')]
            if hex_parts:
                colors.append(hex_parts[0])
    if not colors:
        colors = [mcolors.to_hex(c) for c in plt.get_cmap("tab20").colors]
    return colors


def load_batlow_colors_stride(n_colors):
    """
    Sample n_colors from batlow palette using greedy maxmin (furthest-point) sampling.
    Maximizes minimum distance between selected colors for perceptual distinctness.
    """
    colors = load_batlow_palette()
    if n_colors >= len(colors):
        return colors

    def hex_to_rgb(hex_color):
        """Convert hex color to RGB tuple (0-1)."""
        r = int(hex_color[1:3], 16) / 255.0
        g = int(hex_color[3:5], 16) / 255.0
        b = int(hex_color[5:7], 16) / 255.0
        return (r, g, b)

    def rgb_distance(rgb1, rgb2):
        """Euclidean distance in RGB space."""
        return sum((a - b) ** 2 for a, b in zip(rgb1, rgb2)) ** 0.5

    color_rgbs = [hex_to_rgb(c) for c in colors]
    selected_indices = [0]
    selected_rgbs = [color_rgbs[0]]

    for _ in range(n_colors - 1):
        best_idx = None
        best_min_dist = -1

        for i, rgb in enumerate(color_rgbs):
            if i in selected_indices:
                continue
            min_dist = min(rgb_distance(rgb, selected_rgb) for selected_rgb in selected_rgbs)
            if min_dist > best_min_dist:
                best_min_dist = min_dist
                best_idx = i

        selected_indices.append(best_idx)
        selected_rgbs.append(color_rgbs[best_idx])

    return [colors[i] for i in selected_indices]


def build_devstage_palette(cats_by_count):
    """
    Devstage palette: fixed gradient colors in chronological order.
    Returns dict: category → hex color.
    """
    return dict(zip(DEVSTAGE_ORDER, DEVSTAGE_COLORS))


def build_tissue_palette(cats_by_count):
    """
    Tissue palette: hand-picked Zebrahub-inspired colors, with batlow fallback
    for any category not in TISSUE_COLORS.
    """
    missing = [c for c in cats_by_count if c not in TISSUE_COLORS]
    if missing:
        fallback = load_batlow_colors_stride(len(missing))
        color_map = {c: fallback[i] for i, c in enumerate(missing)}
    else:
        color_map = {}
    color_map.update({c: mcolors.to_hex(mcolors.CSS4_COLORS.get(v, v))
                      if not str(v).startswith("#") else v
                      for c, v in TISSUE_COLORS.items()})
    return color_map


def remap_tissue_categories(adata):
    """
    Merge uninformative tissue categories into 'Multi-tissue/Undetermined'.
    """
    cats_to_merge = {
        "All anatomical structures",
        "Embryo Imprecise",
        "Multi-system",
        "Undetermined",
    }
    key = "tissue_curation_coarse"
    adata.obs[key] = adata.obs[key].astype(str)
    mask = adata.obs[key].isin(cats_to_merge)
    adata.obs.loc[mask, key] = "Multi-tissue/Undetermined"
    adata.obs[key] = adata.obs[key].astype("category")


# =============================================================================
# Load anndata
# =============================================================================

print(f"Loading anndata from {ANNDATA_H5AD} ...")
adata = ad.read_h5ad(ANNDATA_H5AD)
print(f"  {adata.n_obs:,} samples × {adata.n_vars:,} genes")
print(f"  X_umap pre-computed: {'X_umap' in adata.obsm}")

# Remap tissue categories
remap_tissue_categories(adata)

# =============================================================================
# Compute UMAP (only if not pre-computed)
# =============================================================================

if "X_umap" not in adata.obsm:
    print("Computing UMAP from scratch (not found in obsm)...")

    X = adata.X
    if scipy.sparse.issparse(X):
        gene_var = np.asarray(X.power(2).mean(axis=0) - np.square(X.mean(axis=0))).flatten()
    else:
        gene_var = np.var(X, axis=0)

    adata = adata[:, gene_var > 1e-3].copy()
    print(f"  Kept {adata.n_vars:,} genes after low-variance filter")

    X = adata.X
    if scipy.sparse.issparse(X):
        gene_var = np.asarray(X.power(2).mean(axis=0) - np.square(X.mean(axis=0))).flatten()
    else:
        gene_var = np.var(X, axis=0)

    hvg_mask = np.zeros(adata.n_vars, dtype=bool)
    hvg_mask[np.argsort(gene_var)[::-1][:HVG_N_TOP]] = True
    adata.var["highly_variable"] = hvg_mask
    adata = adata[:, adata.var["highly_variable"]].copy()
    print(f"  Using {adata.n_vars} highly variable genes")

    sc.pp.scale(adata, zero_center=True, max_value=10)
    sc.tl.pca(adata, n_comps=N_PCS, svd_solver="arpack")
    sc.pp.neighbors(adata, n_neighbors=N_NEIGHBORS, n_pcs=N_PCS, metric=METRIC)
    sc.tl.umap(adata, min_dist=MIN_DIST, spread=SPREAD)
    print(f"  UMAP shape: {adata.obsm['X_umap'].shape}")

# =============================================================================
# Load metadata category counts (for color rank and legend proportions)
# =============================================================================

_meta = pl.read_csv(
    METADATA_FILE,
    separator="\t",
    columns=["devstage_curation", "tissue_curation_coarse"],
)

# Remap tissue categories in metadata
cats_to_merge = {
    "All anatomical structures",
    "Embryo Imprecise",
    "Multi-system",
    "Undetermined",
}
_meta = _meta.with_columns(
    pl.when(pl.col("tissue_curation_coarse").is_in(cats_to_merge))
    .then(pl.lit("Multi-tissue/Undetermined"))
    .otherwise(pl.col("tissue_curation_coarse"))
    .alias("tissue_curation_coarse")
)

def _count_ordered(col):
    return (
        _meta.filter(pl.col(col).is_not_null())
        .group_by(col)
        .agg(pl.len().alias("n"))
        .sort("n", descending=True)[col]
        .to_list()
    )

devstage_counts_map = (
    _meta.filter(pl.col("devstage_curation").is_not_null())
    .group_by("devstage_curation")
    .agg(pl.len().alias("n"))
    .to_pandas()
    .set_index("devstage_curation")["n"]
    .to_dict()
)
tissue_counts_map = (
    _meta.filter(pl.col("tissue_curation_coarse").is_not_null())
    .group_by("tissue_curation_coarse")
    .agg(pl.len().alias("n"))
    .to_pandas()
    .set_index("tissue_curation_coarse")["n"]
    .to_dict()
)

# =============================================================================
# Build palettes and register in anndata
# =============================================================================

key_devstage = "devstage_curation"
adata.obs[key_devstage] = adata.obs[key_devstage].astype("category")
devstage_cats = list(adata.obs[key_devstage].cat.categories)  # alphabetical
print(f"\nDevstage categories ({len(devstage_cats)}): {devstage_cats}")

devstage_color_map = build_devstage_palette(_count_ordered("devstage_curation"))
adata.uns[key_devstage + "_colors"] = [devstage_color_map.get(c, "#9E9E9E") for c in devstage_cats]

key_tissue = "tissue_curation_coarse"
adata.obs[key_tissue] = adata.obs[key_tissue].astype("category")
tissue_color_map = build_tissue_palette(_count_ordered("tissue_curation_coarse"))
tissue_cats = list(adata.obs[key_tissue].value_counts().index)  # by count desc
adata.obs[key_tissue] = adata.obs[key_tissue].cat.reorder_categories(tissue_cats)
print(f"Tissue categories ({len(tissue_cats)}): {tissue_cats}")
adata.uns[key_tissue + "_colors"] = [tissue_color_map.get(c, "#9E9E9E") for c in tissue_cats]

# =============================================================================
# UMAP figure helpers
# =============================================================================

mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"] = 42
mpl.rcParams["svg.fonttype"] = "none"
mpl.rcParams["figure.figsize"] = [FIG_W, FIG_H]


def _add_umap_axis_inset(ax):
    """
    Add a small bottom-left inset showing UMAP axis orientation via arrows.
    Drawn as annotation arrows + text, no real sub-axes.
    """
    # Position in axes-fraction coords: bottom-left corner
    arrow_len = 0.15   # fraction of axes width/height
    ox, oy = 0.04, 0.04  # origin

    # Up arrow: UMAP 2
    ax.annotate(
        "", xy=(ox, oy + arrow_len), xytext=(ox, oy),
        xycoords="axes fraction", textcoords="axes fraction",
        arrowprops=dict(arrowstyle="-|>", color="black", lw=0.6,
                        mutation_scale=4),
    )
    ax.text(ox + 0.012, oy + arrow_len / 2, "UMAP 2",
            transform=ax.transAxes,
            fontsize=FONT_TICK, va="center", ha="left", color="black",
            rotation=90)

    # Right arrow: UMAP 1
    ax.annotate(
        "", xy=(ox + arrow_len, oy), xytext=(ox, oy),
        xycoords="axes fraction", textcoords="axes fraction",
        arrowprops=dict(arrowstyle="-|>", color="black", lw=0.6,
                        mutation_scale=4),
    )
    ax.text(ox + arrow_len / 2, oy - 0.025, "UMAP 1",
            transform=ax.transAxes,
            fontsize=FONT_TICK, va="top", ha="center", color="black")


def _apply_umap_style(fig):
    """Remove borders/ticks, fill figure bounds, add axis inset to each axes."""
    for ax in fig.axes:
        ax.set_xlabel("")
        ax.set_ylabel("")
        ax.tick_params(left=False, bottom=False, labelleft=False, labelbottom=False)
        for spine in ax.spines.values():
            spine.set_visible(False)
        ax.title.set_fontsize(FONT_TITLE)
        _add_umap_axis_inset(ax)
    fig.tight_layout(pad=0)


# =============================================================================
# UMAP plots
# =============================================================================

print("\nPlotting UMAP: developmental stage...")

adata_devstage = adata[adata.obs[key_devstage].sort_values(ascending=True).index]
fig = sc.pl.umap(
    adata_devstage,
    color=key_devstage,
    title="",
    legend_loc=None,
    add_outline=True,
    outline_width=(0.05, 0.05),
    outline_color=("black", "white"),
    size=3,
    alpha=0.8,
    sort_order=False,
    show=False,
    return_fig=True,
)
_apply_umap_style(fig)
out_devstage = OUTPUT_DIR / "Fig1_D_1_umap_devstage.svg"
fig.savefig(out_devstage, transparent=True)
print(f"Saved {out_devstage}")
plt.close()

print("Plotting UMAP: tissue type...")

adata_tissue = adata[adata.obs[key_tissue].sort_values(ascending=True).index]
fig = sc.pl.umap(
    adata_tissue,
    color=key_tissue,
    title="",
    legend_loc=None,
    add_outline=True,
    outline_width=(0.05, 0.05),
    outline_color=("black", "white"),
    size=3,
    alpha=0.8,
    sort_order=False,
    show=False,
    return_fig=True,
)
_apply_umap_style(fig)
out_tissue = OUTPUT_DIR / "Fig1_D_1_umap_tissue.svg"
fig.savefig(out_tissue, transparent=True)
print(f"Saved {out_tissue}")
plt.close()

# =============================================================================
# Legend figures (color swatch + category label + proportional stacked bar)
# =============================================================================

def _svg_add_text_stroke(svg_path, stroke_color="white", stroke_width=1.15):
    """
    Post-process a saved SVG to give all <text> elements a native SVG stroke
    (paint-order: stroke fill) instead of matplotlib path effects.
    Keeps text as real <text> nodes rather than outlines.
    """
    content = svg_path.read_text(encoding="utf-8")

    def patch_style(m):
        style = m.group(1)
        style = re.sub(r'\bstroke:\s*none\b', f'stroke: {stroke_color}', style)
        if re.search(r'\bstroke-width:', style):
            style = re.sub(r'\bstroke-width:\s*[\d.]+\S*', f'stroke-width: {stroke_width}', style)
        else:
            style += f'; stroke-width: {stroke_width}'
        if 'paint-order' not in style:
            style += '; paint-order: stroke fill'
        if 'stroke:' not in style:
            style += f'; stroke: {stroke_color}'
        return f'style="{style}"'

    def patch_text_el(m):
        return re.sub(r'\bstyle="([^"]*)"', patch_style, m.group(0))

    content = re.sub(r'<text\b[^>]*>.*?</text>', patch_text_el, content, flags=re.DOTALL)
    svg_path.write_text(content, encoding="utf-8")


def _draw_legend_with_bar(color_map, ordered_cats, counts_map, output_path, title=None):
    """
    Draw a legend panel where each row is a colored background bar whose width
    is proportional to that category's share of total samples. The category label
    is drawn in black with a white outline, directly on top of the bar. No
    separate color swatch or stacked-bar column.
    """
    total = sum(counts_map.get(c, 0) for c in ordered_cats)
    counts = [counts_map.get(c, 0) for c in ordered_cats]
    colors = [color_map.get(c, "#9E9E9E") for c in ordered_cats]
    fracs = [c / total if total > 0 else 0 for c in counts]

    fig, ax = plt.subplots(figsize=(LEGEND_W, LEGEND_H))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    n_cats = len(ordered_cats)
    margin_top = 0.04
    margin_bottom = 0.02
    usable = 1.0 - margin_top - margin_bottom
    row_h = usable / n_cats
    gap = row_h * 0.15     # vertical gap between rows
    bar_h = row_h - gap
    text_pad = 0.015       # left padding for label text (axes fraction)

    for i, (cat, color, frac) in enumerate(zip(ordered_cats, colors, fracs)):
        y = 1.0 - margin_top - (i + 1) * row_h + gap / 2
        # Background bar: width = fraction of total
        ax.add_patch(mpatches.Rectangle(
            (0, y), frac, bar_h,
            facecolor=color, edgecolor="none",
            transform=ax.transAxes, clip_on=True,
        ))
        # Format count label: 1 decimal place in thousands if >= 1000, else integer
        n = counts_map.get(cat, 0)
        if n >= 1000:
            n_str = f"(n={n/1000:.1f}k)"
        else:
            n_str = f"(n={n})"
        label = f"{cat} {n_str}"
        # Label: black text — white outline applied via SVG stroke post-processing
        ax.text(
            text_pad, y + bar_h / 2, label,
            transform=ax.transAxes,
            fontsize=FONT_LEGEND, va="center", ha="left",
            color="black",
        )

    if title:
        fig.suptitle(title, fontsize=FONT_TITLE, y=1.01)

    fig.tight_layout(pad=0.2)
    fig.savefig(output_path, transparent=True)
    plt.close()
    if Path(output_path).suffix == ".svg":
        _svg_add_text_stroke(Path(output_path))
    print(f"Saved {output_path}")


def _draw_legend_dot(color_map, ordered_cats, output_path, counts_map=None):
    """
    Draw a legend as colored circles + labels flowing left→right, wrapping to
    a new row when the next item would exceed the figure width.
    Canvas is LEGEND_W × LEGEND_H; vertical overflow is clipped in post.
    If counts_map is provided, appends " (Xk)" sample counts to each label.
    """
    fig, ax = plt.subplots(figsize=(LEGEND_W, LEGEND_H))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    fig_w_pt = LEGEND_W * 72

    # Dot geometry (points)
    r_pt = FONT_LEGEND_DOT * 0.55
    r_x = r_pt / fig_w_pt
    r_y = r_pt / (LEGEND_H * 72)

    char_w_pt  = 0.50 * FONT_LEGEND_DOT
    text_gap_pt = r_pt * 0.8        # dot-to-text gap
    item_gap_pt = r_pt * 2.0        # trailing space after each item

    # Pre-scan: build all labels, find uniform item width from longest
    labels = []
    for cat in ordered_cats:
        if counts_map:
            nv = counts_map.get(cat, 0)
            n_str = f" ({nv/1000:.1f}k)" if nv >= 1000 else f" ({nv})"
        else:
            n_str = ""
        labels.append(cat + n_str)

    max_chars   = max(len(l) for l in labels) if labels else 1
    item_w_pt   = r_pt * 2 + text_gap_pt + char_w_pt * max_chars + item_gap_pt

    row_h   = r_y * 2.65
    margin_x_pt = r_pt * 1.5

    row = 0
    x_pt = margin_x_pt

    for cat, label in zip(ordered_cats, labels):
        # Wrap when next item would overflow (never wrap before first item in row)
        if x_pt + item_w_pt > fig_w_pt and x_pt > margin_x_pt:
            row += 1
            x_pt = margin_x_pt

        cx = (x_pt + r_pt) / fig_w_pt
        cy = 1.0 - (row + 0.5) * row_h

        circ = mpatches.Ellipse(
            (cx, cy), width=r_x * 2, height=r_y * 2,
            facecolor=color_map.get(cat, "#9E9E9E"), edgecolor="none",
            transform=ax.transAxes, clip_on=False,
        )
        ax.add_patch(circ)

        ax.text(
            (x_pt + r_pt * 2 + text_gap_pt) / fig_w_pt, cy, label,
            transform=ax.transAxes,
            fontsize=FONT_LEGEND_DOT, va="center", ha="left",
            color="black", clip_on=False,
        )

        x_pt += item_w_pt

    fig.tight_layout(pad=0)
    fig.savefig(output_path, transparent=True)
    plt.close()
    print(f"Saved {output_path}")


def _is_dark(hex_color):
    r, g, b = int(hex_color[1:3], 16), int(hex_color[3:5], 16), int(hex_color[5:7], 16)
    return 0.299 * r + 0.587 * g + 0.114 * b < 128


print("\nDrawing legends...")

tissue_ordered = _count_ordered("tissue_curation_coarse")

if LEGEND_STYLE == "dot":
    _draw_legend_dot(
        devstage_color_map,
        DEVSTAGE_ORDER,
        OUTPUT_DIR / "Fig1_D_1_umap_devstage_legend.svg",
        counts_map=devstage_counts_map,
    )
    _draw_legend_dot(
        tissue_color_map,
        tissue_ordered,
        OUTPUT_DIR / "Fig1_D_1_umap_tissue_legend.svg",
        counts_map=tissue_counts_map,
    )
else:
    # bar style (old)
    _draw_legend_with_bar(
        devstage_color_map,
        DEVSTAGE_ORDER,
        devstage_counts_map,
        OUTPUT_DIR / "Fig1_D_1_umap_devstage_legend.svg",
    )
    _draw_legend_with_bar(
        tissue_color_map,
        tissue_ordered,
        tissue_counts_map,
        OUTPUT_DIR / "Fig1_D_1_umap_tissue_legend.svg",
    )

print("\nDone!")
