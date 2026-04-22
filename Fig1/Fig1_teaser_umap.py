#!/usr/bin/env python3
"""
Figure 1 teaser strip: tiny UMAP thumbnail

Small UMAP of the 75k zebrafish transcriptome (~61k samples with expression)
drawn at teaser-strip dimensions. Serves as a forward reference to Figure 2,
where the full-size UMAPs live. Points are colored by devstage_curation or
tissue_curation_coarse; no legend, no axes.

Reuses pre-computed X_umap from the anndata object (same embedding as the
full Figure 2 UMAPs under ../Fig2/).
"""

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

ANNDATA_H5AD = Path("data/75k_unstable/75k_unstable_anndata_zfin_aliases_metadata.log2tmmcpm.h5ad")
METADATA_FILE = Path("/hpc/projects/balla_group/sra_experiments/SRA_metadata/dec2025_75k_submitteradded/all_zf_dates_devstage_tissue_tech_curated.tsv")
BATLOW_PALETTE = Path("palette/batlow/DiscretePalettes/batlow100.txt")
OUTPUT_DIR = Path("figures")
OUTPUT_DIR.mkdir(exist_ok=True)

CM = 1 / 2.54
FIG_W = 6.00 * CM   # teaser umap: 6.00 × 6.00 cm
FIG_H = 6.00 * CM

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
    "All anatomical structures": "#BDBDBD",
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
    """Sample n_colors from batlow100 using golden-ratio stride (maximally distinct)."""
    colors = load_batlow_palette()
    stride = int(len(colors) * 0.618)
    return [colors[(i * stride) % len(colors)] for i in range(n_colors)]


def build_devstage_palette(cats_by_count):
    """
    Devstage palette: fixed gradient colors in chronological order.
    Returns dict: category → hex color.
    """
    return dict(zip(DEVSTAGE_ORDER, DEVSTAGE_COLORS))


def build_tissue_palette(cats_by_count):
    """
    Tissue palette: golden-ratio stride through batlow, gray overrides for
    uninformative categories.
    """
    non_gray = [c for c in cats_by_count if c not in GRAY_CATS]
    batlow = load_batlow_colors_stride(len(non_gray))
    color_map = {c: batlow[i] for i, c in enumerate(non_gray)}
    color_map.update(GRAY_CATS)
    return color_map


# =============================================================================
# Load anndata
# =============================================================================

print(f"Loading anndata from {ANNDATA_H5AD} ...")
adata = ad.read_h5ad(ANNDATA_H5AD)
print(f"  {adata.n_obs:,} samples × {adata.n_vars:,} genes")
print(f"  X_umap pre-computed: {'X_umap' in adata.obsm}")

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
    """Remove borders/ticks and fill figure bounds. No axis labels for teaser."""
    for ax in fig.axes:
        ax.set_xlabel("")
        ax.set_ylabel("")
        ax.tick_params(left=False, bottom=False, labelleft=False, labelbottom=False)
        for spine in ax.spines.values():
            spine.set_visible(False)
        ax.title.set_fontsize(FONT_TITLE)
    fig.tight_layout(pad=0)


# =============================================================================
# UMAP plots
# =============================================================================

print("\nPlotting UMAP: developmental stage...")

fig = sc.pl.umap(
    adata,
    color=key_devstage,
    title="",
    legend_loc=None,
    add_outline=True,
    outline_width=(0.05, 0.05),
    outline_color=("black", "white"),
    size=3,
    alpha=0.8,
    show=False,
    return_fig=True,
)
_apply_umap_style(fig)
out_devstage = OUTPUT_DIR / "Fig1_teaser_umap_devstage.svg"
fig.savefig(out_devstage, transparent=True)
print(f"Saved {out_devstage}")
plt.close()

print("Plotting UMAP: tissue type...")

fig = sc.pl.umap(
    adata,
    color=key_tissue,
    title="",
    legend_loc=None,
    add_outline=True,
    outline_width=(0.05, 0.05),
    outline_color=("black", "white"),
    size=3,
    alpha=0.8,
    show=False,
    return_fig=True,
)
_apply_umap_style(fig)
out_tissue = OUTPUT_DIR / "Fig1_teaser_umap_tissue.svg"
fig.savefig(out_tissue, transparent=True)
print(f"Saved {out_tissue}")
plt.close()


print("\nDone!")
