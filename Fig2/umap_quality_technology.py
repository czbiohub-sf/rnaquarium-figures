#!/usr/bin/env python3
"""
Figure 1 Panel D.2: Transcriptome clustering UMAPs — data quality / technology views

Creates UMAP plots colored by:
  1. Seq-detective filtering outcome (BB, BT, TB, TT for PE; B, T for SE)
  2. Sequencing technology category (Smart-seq, 10x Genomics, etc.)

Reuses X_umap pre-computed in the anndata object (same embedding as D.1).
Outcome colors are consistent with Fig1_B_2 OUTCOME_COLORS / Fig1_B_3 BATLOW_COLORS.
Technology colors use the same batlow palette and category mapping as Fig1_B_2.
"""

import numpy as np
import polars as pl
import anndata as ad
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import scanpy as sc
from pathlib import Path

# =============================================================================
# Configuration
# =============================================================================

ANNDATA_H5AD = Path("data/75k_unstable/75k_unstable_anndata_zfin_aliases_metadata.log2tmmcpm.h5ad")
SEQDETECTIVE_FILE = Path("data/75k_unstable/seq-detective-judgement-summary-augmented.txt")
ANNOTATIONS_FILE = Path("data/zf-core-v2-74K_problematic_with_bulk.kmers.allcols.csv")
BATLOW_PALETTE = Path("palette/batlow/DiscretePalettes/batlow100.txt")
OUTPUT_DIR = Path("figures")
OUTPUT_DIR.mkdir(exist_ok=True)

CM = 1 / 2.54
FIG_W = 6.44 * CM   # panel D2: same slot as D1 (6.44 × 7.125 cm)
FIG_H = 7.125 * CM

# Font settings (Arial, 5–7pt per layout spec)
mpl.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['font.sans-serif'] = ['Arial', 'Helvetica Neue', 'DejaVu Sans']
FONT_TITLE = 7
FONT_LABEL = 7
FONT_LEGEND = 5
FONT_TICK = 5

# Seq-detective outcome colors — authoritative source: Fig1_B_3 BATLOW_COLORS
# PE combos match exactly; SE (B/T) are close variants of BB/TT.
OUTCOME_COLORS = {
    "BB":      "#99882c",  # both biological   — olive        (= B-B in Fig1_B_3)
    "TB":      "#426f52",  # M1 tech, M2 bio  — forest green  (= T-B in Fig1_B_3)
    "BT":      "#0b2c5c",  # M1 bio, M2 tech  — dark navy     (= B-T in Fig1_B_3)
    "TT":      "#f29d6c",  # both technical   — coral/peach   (= T-T in Fig1_B_3)
    "B":       "#c4b030",  # SE biological    — lighter olive (close to BB)
    "T":       "#d4603c",  # SE technical     — darker coral  (close to TT)
    "Unknown": "#9E9E9E",  # no seq-detective data
}

# =============================================================================
# Palette loading
# =============================================================================

def load_batlow_colors(n_colors):
    """Load n_colors from batlow100 discrete palette using golden-ratio stride."""
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
    stride = int(len(colors) * 0.618)
    return [colors[(i * stride) % len(colors)] for i in range(n_colors)]

# =============================================================================
# Load anndata
# =============================================================================

print(f"Loading anndata from {ANNDATA_H5AD} ...")
adata = ad.read_h5ad(ANNDATA_H5AD)
print(f"  {adata.n_obs:,} samples × {adata.n_vars:,} genes")
assert "X_umap" in adata.obsm, "X_umap not found — run Fig1_D_1_umap.py first to compute it"

# =============================================================================
# Load seq-detective outcomes
# =============================================================================

seqdetective = pl.read_csv(
    SEQDETECTIVE_FILE,
    separator="\t",
    has_header=False,
    new_columns=["id", "file1", "file2", "grade1", "grade2", "reason"],
).with_columns(
    accession=pl.col("id").str.strip_suffix("_subsample"),
    # PE: grade1+grade2 (e.g. "BB"); SE: grade1 only (null grade2 is ignored)
    outcome=pl.concat_str(pl.col("grade1"), pl.col("grade2"), ignore_nulls=True),
)

print(f"Loaded {len(seqdetective)} seq-detective judgements")

# =============================================================================
# Load technology annotations (same category mapping as Fig1_B_2)
# =============================================================================

annotations = pl.read_csv(
    ANNOTATIONS_FILE,
    columns=["accession", "technology"],
).with_columns(
    tech_category=pl.when(pl.col("technology").str.contains("Smart-seq|SMART-seq"))
        .then(pl.lit("Smart-seq"))
        .when(pl.col("technology").str.contains("DeTCT"))
        .then(pl.lit("DeTCT"))
        .when(pl.col("technology").str.contains("10x"))
        .then(pl.lit("10x Genomics"))
        .when(pl.col("technology").str.contains("CEL-Seq"))
        .then(pl.lit("CEL-Seq"))
        .when(pl.col("technology").str.contains("CORALL"))
        .then(pl.lit("CORALL"))
        .when(pl.col("technology").str.contains("sci-RNA-seq"))
        .then(pl.lit("sci-RNA-seq"))
        .when(pl.col("technology").str.contains("inDrop"))
        .then(pl.lit("inDrop"))
        .when(pl.col("technology").str.contains("CLIP|RIP|ChIP"))
        .then(pl.lit("Protein-RNA"))
        .otherwise(pl.lit("Other"))
)

# =============================================================================
# Join metadata onto adata.obs
# =============================================================================

obs_ids = adata.obs_names.tolist()
obs_base = pl.DataFrame({"accession": obs_ids})

outcome_joined = obs_base.join(
    seqdetective.select("accession", "outcome"),
    on="accession",
    how="left",
).with_columns(pl.col("outcome").fill_null("Unknown"))

tech_joined = obs_base.join(
    annotations.select("accession", "tech_category"),
    on="accession",
    how="left",
).with_columns(pl.col("tech_category").fill_null("Unknown"))

adata.obs["sd_outcome"] = outcome_joined["outcome"].to_list()
adata.obs["sd_outcome"] = adata.obs["sd_outcome"].astype("category")

adata.obs["tech_category"] = tech_joined["tech_category"].to_list()
adata.obs["tech_category"] = adata.obs["tech_category"].astype("category")

# Summary
print("\nSeq-detective outcome distribution (UMAP samples):")
print(adata.obs["sd_outcome"].value_counts().sort_values(ascending=False))
print("\nTechnology distribution (UMAP samples):")
print(adata.obs["tech_category"].value_counts().sort_values(ascending=False))

# =============================================================================
# Assign colors
# =============================================================================

# Outcomes: fixed palette; order categories by count desc
outcome_cats = list(adata.obs["sd_outcome"].value_counts().index)
adata.obs["sd_outcome"] = adata.obs["sd_outcome"].cat.reorder_categories(outcome_cats)
adata.uns["sd_outcome_colors"] = [OUTCOME_COLORS.get(c, "#808080") for c in outcome_cats]

# Technology: batlow by count rank; "Other" and "Unknown" get neutral grays
_tech_vc = adata.obs["tech_category"].value_counts().sort_values(ascending=False)
_tech_ordered = list(_tech_vc.index)
_tech_named = [c for c in _tech_ordered if c not in ("Other", "Unknown")]
_batlow = load_batlow_colors(len(_tech_named))
_tech_color_map = {c: _batlow[i] for i, c in enumerate(_tech_named)}
_tech_color_map["Other"]   = "#BDBDBD"
_tech_color_map["Unknown"] = "#9E9E9E"

tech_cats = _tech_ordered  # already count-desc
adata.obs["tech_category"] = adata.obs["tech_category"].cat.reorder_categories(tech_cats)
adata.uns["tech_category_colors"] = [_tech_color_map.get(c, "#808080") for c in tech_cats]

# =============================================================================
# Figure styling
# =============================================================================

mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"] = 42
mpl.rcParams["svg.fonttype"] = "none"
mpl.rcParams["figure.figsize"] = [FIG_W, FIG_H]

# =============================================================================
# UMAP: seq-detective outcome
# =============================================================================

print("\nPlotting UMAP: seq-detective outcome...")

fig = sc.pl.umap(
    adata,
    color="sd_outcome",
    title="Seq-Detective Filtering Outcome",
    legend_loc=None,
    add_outline=True,
    outline_width=(0.05, 0.05),
    outline_color=("black", "white"),
    size=3,
    alpha=0.8,
    show=False,
    return_fig=True,
)
out_outcomes = OUTPUT_DIR / "Fig1_D_2_umap_sd_outcomes.svg"
fig.savefig(out_outcomes, transparent=True)
print(f"Saved {out_outcomes}")
plt.close()

# =============================================================================
# UMAP: sequencing technology
# =============================================================================

print("Plotting UMAP: sequencing technology...")

fig = sc.pl.umap(
    adata,
    color="tech_category",
    title="Sequencing Technology",
    legend_loc=None,
    add_outline=True,
    outline_width=(0.05, 0.05),
    outline_color=("black", "white"),
    size=3,
    alpha=0.8,
    show=False,
    return_fig=True,
)
out_tech = OUTPUT_DIR / "Fig1_D_2_umap_technology.svg"
fig.savefig(out_tech, transparent=True)
print(f"Saved {out_tech}")
plt.close()

print("\nDone!")
