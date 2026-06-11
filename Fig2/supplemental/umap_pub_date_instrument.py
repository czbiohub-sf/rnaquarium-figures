#!/usr/bin/env python3
"""
UMAP panels colored by publication date and sequencing instrument,
plus DeTCT separation analysis.

Panels produced:
  1. UMAP colored by publication year (continuous, sequential palette)
  2. UMAP colored by sequencing instrument model (batlow discrete)
  3. UMAP with DeTCT studies highlighted
  4. Top genes / QC features distinguishing DeTCT clusters (bar chart)

DeTCT studies are identified by technology == 'DeTCT' in the annotations file.
The user-supplied PRJEB IDs (PRJEB6584, PRJEB7799, PRJEB8827, PRJEB7614,
PRJEB7244) are BioProject aliases for the same ERP studies.
"""

import numpy as np
import polars as pl
import anndata as ad
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.cm as cm
import scanpy as sc
import pandas as pd
from pathlib import Path

# =============================================================================
# Configuration
# =============================================================================

ANNDATA_H5AD = Path("../Fig1/data/75k_unstable/75k_unstable_anndata_zfin_aliases_metadata.log2tmmcpm.h5ad")
ANNOTATIONS_FILE = Path("../Fig1/data/zf-core-v2-74K_problematic_with_bulk.kmers.allcols.csv")
BATLOW_PALETTE = Path("../Fig1/palette/batlow/DiscretePalettes/batlow100.txt")
OUTPUT_DIR = Path("figures")
OUTPUT_DIR.mkdir(exist_ok=True)

CM = 1 / 2.54
FIG_W = 6.44 * CM
FIG_H = 7.125 * CM

mpl.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['font.sans-serif'] = ['Arial', 'Helvetica Neue', 'DejaVu Sans']
mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"] = 42
mpl.rcParams["svg.fonttype"] = "none"
FONT_TITLE = 7
FONT_LABEL = 7
FONT_LEGEND = 5
FONT_TICK = 5

# Top instruments to show individually; rest collapsed to "Other"
TOP_N_INSTRUMENTS = 8

# =============================================================================
# Helpers
# =============================================================================

def load_batlow_colors(n_colors):
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


def save_umap(adata, color_key, title, outpath, **kwargs):
    fig = sc.pl.umap(
        adata,
        color=color_key,
        title=title,
        legend_loc=None,
        add_outline=True,
        outline_width=(0.05, 0.05),
        outline_color=("black", "white"),
        size=3,
        alpha=0.8,
        show=False,
        return_fig=True,
        **kwargs,
    )
    fig.savefig(outpath, transparent=True)
    print(f"Saved {outpath}")
    plt.close()

# =============================================================================
# Load data
# =============================================================================

print(f"Loading anndata from {ANNDATA_H5AD} ...")
adata = ad.read_h5ad(ANNDATA_H5AD)
print(f"  {adata.n_obs:,} samples × {adata.n_vars:,} genes")
assert "X_umap" in adata.obsm, "X_umap not found — run umap compute script first"

annotations = pl.read_csv(
    ANNOTATIONS_FILE,
    null_values=["#DIV/0!", "", "/"],
    ignore_errors=True,
).select(["accession", "study", "technology", "instrument_model", "published_date"])

obs_ids = adata.obs_names.tolist()
obs_base = pl.DataFrame({"accession": obs_ids})
ann_joined = obs_base.join(annotations, on="accession", how="left")

# =============================================================================
# Publication year
# =============================================================================

pub_year = (
    ann_joined
    .with_columns(
        pub_year=pl.col("published_date")
        .str.slice(0, 4)          # "YYYY"
        .cast(pl.Int32, strict=False)
    )
    ["pub_year"]
    .to_list()
)

adata.obs["pub_year"] = pd.array(pub_year, dtype="Int32").astype("float")  # NaN for unknowns

year_min = int(np.nanmin(adata.obs["pub_year"]))
year_max = int(np.nanmax(adata.obs["pub_year"]))
print(f"Publication year range: {year_min}–{year_max}")

# =============================================================================
# Instrument model
# =============================================================================

top_instruments = (
    ann_joined["instrument_model"]
    .drop_nulls()
    .value_counts()
    .sort("count", descending=True)
    .head(TOP_N_INSTRUMENTS)["instrument_model"]
    .to_list()
)

instrument_cat = (
    ann_joined
    .with_columns(
        instrument=pl.when(pl.col("instrument_model").is_in(top_instruments))
        .then(pl.col("instrument_model"))
        .otherwise(pl.lit("Other / Unknown"))
    )
    ["instrument"]
    .fill_null("Other / Unknown")
    .to_list()
)

adata.obs["instrument"] = pd.Categorical(instrument_cat)

# Assign colors: named instruments get batlow, Other/Unknown gets gray
_inst_named = [c for c in top_instruments]
_batlow = load_batlow_colors(len(_inst_named))
_inst_color_map = {c: _batlow[i] for i, c in enumerate(_inst_named)}
_inst_color_map["Other / Unknown"] = "#BDBDBD"

inst_cats = list(adata.obs["instrument"].value_counts().index)
adata.obs["instrument"] = adata.obs["instrument"].cat.reorder_categories(
    sorted(inst_cats, key=lambda c: -adata.obs["instrument"].value_counts()[c])
)
adata.uns["instrument_colors"] = [
    _inst_color_map.get(c, "#9E9E9E") for c in adata.obs["instrument"].cat.categories
]

print("Instrument distribution:")
print(adata.obs["instrument"].value_counts().sort_values(ascending=False))

# =============================================================================
# DeTCT flag
# =============================================================================

is_detct = (
    ann_joined
    .with_columns(
        detct=pl.col("technology").str.to_lowercase().str.contains("detct")
    )
    ["detct"]
    .fill_null(False)
    .to_list()
)

adata.obs["is_detct"] = pd.Categorical(
    ["DeTCT" if v else "Other" for v in is_detct]
)
adata.uns["is_detct_colors"] = ["#e63946", "#DDDDDD"]  # red for DeTCT, light gray for Other

n_detct = sum(is_detct)
print(f"\nDeTCT samples in UMAP: {n_detct:,} / {adata.n_obs:,}")

# =============================================================================
# UMAP 1: publication year (continuous)
# =============================================================================

print("\nPlotting UMAP: publication year...")
mpl.rcParams["figure.figsize"] = [FIG_W, FIG_H]

fig = sc.pl.umap(
    adata,
    color="pub_year",
    title="Publication Year",
    color_map="viridis",
    vmin=year_min,
    vmax=year_max,
    legend_loc=None,
    add_outline=True,
    outline_width=(0.05, 0.05),
    outline_color=("black", "white"),
    size=3,
    alpha=0.8,
    show=False,
    return_fig=True,
)
fig.savefig(OUTPUT_DIR / "umap_publication_year.svg", transparent=True)
print(f"Saved {OUTPUT_DIR / 'umap_publication_year.svg'}")
plt.close()

# =============================================================================
# UMAP 2: sequencing instrument
# =============================================================================

print("Plotting UMAP: sequencing instrument...")
save_umap(adata, "instrument", "Sequencing Instrument",
          OUTPUT_DIR / "umap_instrument_model.svg")

# =============================================================================
# UMAP 3: DeTCT highlight
# =============================================================================

print("Plotting UMAP: DeTCT highlight...")
# Plot Other first (gray), DeTCT on top (red)
adata.obs["is_detct"] = adata.obs["is_detct"].cat.reorder_categories(["Other", "DeTCT"])
adata.uns["is_detct_colors"] = ["#DDDDDD", "#e63946"]

save_umap(adata, "is_detct", "DeTCT Studies",
          OUTPUT_DIR / "umap_detct_highlight.svg")

# =============================================================================
# DeTCT feature analysis
# =============================================================================
# Use scanpy's rank_genes_groups (Wilcoxon) to find genes that best separate
# DeTCT from everything else.  Then also compare a handful of QC-level numeric
# metadata features.

print("\nRunning DeTCT vs Other differential expression (Wilcoxon)...")

sc.tl.rank_genes_groups(
    adata,
    groupby="is_detct",
    groups=["DeTCT"],
    reference="Other",
    method="wilcoxon",
    n_genes=adata.n_vars,  # all genes so we can find both up and down
    key_added="detct_vs_other",
)

# Extract full results sorted by absolute score
rgg = adata.uns["detct_vs_other"]
all_genes = pd.DataFrame({
    "gene":     rgg["names"]["DeTCT"],
    "score":    rgg["scores"]["DeTCT"],
    "log2fc":   rgg["logfoldchanges"]["DeTCT"],
    "pval_adj": rgg["pvals_adj"]["DeTCT"],
})
all_genes["abs_score"] = all_genes["score"].abs()

# Top 50 by absolute score for the bar chart / CSV
top_genes = all_genes.nlargest(50, "abs_score")

out_csv = OUTPUT_DIR / "detct_top_genes.csv"
top_genes.to_csv(out_csv, index=False)
print(f"Saved top genes to {out_csv}")

# Bar chart: top 20 genes by Wilcoxon score
print("Plotting DeTCT top-gene bar chart...")
top20 = top_genes.head(20).sort_values("score")

fig, ax = plt.subplots(figsize=(8 * CM, 10 * CM))
colors = ["#e63946" if s > 0 else "#457b9d" for s in top20["score"]]
ax.barh(top20["gene"], top20["score"], color=colors)
ax.axvline(0, color="black", linewidth=0.5)
ax.set_xlabel("Wilcoxon score (DeTCT vs Other)", fontsize=FONT_LABEL)
ax.set_title("Top genes distinguishing DeTCT", fontsize=FONT_TITLE)
ax.tick_params(labelsize=FONT_TICK)
plt.tight_layout()
out_bar = OUTPUT_DIR / "detct_top_genes_bar.svg"
fig.savefig(out_bar, transparent=True)
print(f"Saved {out_bar}")
plt.close()

# Pick top up- and down-regulated genes from the full ranked list
_up_rows = all_genes[all_genes["log2fc"] > 0].nlargest(1, "score")
_dn_rows = all_genes[all_genes["log2fc"] < 0].nsmallest(1, "score")

top_up = _up_rows["gene"].iloc[0] if not _up_rows.empty else None
top_dn = _dn_rows["gene"].iloc[0] if not _dn_rows.empty else None

if top_up:
    print(f"\nTop upregulated in DeTCT:   {top_up}  (log2FC={_up_rows['log2fc'].values[0]:.2f})")
if top_dn:
    print(f"Top downregulated in DeTCT: {top_dn}  (log2FC={_dn_rows['log2fc'].values[0]:.2f})")
else:
    print("No downregulated genes found in DeTCT vs Other comparison.")

for gene, label in [(top_up, "up"), (top_dn, "dn")]:
    if gene is None:
        continue
    if gene not in adata.var_names:
        print(f"  {gene} not in adata.var_names, skipping UMAP")
        continue
    fig = sc.pl.umap(
        adata, color=gene,
        title=f"{gene} (DeTCT {label}reg)",
        color_map="RdBu_r",
        add_outline=True,
        outline_width=(0.05, 0.05),
        outline_color=("black", "white"),
        size=3, alpha=0.8,
        show=False, return_fig=True,
    )
    out = OUTPUT_DIR / f"umap_detct_{label}reg_{gene}.svg"
    fig.savefig(out, transparent=True)
    print(f"Saved {out}")
    plt.close()

# Print summary table
print("\n=== Top 20 DeTCT-distinguishing genes ===")
print(top20[["gene", "log2fc", "score", "pval_adj"]].to_string(index=False))

print("\nDone!")
