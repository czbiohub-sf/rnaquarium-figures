#!/usr/bin/env python3
"""
UMAP colored by max(r1_length, r2_length) per sample — proxy for the biological
read length.  For DeTCT the short oligo-dT tag is read 1; max picks the longer
mate 2 automatically.  For standard SE runs r2 is null so max = r1.
"""

import numpy as np
import pandas as pd
import anndata as ad
import matplotlib as mpl
import matplotlib.pyplot as plt
import scanpy as sc
from pathlib import Path

ANNDATA_H5AD = Path("../Fig1/data/75k_unstable/75k_unstable_anndata_zfin_aliases_metadata.log2tmmcpm.h5ad")
OUTPUT_DIR      = Path("figures")
OUTPUT_DIR.mkdir(exist_ok=True)

CM = 1 / 2.54
mpl.rcParams['font.family']     = 'sans-serif'
mpl.rcParams['font.sans-serif'] = ['Arial', 'Helvetica Neue', 'DejaVu Sans']
mpl.rcParams["pdf.fonttype"]    = 42
mpl.rcParams["svg.fonttype"]    = "none"
FONT_TITLE = 7
FONT_LABEL = 6
FONT_TICK  = 5

print(f"Loading {ANNDATA_H5AD} ...")
adata = ad.read_h5ad(ANNDATA_H5AD)
print(f"  {adata.n_obs:,} samples × {adata.n_vars:,} genes")
assert "X_umap" in adata.obsm

# r1/r2 lengths are already in obs — coerce to numeric
r1 = pd.to_numeric(adata.obs["run.r1_length"], errors="coerce")
r2 = pd.to_numeric(adata.obs["run.r2_length"], errors="coerce")

print(f"r1 length — median: {r1.median():.0f}  range: {r1.min():.0f}–{r1.max():.0f}  n_null: {r1.isna().sum()}")
print(f"r2 length — median: {r2.median():.0f}  range: {r2.min():.0f}–{r2.max():.0f}  n_null: {r2.isna().sum()}")

max_len = pd.concat([r1, r2], axis=1).max(axis=1)  # NaN-safe: ignores nulls
adata.obs["max_read_length"] = max_len.values

print(f"max length — median: {max_len.median():.0f}  range: {max_len.min():.0f}–{max_len.max():.0f}")

vmin = float(np.nanpercentile(max_len.dropna(), 2))
vmax = float(np.nanpercentile(max_len.dropna(), 98))

mpl.rcParams["figure.figsize"] = [6.44 * CM, 7.125 * CM]
fig = sc.pl.umap(
    adata,
    color="max_read_length",
    title="Max read length (bp)",
    color_map="plasma",
    vmin=vmin,
    vmax=vmax,
    add_outline=True,
    outline_width=(0.05, 0.05),
    outline_color=("black", "white"),
    size=3,
    alpha=0.8,
    show=False,
    return_fig=True,
)
out = OUTPUT_DIR / "umap_max_readlength.svg"
fig.savefig(out, transparent=True)
print(f"Saved {out}")
plt.close()

print("\nDone!")
