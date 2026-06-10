#!/usr/bin/env python3
"""
UMAP colored by number of detected genes per sample.

Gene detection threshold: per-sample Otsu's method applied to the log2 TMM-CPM
histogram.  Otsu maximises between-class variance across two populations (noise
lobe near 0, signal lobe above) and is analytically equivalent to a 2-component
equal-variance GMM decision boundary — but runs in microseconds per sample via
fully vectorised numpy operations.

Algorithm
---------
1. Build per-sample histograms in chunks (dense slice → searchsorted → bincount),
   fully vectorised — no per-sample Python loop.
2. Apply Otsu's formula across all 61k samples in one numpy broadcast.
3. Count genes above each sample's own threshold.

Outputs
-------
  figures/umap_ngenes_detected.svg          UMAP panel
  figures/gmm_threshold_diagnostics.svg     threshold distribution + 9 example fits
  figures/ngenes_threshold_comparison.svg   naive (>0) vs Otsu count scatter
"""

import numpy as np
import scipy.sparse as sp
import anndata as ad
import matplotlib as mpl
import matplotlib.pyplot as plt
import scanpy as sc
from pathlib import Path

ANNDATA_H5AD = Path("../Fig1/data/75k_unstable/75k_unstable_anndata_zfin_aliases_metadata.log2tmmcpm.h5ad")
OUTPUT_DIR = Path("figures")
OUTPUT_DIR.mkdir(exist_ok=True)

CM = 1 / 2.54
FIG_W  = 6.44 * CM
FIG_H  = 7.125 * CM
FONT_TITLE = 7
FONT_LABEL = 7
FONT_TICK  = 5

mpl.rcParams['font.family']     = 'sans-serif'
mpl.rcParams['font.sans-serif'] = ['Arial', 'Helvetica Neue', 'DejaVu Sans']
mpl.rcParams["pdf.fonttype"]    = 42
mpl.rcParams["ps.fonttype"]     = 42
mpl.rcParams["svg.fonttype"]    = "none"

N_BINS  = 300   # histogram resolution — fine enough for log2CPM range ~0–15
CHUNK   = 2000  # samples per dense batch (~160 MB peak per chunk)

# =============================================================================
# Load
# =============================================================================

print(f"Loading {ANNDATA_H5AD} ...")
adata = ad.read_h5ad(ANNDATA_H5AD)
print(f"  {adata.n_obs:,} samples × {adata.n_vars:,} genes")
assert "X_umap" in adata.obsm, "X_umap not found"

X = adata.X
if sp.issparse(X):
    if not sp.isspmatrix_csr(X):
        X = X.tocsr()
    is_sparse = True
else:
    X = np.asarray(X)
    is_sparse = False

n_samples, n_genes = X.shape

# =============================================================================
# Step 1 — vectorised per-sample histograms
# =============================================================================

# Bin edges span 0 → max expressed value; zeros (unexpressed) fall in bin 0.
if is_sparse:
    expr_max = float(X.data.max()) if X.nnz else 15.0
    nnz_per_row     = np.diff(X.indptr)
else:
    expr_max = float(X.max())
    nnz_per_row     = (X > 0).sum(axis=1).astype(np.int32)

bin_edges   = np.linspace(0.0, expr_max, N_BINS + 1)
bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2.0

n_zeros_per_row = n_genes - nnz_per_row    # structural zeros → bin 0

histograms = np.zeros((n_samples, N_BINS), dtype=np.int32)

print(f"Building histograms in chunks of {CHUNK:,} ...")
for start in range(0, n_samples, CHUNK):
    end = min(start + CHUNK, n_samples)
    sz  = end - start

    # Densify chunk: (sz, n_genes) float32
    chunk = X[start:end].toarray().astype(np.float32) if is_sparse else X[start:end].astype(np.float32)

    # Bin every value: searchsorted on bin_edges[1:] gives 0-based bin index
    bin_idx = np.searchsorted(bin_edges[1:], chunk, side='left')   # (sz, n_genes)
    bin_idx = np.clip(bin_idx, 0, N_BINS - 1).astype(np.int32)

    # Vectorised 2-D bincount via flat index trick (no per-row Python loop)
    sample_offsets = np.arange(sz, dtype=np.int32)[:, None] * N_BINS
    flat_idx = (bin_idx + sample_offsets).ravel()
    counts   = np.bincount(flat_idx, minlength=sz * N_BINS)
    histograms[start:end] = counts.reshape(sz, N_BINS)

    print(f"  {end:,} / {n_samples:,}", flush=True)

# Structural zeros all have value 0.0 → land in bin 0
histograms[:, 0] += n_zeros_per_row

# =============================================================================
# Step 2 — vectorised Otsu thresholds
# =============================================================================

print("Computing Otsu thresholds ...")

P = histograms.astype(np.float64)
P /= P.sum(axis=1, keepdims=True)              # normalise to probabilities

omega   = np.cumsum(P, axis=1)                 # cumulative weight  (n, B)
mu_cum  = np.cumsum(P * bin_centers, axis=1)   # cumulative mean    (n, B)
mu_T    = mu_cum[:, -1:]                       # total mean         (n, 1)

with np.errstate(divide='ignore', invalid='ignore'):
    sigma_b2 = (mu_T * omega - mu_cum) ** 2 / (omega * (1.0 - omega))
sigma_b2 = np.nan_to_num(sigma_b2, nan=0.0, posinf=0.0)

best_bin   = np.argmax(sigma_b2, axis=1)       # (n_samples,)
thresholds = bin_centers[best_bin]             # (n_samples,)

# =============================================================================
# Step 3 — count genes above per-sample threshold
# =============================================================================

print("Counting detected genes ...")
n_genes_otsu  = np.empty(n_samples, dtype=np.int32)
n_genes_naive = nnz_per_row.copy()   # genes with any stored (non-zero) value

# threshold_bin[i] = first histogram bin ABOVE threshold for sample i
threshold_bin = best_bin  # genes in bins >= best_bin are "detected"
# sum histograms from threshold_bin onward (excluding the threshold bin itself)
cumH_rev = np.cumsum(histograms[:, ::-1], axis=1)[:, ::-1]   # suffix sums
for i in range(n_samples):
    n_genes_otsu[i] = int(cumH_rev[i, best_bin[i] + 1]) if best_bin[i] + 1 < N_BINS else 0

# Vectorised version of the loop above (avoids 61k Python iters)
# fancy-index each row's suffix sum at (best_bin + 1), clamped to N_BINS - 1
idx_clamp = np.minimum(best_bin + 1, N_BINS - 1)
n_genes_otsu = cumH_rev[np.arange(n_samples), idx_clamp]
# Edge case: best_bin == N_BINS - 1  →  nothing detected (set to 0)
n_genes_otsu[best_bin == N_BINS - 1] = 0

print(f"Otsu threshold  — median: {np.median(thresholds):.3f}  "
      f"p5: {np.percentile(thresholds, 5):.3f}  "
      f"p95: {np.percentile(thresholds, 95):.3f}")
print(f"n_genes (Otsu)  — min: {n_genes_otsu.min():,}  "
      f"median: {int(np.median(n_genes_otsu)):,}  "
      f"max: {n_genes_otsu.max():,}")

adata.obs["n_genes_otsu"]  = n_genes_otsu.astype(float)
adata.obs["otsu_threshold"] = thresholds.astype(float)
adata.obs["n_genes_naive"]  = n_genes_naive.astype(float)

# =============================================================================
# UMAP
# =============================================================================

print("\nPlotting UMAP: Otsu-detected genes ...")
mpl.rcParams["figure.figsize"] = [FIG_W, FIG_H]

fig = sc.pl.umap(
    adata,
    color="n_genes_otsu",
    title="Genes detected (Otsu threshold)",
    color_map="viridis",
    add_outline=True,
    outline_width=(0.05, 0.05),
    outline_color=("black", "white"),
    size=3,
    alpha=0.8,
    show=False,
    return_fig=True,
)
out_umap = OUTPUT_DIR / "umap_ngenes_detected.svg"
fig.savefig(out_umap, transparent=True)
print(f"Saved {out_umap}")
plt.close()

# =============================================================================
# Diagnostic 1: threshold distribution + 9 example histograms with threshold
# =============================================================================

print("Plotting diagnostics ...")
rng = np.random.default_rng(42)
diag_idx = rng.choice(n_samples, size=9, replace=False)

fig = plt.figure(figsize=(18 * CM, 14 * CM))
gs  = fig.add_gridspec(4, 3, hspace=0.55, wspace=0.4)

ax_hist = fig.add_subplot(gs[0, :])
ax_hist.hist(thresholds, bins=80, color="#457b9d", edgecolor="none", alpha=0.85)
ax_hist.axvline(np.median(thresholds), color="#e63946", linewidth=1, linestyle="--",
                label=f"median = {np.median(thresholds):.2f}")
ax_hist.set_xlabel("Otsu threshold (log2 TMM-CPM)", fontsize=FONT_LABEL)
ax_hist.set_ylabel("Samples", fontsize=FONT_LABEL)
ax_hist.set_title("Per-sample Otsu detection threshold", fontsize=FONT_TITLE)
ax_hist.tick_params(labelsize=FONT_TICK)
ax_hist.legend(fontsize=5, frameon=False)

for plot_i, samp_i in enumerate(diag_idx):
    r, c = 1 + plot_i // 3, plot_i % 3
    ax = fig.add_subplot(gs[r, c])
    ax.bar(bin_centers, histograms[samp_i], width=bin_centers[1] - bin_centers[0],
           color="#adb5bd", edgecolor="none", alpha=0.7)
    ax.axvline(thresholds[samp_i], color="#e63946", linewidth=1, linestyle="--")
    ax.set_title(f"sample {samp_i}  thr={thresholds[samp_i]:.2f}", fontsize=5)
    ax.tick_params(labelsize=4)
    ax.set_xlabel("log2 TMM-CPM", fontsize=4)

out_diag = OUTPUT_DIR / "otsu_threshold_diagnostics.svg"
fig.savefig(out_diag, transparent=True)
print(f"Saved {out_diag}")
plt.close()

# =============================================================================
# Diagnostic 2: naive vs Otsu scatter
# =============================================================================

fig, ax = plt.subplots(figsize=(7 * CM, 7 * CM))
ax.scatter(n_genes_naive, n_genes_otsu, s=0.5, alpha=0.15,
           color="#457b9d", linewidths=0)
lim = max(n_genes_naive.max(), n_genes_otsu.max())
ax.plot([0, lim], [0, lim], color="#e63946", linewidth=0.8,
        linestyle="--", label="y = x")
ax.set_xlabel("Naive count  (log2CPM > 0)", fontsize=FONT_LABEL)
ax.set_ylabel("Otsu count  (above noise lobe)", fontsize=FONT_LABEL)
ax.set_title("Naive vs Otsu gene detection", fontsize=FONT_TITLE)
ax.tick_params(labelsize=FONT_TICK)
ax.legend(fontsize=5, frameon=False)
plt.tight_layout()
out_scatter = OUTPUT_DIR / "ngenes_threshold_comparison.svg"
fig.savefig(out_scatter, transparent=True)
print(f"Saved {out_scatter}")
plt.close()

print("\nDone!")
