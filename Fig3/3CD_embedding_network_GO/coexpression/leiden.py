"""
df_to_leiden_parquet_corr_knn.py  (robustified)

Additions:
  - corr_type: 'pearson' | 'spearman' | 'robust'
      * 'spearman': Pearson on per-gene ranks (rankdata), robust to magnitude outliers.
      * 'robust'  : median/MAD pre-transform, then Pearson on the transformed values
                    (implemented as: robust-scale -> standard z-score -> Pearson).
  - winsorize_q: optional per-gene winsorization BEFORE the chosen transform.
      * Example: --winsorize-q 0.01 clips each gene at its 1st/99th percentiles.

Notes:
  - 'knn' and 'pynndescent' compute cosine on per-gene standardized rows; for
    corr_type='spearman' this equals Spearman correlation; for 'robust' it is
    Pearson on the robust-transformed data.
  - 'dense' computes the dense correlation matrix from Z and then selects top-k.
  - 'corr' returns the dense correlation matrix (respecting corr_type and winsorization).

Dependencies:
  pip install pandas pyarrow scanpy igraph leidenalg scikit-learn scipy
"""

from __future__ import annotations

import warnings
from typing import Dict, Any, Optional, Tuple
import time, os

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import rankdata

# Core deps
try:
    import scanpy as sc
except Exception as e:
    raise ImportError("Scanpy is required. Install with: pip install scanpy") from e

# Leiden deps (required by scanpy.tl.leiden)
try:
    import igraph  # noqa: F401
    import leidenalg  # noqa: F401
except Exception as e:
    raise ImportError(
        "Leiden requires 'igraph' and 'leidenalg'. Install with: pip install igraph leidenalg"
    ) from e

from sklearn.neighbors import NearestNeighbors


# ---------- IO & basic helpers ----------

def _ensure_numeric_df(df: pd.DataFrame) -> pd.DataFrame:
    df_num = df.select_dtypes(include=[np.number]).copy()
    if df_num.shape[1] == 0:
        raise ValueError("No numeric columns found in the DataFrame.")
    return df_num


def _handle_nas(df_num: pd.DataFrame, strategy: str) -> pd.DataFrame:
    if strategy == "drop":
        return df_num.dropna(axis=0, how="any")
    if strategy == "zero":
        return df_num.fillna(0.0)
    if strategy == "mean":
        return df_num.fillna(df_num.mean(numeric_only=True))
    raise ValueError("na_strategy must be one of: 'drop', 'zero', 'mean'.")


def _read_parquet(parquet_path: str, index_col: Optional[str]) -> pd.DataFrame:
    df = pd.read_parquet(parquet_path)
    if index_col is not None:
        if index_col in df.columns:
            df = df.set_index(index_col)
        else:
            warnings.warn(
                f"index_col='{index_col}' not found in columns; using existing DataFrame index."
            )
    return df


# ---------- Row-wise transforms ----------

def _winsorize_rows(G: np.ndarray, q: float) -> np.ndarray:
    """Per-gene winsorization at quantiles [q, 1-q]."""
    if q <= 0.0:
        return G
    if not (0.0 < q < 0.5):
        raise ValueError("winsorize_q must be in (0, 0.5).")
    lo = np.quantile(G, q, axis=1, keepdims=True)
    hi = np.quantile(G, 1.0 - q, axis=1, keepdims=True)
    return np.clip(G, lo, hi)


def _zscore_rows_ddof1(X: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """
    Standard z-score per row using sample std (ddof=1).
    Returns:
      Z: standardized rows
      keep: boolean mask of rows with non-zero std
    """
    n = X.shape[1]
    if n < 2:
        raise ValueError("Need at least 2 samples to compute Pearson correlation.")
    means = X.mean(axis=1, keepdims=True)
    centered = X - means
    ss = np.sum(centered * centered, axis=1, keepdims=False)
    denom = np.sqrt(ss / (n - 1))
    keep = denom > 0
    if not np.all(keep):
        dropped = np.count_nonzero(~keep)
        warnings.warn(f"Dropping {dropped} genes with zero variance after transform.")
    Z = centered[keep] / denom[keep, None]
    return Z, keep


def _prepare_Z(
    G: np.ndarray,
    corr_type: str = "pearson",
    winsorize_q: float = 0.0,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Produce Z (genes x samples) suitable for correlation:
      - If corr_type='pearson': optional winsorization -> z-score (mean/std).
      - If corr_type='spearman': optional winsorization -> per-row ranks -> z-score.
      - If corr_type='robust':  optional winsorization -> (x - median)/MAD -> z-score.
    Returns (Z, keep_mask) where keep_mask refers to original rows of G.
    """
    corr_type = corr_type.lower()
    if corr_type not in ("pearson", "spearman", "robust"):
        raise ValueError("corr_type must be 'pearson', 'spearman', or 'robust'.")

    # Optional winsorization first (applies to all corr types)
    X0 = _winsorize_rows(G, winsorize_q) if winsorize_q > 0 else G

    if corr_type == "spearman":
        # Rank per gene (row). Ties -> average ranks.
        # np.apply_along_axis is fine here; for very large N you can consider numba.
        ranks = np.apply_along_axis(rankdata, 1, X0, method="average").astype(np.float32, copy=False)
        Z, keep = _zscore_rows_ddof1(ranks)
        return Z, keep

    if corr_type == "robust":
        # Robust scale per row via MAD
        med = np.median(X0, axis=1, keepdims=True).astype(np.float64)
        mad = (1.4826 * np.median(np.abs(X0 - med), axis=1, keepdims=True)).astype(np.float64)

        eps = 1e-12
        use_std = (mad[:, 0] < eps)

        # winsorized std as fallback scale (center still by median)
        std_win = np.std(X0, axis=1, ddof=1, keepdims=True).astype(np.float64)

        scale = np.where(use_std[:, None], std_win, mad)
        keep_scale = (scale[:, 0] > eps)
        if not np.all(keep_scale):
            warnings.warn(f"Using std fallback for {use_std.sum()} genes; dropping {np.count_nonzero(~keep_scale)} with zero scale.")

        Xr = (X0[keep_scale] - med[keep_scale]) / scale[keep_scale]
        # now safety in z-score step: treat 'near-zero' as zero
        def _z_ddof1_tol(X):
            n = X.shape[1]
            m = X.mean(axis=1, keepdims=True)
            C = X - m
            sd = np.sqrt((C * C).sum(axis=1, keepdims=False) / max(n - 1, 1))
            keep = sd > 1e-12
            Z = C[keep] / sd[keep, None]
            return Z, keep
        Z_sub, keep2 = _z_ddof1_tol(Xr)

        keep = np.zeros(G.shape[0], dtype=bool)
        idx = np.flatnonzero(keep_scale)
        keep[idx[keep2]] = True
        return Z_sub.astype(np.float32, copy=False), keep

    # Pearson
    Z, keep = _zscore_rows_ddof1(X0.astype(np.float32, copy=False))
    return Z, keep


# ---------- Build graphs ----------

def _build_knn_from_Z_cosine_as_corr(
    Z: np.ndarray,
    n_neighbors: int,
    min_corr: float = 0.0,
    n_jobs: int = -1,
) -> Tuple[sparse.csr_matrix, sparse.csr_matrix]:
    """
    Build a KNN graph from standardized Z (genes x samples) using cosine distance.
    For corr_type='pearson' or 'spearman', this equals Pearson/Spearman correlation.
    For corr_type='robust', this is Pearson on robust-transformed data.
    Returns (connectivities, distances) as CSR.
    """
    n_genes = Z.shape[0]
    if n_neighbors >= n_genes:
        raise ValueError("n_neighbors must be smaller than the number of genes.")

    nn = NearestNeighbors(
        n_neighbors=n_neighbors + 1,  # include self, will drop
        metric="cosine",
        algorithm="auto",
        n_jobs=n_jobs,
    ).fit(Z)

    dists, idxs = nn.kneighbors(Z, return_distance=True)
    # drop self neighbor at column 0 (distance=0)
    dists = dists[:, 1:]
    idxs = idxs[:, 1:]

    # Convert cosine distance to correlation weight: r = 1 - d
    r_vals = 1.0 - dists
    # Apply threshold
    mask = r_vals >= min_corr
    rows, cols, data = [], [], []
    for i in range(n_genes):
        mi = mask[i]
        if np.any(mi):
            rows.append(np.full(mi.sum(), i, dtype=np.int32))
            cols.append(idxs[i, mi].astype(np.int32))
            data.append(r_vals[i, mi].astype(np.float64))

    if len(data) == 0:
        raise ValueError("No edges remain after applying min_corr; try lowering --min-corr.")

    rows = np.concatenate(rows)
    cols = np.concatenate(cols)
    data = np.concatenate(data)

    A = sparse.csr_matrix((data, (rows, cols)), shape=(n_genes, n_genes), dtype=np.float32)
    # Symmetrize by max to get an undirected graph
    A = A.maximum(A.T).tocsr()
    # Distances derived from correlation weights
    D = A.copy().tocsr()
    D.data = 1.0 - D.data
    D.data = np.clip(D.data, 0.0, None)
    D.eliminate_zeros()
    return A, D


def _build_knn_from_dense_corr(
    Z: np.ndarray,
    n_neighbors: int,
    use_abs: bool = False,
    min_corr: float = 0.0,
) -> Tuple[sparse.csr_matrix, sparse.csr_matrix]:
    """
    Explicitly compute dense correlation matrix from Z and keep top-k per row.
    Supports absolute correlations (use_abs=True).
    """
    n_genes, n_samples = Z.shape
    # Because rows of Z are standardized, dot(Z, Z.T) / (n-1) is the correlation
    R = (Z @ Z.T) / (n_samples - 1)
    np.fill_diagonal(R, -np.inf)  # exclude self from top-k

    S = np.abs(R) if use_abs else R
    if min_corr > 0:
        S_mask = S >= min_corr
        if not S_mask.any():
            raise ValueError("No correlations meet min_corr; try lowering it.")

    # Top-k indices per row
    k = min(n_neighbors, n_genes - 1)
    part = np.argpartition(S, -k, axis=1)[:, -k:]  # (n_genes, k)

    # For each row, sort selected top-k indices by S descending
    row_idx = np.arange(n_genes)[:, None]
    top_vals = S[row_idx, part]
    order = np.argsort(-top_vals, axis=1)
    top_cols = part[row_idx, order]
    top_scores = S[row_idx, top_cols]  # noqa: F841 (kept for clarity)

    if use_abs:
        weights = R[row_idx, top_cols]
        weights = np.where(np.abs(weights) >= min_corr, weights, 0.0)
    else:
        weights = R[row_idx, top_cols]
        if min_corr > 0.0:
            weights = np.where(weights >= min_corr, weights, 0.0)

    keep_mask = (weights != 0.0)
    rows = np.repeat(np.arange(n_genes), k)
    cols = top_cols.ravel()
    data = weights.ravel()
    rows, cols, data = rows[keep_mask.ravel()], cols[keep_mask.ravel()], data[keep_mask.ravel()]

    if data.size == 0:
        raise ValueError("No edges left after thresholding/top-k selection.")

    A = sparse.csr_matrix((data, (rows, cols)), shape=(n_genes, n_genes), dtype=np.float32)
    # Symmetrize
    A = A.maximum(A.T).tocsr()
    D = (1.0 - A).tocsr()
    D.data = np.clip(D.data, 0.0, None)
    return A, D


def _knn_with_pynndescent(Z, n_neighbors, min_corr=0.0, random_state=0, n_jobs=-1):
    from pynndescent import NNDescent
    index = NNDescent(
        Z, metric="cosine",
        n_neighbors=n_neighbors + 1,
        random_state=random_state,
        n_jobs=n_jobs
    )
    idxs, dists = index.query(Z, k=n_neighbors + 1)
    dists, idxs = dists[:, 1:], idxs[:, 1:]        # drop self
    r_vals = 1.0 - dists
    mask = r_vals >= min_corr
    rows, cols, data = [], [], []
    for i in range(Z.shape[0]):
        mi = mask[i]
        if mi.any():
            rows.append(np.full(mi.sum(), i, dtype=np.int32))
            cols.append(idxs[i, mi].astype(np.int32))
            data.append(r_vals[i, mi].astype(np.float32))
    rows = np.concatenate(rows); cols = np.concatenate(cols); data = np.concatenate(data)
    A = sparse.csr_matrix((data, (rows, cols)), shape=(Z.shape[0], Z.shape[0]))
    A = A.maximum(A.T).tocsr()
    D = A.copy(); D.data = 1.0 - D.data; D.eliminate_zeros()
    return A, D


# ---------- Main API ----------

def leiden_from_parquet_corr_knn(
    parquet_path: str,
    *,
    index_col: Optional[str] = None,
    axis: str = "samples_by_genes",   # rows=samples, cols=genes (default)
    n_neighbors: int = 15,
    method: str = "knn",              # 'knn' | 'dense' | 'pynndescent' | 'corr'
    use_abs: bool = False,            # only used in method='dense'
    min_corr: float = 0.1,
    resolution: float = 1.0,
    random_state: int = 0,
    na_strategy: str = "drop",
    key_added: str = "leiden",
    n_jobs: int = -1,
    corr_type: str = "pearson",       # 'pearson' | 'spearman' | 'robust'
    winsorize_q: float = 0.0,         # e.g., 0.01 clips to 1st/99th per gene
) -> Dict[str, Any]:
    """
    Build a gene–gene KNN graph with correlation weights and run Leiden.

    Returns dict with:
      'labels'         : pd.Series (index = kept gene names)
      'connectivities' : csr_matrix
      'distances'      : csr_matrix
      'adata'          : AnnData (obs = genes)
      'dropped_genes'  : list of genes removed due to zero variance/MAD
      (if method='corr')
      'corr'           : np.ndarray (dense correlation on chosen corr_type)
      'kept_genes'     : np.ndarray of gene names
    """
    df = _read_parquet(parquet_path, index_col=index_col)
    if df.empty:
        raise ValueError("Input DataFrame is empty after reading the Parquet file.")

    df_num = _ensure_numeric_df(df)
    df_num = _handle_nas(df_num, na_strategy)

    # Build genes x samples matrix
    axis = axis.lower()
    if axis not in ("samples_by_genes", "genes_by_samples"):
        raise ValueError("axis must be 'samples_by_genes' or 'genes_by_samples'.")
    if axis == "samples_by_genes":
        gene_names_all = df_num.columns.astype(str).to_numpy()
        G = df_num.to_numpy(dtype=np.float64, copy=False).T  # genes x samples
    else:
        gene_names_all = df_num.index.astype(str).to_numpy()
        G = df_num.to_numpy(dtype=np.float64, copy=False)    # genes x samples

    print(f"Genes x Samples matrix: {G.shape}")
    print(f"corr_type={corr_type}, winsorize_q={winsorize_q}")

    # Prepare Z per selected correlation type; drop zero-variance/MAD rows
    Z, keep_mask = _prepare_Z(G, corr_type=corr_type, winsorize_q=winsorize_q)
    kept_genes = gene_names_all[keep_mask]
    dropped_genes = gene_names_all[~keep_mask]
    print(f"Kept genes: {Z.shape[0]} | Dropped (zero scale): {dropped_genes.size}")

    # Correlation-only mode: compute dense correlation and return, skip Leiden
    t0 = time.time()
    method = method.lower()
    if method == "corr":
        n_genes, n_samples = Z.shape
        R = (Z @ Z.T) / (n_samples - 1)
        np.fill_diagonal(R, 1.0)
        t1 = time.time()
        print(f"Dense correlation computed in {t1 - t0:.2f} seconds; shape={R.shape}")
        return {
            "corr": R,
            "kept_genes": kept_genes,
            "dropped_genes": dropped_genes.tolist(),
            "labels": None,
            "connectivities": None,
            "distances": None,
            "adata": None,
        }

    # Build KNN correlation graph
    t0 = time.time()
    if method == "pynndescent":
        print("Using pynndescent for KNN graph construction")
        A, D = _knn_with_pynndescent(
            Z,
            n_neighbors=n_neighbors,
            min_corr=min_corr,
            random_state=random_state,
            n_jobs=n_jobs,
        )
    elif method == "knn":
        if use_abs:
            warnings.warn("use_abs=True is ignored for method='knn'. Use method='dense' for |r| support.")
        print("Using scikit-learn for KNN graph construction")
        A, D = _build_knn_from_Z_cosine_as_corr(
            Z,
            n_neighbors=n_neighbors,
            min_corr=min_corr,
            n_jobs=n_jobs,
        )
    elif method == "dense":
        A, D = _build_knn_from_dense_corr(
            Z,
            n_neighbors=n_neighbors,
            use_abs=use_abs,
            min_corr=min_corr,
        )
    else:
        raise ValueError("method must be 'knn', 'dense', 'pynndescent', or 'corr'.")
    t1 = time.time()
    print(f"KNN graph construction time: {t1-t0:.2f} seconds")
    print(f"connectivities shape: {A.shape}, nnz={A.nnz}")
    print(f"distances shape: {D.shape}, nnz={D.nnz}")

    # Pack into AnnData and run Leiden
    adata = sc.AnnData(X=np.zeros((A.shape[0], 1), dtype=np.float32))
    adata.obs_names = kept_genes
    adata.obsp["connectivities"] = A.tocsr()
    adata.obsp["distances"] = D.tocsr()
    adata.uns["neighbors"] = {
        "connectivities_key": "connectivities",
        "distances_key": "distances",
        "params": {
            "n_neighbors": n_neighbors,
            "metric": "cosine_on_standardized_rows",
            "method": method,
            "min_corr": float(min_corr),
            "use_abs": bool(use_abs),
            "corr_type": corr_type,
            "winsorize_q": float(winsorize_q),
        },
    }

    sc.tl.leiden(
        adata,
        resolution=resolution,
        random_state=random_state,
        key_added=key_added,
        adjacency=adata.obsp["connectivities"],
    )
    print("leiden done")

    labels = adata.obs[key_added].copy()
    labels.name = key_added

    return {
        "labels": labels,
        "connectivities": adata.obsp["connectivities"],
        "distances": adata.obsp["distances"],
        "adata": adata,
        "dropped_genes": dropped_genes.tolist(),
    }


# ---------------- CLI ----------------
if __name__ == "__main__":
    import argparse

    p = argparse.ArgumentParser(
        description="Read Parquet, build a gene–gene KNN graph from correlation (Pearson/Spearman/Robust), and run Leiden."
    )
    p.add_argument("--parquet", required=True, help="Path to Parquet table.")
    p.add_argument("--index-col", type=str, default=None, help="Column name to use as index (optional).")
    p.add_argument("--axis", type=str, default="samples_by_genes",
                   choices=["samples_by_genes", "genes_by_samples"],
                   help="Orientation of the table. Default: rows=samples, cols=genes.")
    p.add_argument("--n-neighbors", type=int, default=15, help="k for KNN (default: 15).")
    p.add_argument("--method", type=str, default="pynndescent",
                   choices=["knn", "dense", "pynndescent", "corr"],
                   help="KNN construction method. 'knn' avoids dense corr matrix. 'pynndescent' uses approximation and can be faster. 'dense' builds from dense correlation. 'corr' computes the dense correlation only and skips KNN & Leiden.")
    p.add_argument("--abs", dest="use_abs", action="store_true",
                   help="Use absolute correlation |r| (only with --method dense).")
    p.add_argument("--min-corr", type=float, default=0.1,
                   help="Drop edges with correlation < min_corr (default: 0.1).")
    p.add_argument("--resolution", type=float, default=1.0, help="Leiden resolution (default: 1.0).")
    p.add_argument("--na-strategy", choices=["drop", "zero", "mean"], default="drop",
                   help="NaN handling for the input table (default: drop rows with any NaN).")
    p.add_argument("--key", type=str, default="leiden", help="Column name for labels (default: 'leiden').")
    p.add_argument("--out", type=str, default=None, help="Write labels to this CSV path (optional).")
    p.add_argument("--npz-prefix", type=str, default=None,
                   help="If set, save *_connectivities.npz and *_distances.npz.")
    p.add_argument("--corr-out", type=str, default=None,
                   help="If method='corr', write dense correlation to this .npz file (arrays: 'corr', 'genes').")
    p.add_argument("--n-jobs", type=int, default=-1,
                   help="Number of jobs for knn/pynndescent (default: -1, all cores).")
    p.add_argument("--corr-type", type=str, default="pearson",
                   choices=["pearson", "spearman", "robust"],
                   help="Correlation type: 'pearson' (default), 'spearman' (rank-based), or 'robust' (median/MAD pre-transform).")
    p.add_argument("--winsorize-q", type=float, default=0.0,
                   help="Optional per-gene winsorization quantile (0.0 disables). Example: 0.01 clips to 1st/99th percentiles.")

    args = p.parse_args()

    res = leiden_from_parquet_corr_knn(
        parquet_path=args.parquet,
        index_col=args.index_col,
        axis=args.axis,
        n_neighbors=args.n_neighbors,
        method=args.method,
        use_abs=args.use_abs,
        min_corr=args.min_corr,
        resolution=args.resolution,
        random_state=0,
        na_strategy=args.na_strategy,
        key_added=args.key,
        n_jobs=args.n_jobs,
        corr_type=args.corr_type,
        winsorize_q=args.winsorize_q,
    )

    if args.method == "corr":
        if not args.corr_out:
            raise SystemExit("When --method corr, you must provide --corr-out <path>.npz")
        os.makedirs(os.path.dirname(args.corr_out), exist_ok=True)
        np.savez_compressed(args.corr_out, corr=res["corr"], genes=res["kept_genes"])
        print(f"Wrote correlation matrix to: {args.corr_out} (arrays: 'corr', 'genes')")
        print("Done.")
    else:
        labels = res["labels"]
        if args.out:
            labels.to_csv(args.out, header=True)
            print(f"Wrote labels to: {args.out}")

        if args.npz_prefix:
            sparse.save_npz(f"{args.npz_prefix}_connectivities.npz", res["connectivities"])
            sparse.save_npz(f"{args.npz_prefix}_distances.npz", res["distances"])
            print(f"Wrote graphs to: {args.npz_prefix}_connectivities.npz/_distances.npz")

        print("Done.")
        print(f"  Genes kept: {labels.shape[0]}")
        print(f"  Unique Leiden clusters: {labels.nunique()}")
