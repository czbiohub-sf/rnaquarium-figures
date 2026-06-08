#!/usr/bin/env python3
from __future__ import annotations

import os,sys,re
import warnings
import argparse
import multiprocessing as mp
from typing import Dict, Mapping, Sequence, Tuple, Optional, List

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scanpy as sc
from concurrent.futures import ProcessPoolExecutor, as_completed
import time
from gprofiler import GProfiler
from matplotlib.lines import Line2D

from itertools import cycle, islice
import matplotlib.colors as mcolors

# -----------------------------
# Config
# -----------------------------
DEFAULT_GO_SOURCES = ["GO:BP", "GO:MF", "GO:CC"]  # biological process, molecular function, cellular component


# -----------------------------
# PCA explained variance helpers
# -----------------------------
def _pca_explained_variance_series(
    emb_df: pd.DataFrame,
    max_pcs: int = 100,
    scale: bool = False,
    pca_svd_solver: str = "arpack",
) -> pd.Series:
    """
    Compute per-PC explained variance ratio for a single embedding set.

    Returns:
        pd.Series indexed by PC number (1..k) with values = variance_ratio.
    """
    if not isinstance(emb_df.index, pd.Index):
        raise ValueError("emb_df must be a DataFrame with gene IDs as index.")

    adata = sc.AnnData(emb_df.values.astype(np.float32, copy=False))
    adata.obs_names = emb_df.index.astype(str)
    adata.var_names = [f"dim_{i}" for i in range(emb_df.shape[1])]

    if scale:
        sc.pp.scale(adata, zero_center=True, max_value=None)

    # n_comps cannot exceed rank; guard with n_obs-1 and n_vars
    n_comps = int(max(1, min(max_pcs, adata.n_vars, max(1, adata.n_obs - 1))))
    sc.tl.pca(adata, n_comps=n_comps, use_highly_variable=None, svd_solver=pca_svd_solver)

    vr = np.asarray(adata.uns["pca"]["variance_ratio"]).ravel()
    idx = np.arange(1, len(vr) + 1)
    return pd.Series(vr, index=idx, name="variance_ratio")


def compute_pca_explained_variance_across_sets(
    embedding_sets: Mapping[str, pd.DataFrame],
    max_pcs: int = 100,
    scale: bool = False,
    pca_svd_solver: str = "arpack",
) -> pd.DataFrame:
    """
    Returns a tidy DataFrame with columns:
      set, pc, variance_ratio, cumulative_variance_ratio
    """
    rows = []
    for set_name, df in embedding_sets.items():
        s = _pca_explained_variance_series(df, max_pcs=max_pcs, scale=scale, pca_svd_solver=pca_svd_solver)
        cum = s.cumsum()
        for pc in s.index:
            rows.append(
                {
                    "set": set_name,
                    "pc": int(pc),
                    "variance_ratio": float(s.loc[pc]),
                    "cumulative_variance_ratio": float(cum.loc[pc]),
                }
            )
    return pd.DataFrame(rows)


def plot_pca_explained_variance(
    ev_df: pd.DataFrame,
    out_png: Optional[str] = None,
    mode: str = "cumulative",
    vline_at: Optional[int] = None,
    title: str = "PCA explained variance across embedding sets",
    custom_palette: Optional[object] = None,   # <-- NEW
    include_sets: Optional[Sequence[str]] = None,  # <-- NEW: specify series to plot
):
    """
    Makes a single-panel line plot (no seaborn, no explicit colors).
    """
    if ev_df is None or ev_df.empty:
        warnings.warn("No explained-variance data to plot.")
        return

    plot_df = ev_df.copy()
    if include_sets:
        requested = [s for s in include_sets if s in plot_df["set"].unique()]
        if not requested:
            warnings.warn("None of the requested series were found in ev_df; nothing to plot.")
            return
        plot_df = plot_df[plot_df["set"].isin(requested)]
        labels = requested
    else:
        labels = list(plot_df["set"].unique())

    pal = _build_palette_for_labels(labels, custom_palette)

    plt.figure(figsize=(8.0, 5.0))
    for set_name in labels:
        sub = plot_df[plot_df["set"] == set_name].sort_values("pc")
        y = sub["cumulative_variance_ratio"] if mode == "cumulative" else sub["variance_ratio"]
        plt.plot(
            sub["pc"].values,
            y.values,
            marker="o",
            label=set_name,
            color=pal[set_name],               # <-- use manual color
        )

    if vline_at is not None and np.isfinite(vline_at):
        plt.axvline(int(vline_at), linestyle="--", alpha=0.6)

    plt.xlabel("Principal component #")
    plt.ylabel("Cumulative explained variance" if mode == "cumulative" else "Explained variance (per PC)")
    plt.title(title)
    plt.grid(True, linestyle="--", alpha=0.4)
    plt.legend(title="Embedding set", frameon=False)
    plt.tight_layout()
    if out_png:
        plt.savefig(out_png, dpi=150)
    plt.close()


# -----------------------------
# Clustering helpers
# -----------------------------
def _build_representation(
    emb_df: pd.DataFrame,
    cluster_space: str = "embedding",
    n_pcs: int = 0,
    scale_pca: bool = False,
    pca_svd_solver: str = "arpack",  # "arpack" or "randomized"
) -> sc.AnnData:
    """
    Create an AnnData with the representation to be used for the neighbor graph.

    cluster_space:
      - "embedding": use the provided embedding dims directly (X)
      - "pca": compute PCA and use X_pca as representation

    Returns:
      AnnData with:
        - X = original embeddings
        - obsm["X_pca"] present if cluster_space == "pca"
    """
    if not isinstance(emb_df.index, pd.Index):
        raise ValueError("emb_df must be a DataFrame with gene IDs (ENSDARG...) as index.")

    adata = sc.AnnData(emb_df.values.astype(np.float32, copy=False))
    adata.obs_names = emb_df.index.astype(str)
    adata.var_names = [f"dim_{i}" for i in range(emb_df.shape[1])]

    if cluster_space == "pca":
        if scale_pca:
            sc.pp.scale(adata, zero_center=True, max_value=None)
        n_comps = min(n_pcs, adata.n_vars)
        sc.tl.pca(adata, n_comps=n_comps, use_highly_variable=None, svd_solver=pca_svd_solver)
    elif cluster_space == "embedding":
        pass
    else:
        raise ValueError("cluster_space must be 'embedding' or 'pca'.")

    return adata


def leiden_clusters_from_embeddings(
    emb_df: pd.DataFrame,
    resolutions: Sequence[float],
    cluster_space: str = "embedding",   # {"embedding", "pca"}
    n_pcs: int = 50,
    scale_pca: bool = False,
    n_neighbors: int = 20,
    metric: Optional[str] = None,       # default chosen by space if None
    random_state: int = 0,
    pca_svd_solver: str = "arpack",
) -> Dict[float, pd.Series]:
    """
    Build a kNN graph on either the raw embedding dims or PCA scores and run Leiden.

    Defaults:
      - metric = "cosine" when cluster_space=="embedding"
      - metric = "euclidean" when cluster_space=="pca"
    """
    if metric is None:
        metric = "cosine" if cluster_space == "embedding" else "euclidean"

    adata = _build_representation(
        emb_df=emb_df,
        cluster_space=cluster_space,
        n_pcs=n_pcs,
        scale_pca=scale_pca,
        pca_svd_solver=pca_svd_solver,
    )

    if cluster_space == "pca":
        sc.pp.neighbors(adata, n_neighbors=n_neighbors, use_rep="X_pca", metric=metric)
    else:
        sc.pp.neighbors(adata, n_neighbors=n_neighbors, use_rep='X', metric=metric)

    out = {}
    for res in resolutions:
        key = f"leiden_{res:g}"
        sc.tl.leiden(adata, resolution=float(res), key_added=key, random_state=random_state)
        out[float(res)] = pd.Series(
            adata.obs[key].astype(str).values, index=adata.obs_names, name="cluster"
        )
    return out


# -----------------------------
# g:Profiler enrichment (parallelized)
# -----------------------------
def _worker_gprofiler(task: dict) -> pd.DataFrame:
    """
    One (resolution, cluster) enrichment job executed in a separate process.

    Returns a DataFrame with all significant rows for this job, including:
      - 'term' (GO term description)
      - 'go_id'
      - 'pvalue_adj' (adjusted p from g:Profiler)
      - 'pvalue' (raw p if present; else NaN)
      - 'term_size' (genes annotated to the GO term)
      - 'query_size' (genes in the cluster/query)
      - 'intersection_size' (overlap count)
      - 'intersection_genes' (overlap gene list as comma-separated string)
    """
    set_name = task["set_name"]
    res = float(task["resolution"])
    cl = str(task["cluster"])
    genes = task["genes"]
    go_sources = list(task["go_sources"])
    alpha = float(task["alpha"])
    organism = task["organism"]
    correction = task["correction"]
    domain_scope = task["domain_scope"]
    exclude_iea = bool(task["exclude_iea"])
    max_retries = int(task["max_retries"])
    retry_backoff = float(task["retry_backoff"])

    gp = GProfiler(return_dataframe=True)

    last_exc: Optional[Exception] = None
    enr = None
    for attempt in range(max_retries):
        try:
            enr = gp.profile(
                organism=organism,
                query=genes,
                sources=go_sources,
                user_threshold=alpha,
                significance_threshold_method=correction,
                domain_scope=domain_scope,
                no_evidences=exclude_iea,   # legacy naming, ok to leave as-is
                all_results=True,
            )
            last_exc = None
            break
        except Exception as e:
            last_exc = e
            time.sleep(retry_backoff * (2 ** attempt))

    if last_exc is not None or enr is None or len(enr) == 0:
        return pd.DataFrame(columns=[
            "set","resolution","cluster","library","term","go_id",
            "pvalue","pvalue_adj","term_size","query_size",
            "intersection_size","intersection_genes","term_key"
        ])

    df = enr.copy()

    # Robust column detection
    go_id_col = "native" if "native" in df.columns else ("term_id" if "term_id" in df.columns else None)
    name_col  = "name"   if "name"   in df.columns else ("description" if "description" in df.columns else None)
    p_adj_col = "p_value" if "p_value" in df.columns else None  # g:Profiler's p_value is adjusted
    # Try to find an unadjusted/raw p if g:Profiler exposes one in this environment
    p_raw_col_candidates = ["p_value_raw", "pval_raw", "p_value_unadjusted", "p_raw"]
    p_raw_col = next((c for c in p_raw_col_candidates if c in df.columns), None)

    # Intersection genes col (string with comma-separated genes in the overlap)
    inter_col_candidates = ["intersection", "intersections"]
    inter_col = next((c for c in inter_col_candidates if c in df.columns), None)

    if p_adj_col is None or go_id_col is None or name_col is None:
        return pd.DataFrame(columns=[
            "set","resolution","cluster","library","term","go_id",
            "pvalue","pvalue_adj","term_size","query_size",
            "intersection_size","intersection_genes","term_key"
        ])

    if "source" in df.columns:
        df = df[df["source"].isin(go_sources)].copy()

    # Keep terms passing the chosen alpha on the (adjusted) p-values
    df = df[np.isfinite(df[p_adj_col]) & (df[p_adj_col] <= alpha)].copy()
    if df.empty:
        return pd.DataFrame(columns=[
            "set","resolution","cluster","library","term","go_id",
            "pvalue","pvalue_adj","term_size","query_size",
            "intersection_size","intersection_genes","term_key"
        ])

    # Build output
    out = pd.DataFrame({
        "set": set_name,
        "resolution": res,
        "cluster": cl,
        "library": df.get("source", "GO").astype(str).values,
        "term": df[name_col].astype(str).values,   # description
        "go_id": df[go_id_col].astype(str).where(df[go_id_col].notna(), None),
        "pvalue_adj": df[p_adj_col].astype(float).values,
        "pvalue": pd.to_numeric(df[p_raw_col], errors="coerce").values if p_raw_col else np.nan,
        "term_size": pd.to_numeric(df.get("term_size"), errors="coerce"),
        "query_size": pd.to_numeric(df.get("query_size"), errors="coerce"),
        "intersection_size": pd.to_numeric(df.get("intersection_size"), errors="coerce"),
        "intersection_genes": (
            df[inter_col].astype(str).str.replace(";", ",") if inter_col in df.columns else ""
        ),
    })

    # Ensure integer types where applicable
    for col in ("term_size", "query_size", "intersection_size"):
        if col in out.columns:
            out[col] = pd.to_numeric(out[col], errors="coerce").astype("Int64")

    # Helper key for per-resolution uniqueness if needed downstream
    tk = out["go_id"].where(out["go_id"].notna(), out["library"].astype(str) + "::" + out["term"].astype(str))
    out["term_key"] = tk.astype(str)

    return out


def gprofiler_enrich_clusters_for_set(
    clusters_by_res: Dict[float, pd.Series],
    go_sources: Sequence[str],
    alpha: float = 0.05,
    min_cluster_size: int = 10,
    organism: str = "drerio",
    correction: str = "fdr",          # {"fdr","bonferroni","g_SCS"}
    domain_scope: str = "annotated",  # {"annotated","known"}
    exclude_iea: bool = False,        # kept for compatibility with original code
    set_name: str = "set",
    parallel: bool = True,
    max_workers: int = 4,
    max_retries: int = 3,
    retry_backoff: float = 1.5,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Run g:Profiler GO enrichment for each cluster at each resolution.

    Returns:
      details_df: rows = significant terms (p_value <= alpha under chosen correction)
      summary_df: counts per (set, resolution) — total/unique terms AND number of clusters
                  (column 'n_clusters' counts all clusters produced by Leiden at that resolution,
                   irrespective of min_cluster_size filtering for enrichment).
    """
    # Precompute number of clusters per resolution (no min size filtering)
    cluster_counts = {
        float(res): int(clusters.astype(str).nunique())
        for res, clusters in clusters_by_res.items()
    }

    # Build tasks (one per resolution × cluster that passes min_cluster_size)
    tasks: List[dict] = []
    for res, clusters in clusters_by_res.items():
        for cl in sorted(clusters.unique(), key=str):
            genes = clusters.index[clusters == cl].tolist()
            if len(genes) < min_cluster_size:
                continue
            tasks.append({
                "set_name": set_name,
                "resolution": float(res),
                "cluster": str(cl),
                "genes": genes,
                "go_sources": list(go_sources),
                "alpha": alpha,
                "organism": organism,
                "correction": correction,
                "domain_scope": domain_scope,
                "exclude_iea": exclude_iea,
                "max_retries": max_retries,
                "retry_backoff": retry_backoff,
            })

    if not tasks:
        empty_details = pd.DataFrame(columns=[
            "set", "resolution", "cluster", "library", "term", "go_id",
            "pvalue_adj", "term_size", "query_size", "intersection_size"
        ])
        # Include n_clusters even if no tasks (still have clusters from Leiden)
        summary_rows = []
        for res in sorted(cluster_counts.keys()):
            summary_rows.append({
                "set": set_name,
                "resolution": float(res),
                "n_sig_terms_total": 0,
                "n_sig_terms_unique": 0,
                "n_clusters": int(cluster_counts.get(float(res), 0)),
            })
        empty_summary = pd.DataFrame.from_records(summary_rows).sort_values(
            ["set", "resolution"], ignore_index=True
        )
        return empty_details, empty_summary

    # Execute tasks
    results: List[pd.DataFrame] = []
    if parallel and len(tasks) > 1:
        ctx = mp.get_context("spawn")
        with ProcessPoolExecutor(max_workers=max_workers, mp_context=ctx) as ex:
            futures = [ex.submit(_worker_gprofiler, t) for t in tasks]
            for fut in as_completed(futures):
                df = fut.result()
                if df is not None and not df.empty:
                    results.append(df)
    else:
        for t in tasks:
            df = _worker_gprofiler(t)
            if df is not None and not df.empty:
                results.append(df)

    if results:
        details_with_keys = pd.concat(results, axis=0, ignore_index=True)
    else:
        details_with_keys = pd.DataFrame(columns=[
            "set", "resolution", "cluster", "library", "term", "go_id",
            "pvalue_adj", "term_size", "query_size", "intersection_size", "term_key"
        ])

    # Build summary per resolution (include all resolutions present in clusters_by_res)
    summary_rows: List[dict] = []
    for res in sorted(float(r) for r in clusters_by_res.keys()):
        sub = details_with_keys[details_with_keys["resolution"] == float(res)]
        n_total = int(len(sub))
        n_unique = int(sub["term_key"].nunique()) if "term_key" in sub.columns else 0
        summary_rows.append({
            "set": set_name,
            "resolution": float(res),
            "n_sig_terms_total": n_total,
            "n_sig_terms_unique": n_unique,
            "n_clusters": int(cluster_counts.get(float(res), 0)),
        })

    # Finalize outputs: drop helper column
    details_df = details_with_keys.drop(columns=["term_key"], errors="ignore").sort_values(
        ["set", "resolution", "cluster", "library", "pvalue_adj"], ignore_index=True
    )
    summary_df = pd.DataFrame.from_records(summary_rows).sort_values(
        ["set", "resolution"], ignore_index=True
    )

    return details_df, summary_df

def _write_all_enrichments(details_df: pd.DataFrame, out_prefix: str, out_dir: str) -> Optional[str]:
    out_path = os.path.join(out_dir, f"{out_prefix}_enriched_terms.tsv")
    cols = [
        "set","resolution","cluster","library",
        "go_id","term","intersection_size","intersection_genes",
        "pvalue","pvalue_adj","term_size","query_size",
    ]
    if details_df is None or details_df.empty:
        pd.DataFrame(columns=cols).to_csv(out_path, index=False, sep="\t")
        print(f"[INFO] Wrote (empty): {out_path}")
        return out_path

    present = [c for c in cols if c in details_df.columns]
    df = (details_df[present]
          .sort_values(["set","resolution","cluster","library","pvalue_adj","term","go_id"],
                       kind="mergesort", na_position="last"))
    df.to_csv(out_path, index=False, sep="\t")
    print(f"[INFO] Wrote combined enriched terms: {out_path}")
    return out_path

# -----------------------------
# Plotting
# -----------------------------
from matplotlib import colors as mcolors

# --- Palette helper ----------------------------------------------------------
# Accepts: dict {label: color}, list of colors, or a matplotlib palette name.
# Falls back to a named palette (default "tab10") for any missing labels.
def _build_palette_for_labels(
    labels: List[str],
    palette_spec: Optional[object] = None,
    fallback_name: str = "tab10",
) -> Dict[str, str]:
    labels = list(labels)
    out: Dict[str, str] = {}

    # Build a fallback list
    fb = plt.colormaps.get(fallback_name)
    if fb is None:
        fb = plt.colormaps.get("tab10")
    fb_list = [fb(i / max(1, len(labels) - 1)) for i in range(len(labels))]

    def _maybe_base(label: str) -> str:
        # convenience: allow mapping just "Set" for a label like "Set [PCA50]"
        if " [" in label:
            return label.split(" [", 1)[0]
        return label

    if isinstance(palette_spec, dict):
        # dict: exact (or base) label matches first; missing → fallback cycle
        for i, lab in enumerate(labels):
            if lab in palette_spec:
                out[lab] = palette_spec[lab]
            else:
                base = _maybe_base(lab)
                out[lab] = palette_spec.get(base, fb_list[i])
        return out

    if isinstance(palette_spec, (list, tuple)):
        # list: assign in order; cycle if fewer colors than labels
        pal_list = list(palette_spec)
        for i, lab in enumerate(labels):
            out[lab] = pal_list[i % len(pal_list)]
        return out

    if isinstance(palette_spec, str):
        # string: treat as a matplotlib palette name
        cm = plt.colormaps.get(palette_spec)
        if cm is not None:
            pal_list = [cm(i / max(1, len(labels) - 1)) for i in range(len(labels))]
            for i, lab in enumerate(labels):
                out[lab] = pal_list[i]
            return out
        # otherwise, treat as a single color string (rare, but safe)
        for lab in labels:
            out[lab] = palette_spec
        return out

    # default fallback
    for i, lab in enumerate(labels):
        out[lab] = fb_list[i]
    return out


def _save_dual_formats(fig, out_path: Optional[str]):
    """
    Save figure to both .svg and .png. If out_path has an extension, it's stripped to a base.
    """
    if not out_path:
        return None, None
    base, ext = os.path.splitext(out_path)
    base = base if ext.lower() in (".svg", ".png") else out_path
    svg_path = base + ".svg"
    png_path = base + ".png"
    fig.savefig(svg_path, bbox_inches="tight")
    fig.savefig(png_path, dpi=150, bbox_inches="tight")
    return svg_path, png_path

def _coerce_summary_df(df: pd.DataFrame) -> pd.DataFrame:
    required_cols = {"set", "resolution", "n_sig_terms_unique"}
    missing = required_cols.difference(df.columns)
    if missing:
        raise ValueError(f"summary_df is missing required columns: {sorted(missing)}")

    df = df.copy()
    df["resolution"] = pd.to_numeric(df["resolution"], errors="coerce")
    df["n_sig_terms_unique"] = pd.to_numeric(df["n_sig_terms_unique"], errors="coerce")
    if "n_clusters" in df.columns:
        df["n_clusters"] = pd.to_numeric(df["n_clusters"], errors="coerce")

    df = df[np.isfinite(df["resolution"]) & np.isfinite(df["n_sig_terms_unique"])].copy()
    df["set"] = df.get("set", "series").astype(str)
    return df

def _legend_no_markers(ax, title: str):
    """Place legend outside on the right, with handles that have NO markers."""
    handles, labels = ax.get_legend_handles_labels()
    if not handles:
        return
    # Create proxy handles without markers, matching line color/linestyle/linewidth
    proxy = [
        Line2D([0], [0],
               color=h.get_color(),
               linestyle=h.get_linestyle() if hasattr(h, "get_linestyle") else "-",
               linewidth=h.get_linewidth() if hasattr(h, "get_linewidth") else 1.5,
               marker=None)
        for h in handles
    ]
    ax.legend(proxy, labels, title=title, frameon=False,
              loc="center left", bbox_to_anchor=(1.15, 0.5), borderaxespad=0)


def plot_summary_from_df(
    summary_df: pd.DataFrame,
    out_base: Optional[str] = None,
    title_suffix: Optional[str] = None,
    right_margin: float = 0.75,
    custom_palette: Optional[object] = None,   # <-- NEW
    include_sets: Optional[Sequence[str]] = None,  # <-- NEW: specify series to plot
) -> Tuple[Optional[str], Optional[str]]:
    """
    Dual y-axes plot: left = unique GO terms (solid + filled circles), right = # clusters (dashed + open circles).
    Saves both .svg and .png using out_base as the filename base (no extension needed).
    Legend is outside and rendered without markers.
    """
    df = _coerce_summary_df(summary_df)
    if df.empty:
        print("WARNING: empty summary_df; nothing to plot.", file=sys.stderr)
        return None, None

    if include_sets:
        requested = [s for s in include_sets if s in df["set"].unique()]
        if not requested:
            print("WARNING: none of the requested series present in summary_df; nothing to plot.", file=sys.stderr)
            return None, None
        df = df[df["set"].isin(requested)]
        labels: List[str] = requested
    else:
        labels: List[str] = df["set"].unique().tolist()

    pal = _build_palette_for_labels(labels, custom_palette)

    fig, ax = plt.subplots(figsize=(8.6, 5.2))
    fig.subplots_adjust(right=right_margin)

    has_clusters = "n_clusters" in df.columns and df["n_clusters"].notna().any()
    ax2 = ax.twinx() if has_clusters else None

    for label in labels:
        sub = df[df["set"] == label].sort_values("resolution")

        line, = ax.plot(
            sub["resolution"].values,
            sub["n_sig_terms_unique"].values,
            marker="o",
            label=label,
            color=pal[label],                  # <-- use manual color
        )

        if has_clusters and ax2 is not None:
            ax2.plot(
                sub["resolution"].values,
                sub["n_clusters"].values,
                linestyle="--",
                linewidth=1.0,
                alpha=0.5,
                marker="o",
                markerfacecolor="none",
                markeredgecolor=pal[label],
                markeredgewidth=1.0,
                color=pal[label],
            )

    ax.set_xlabel("Leiden resolution")
    ax.set_ylabel("Unique significant GO terms (q ≤ threshold)\n—●—")
    if has_clusters and ax2 is not None:
        ax2.set_ylabel("Number of clusters\n– – ○ – – ")

    title_main = "GO enrichment vs. clustering resolution (g:Profiler)"
    if title_suffix:
        title_main += f"\n{title_suffix}"
    ax.set_title(title_main)
    ax.grid(True, linestyle="--", alpha=0.4)

    _legend_no_markers(ax, title="Series")

    paths = _save_dual_formats(fig, out_base)
    plt.close(fig)
    return paths

# -----------------------------
# umap
# -----------------------------
def _sanitize_filename(s: str) -> str:
    return re.sub(r"[^A-Za-z0-9._-]+", "_", str(s)).strip("_")

def compute_umap_embedding(
    emb_df: pd.DataFrame,
    cluster_space: str = "embedding",
    n_pcs: int = 0,
    scale_pca: bool = False,
    n_neighbors: int = 20,
    metric_umap: str = "cosine",
    random_state: int = 42,
    pca_svd_solver: str = "arpack",
) -> pd.DataFrame:
    """
    Build representation (embedding or PCA), compute neighbors with COSINE, then UMAP.
    Returns a DataFrame with index=genes, columns=['UMAP1','UMAP2'].
    """
    adata = _build_representation(
        emb_df=emb_df,
        cluster_space=cluster_space,
        n_pcs=n_pcs,
        scale_pca=scale_pca,
        pca_svd_solver=pca_svd_solver,
    )
    use_rep = "X_pca" if cluster_space == "pca" else "X"
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, use_rep=use_rep, metric=metric_umap)
    sc.tl.umap(adata, random_state=random_state, n_components=2, min_dist=1.5, spread=3.0, gamma=5, maxiter=3000, negative_sample_rate=20)
    coords = np.asarray(adata.obsm["X_umap"], dtype=float)
    return pd.DataFrame(coords, index=adata.obs_names, columns=["UMAP1", "UMAP2"])


def _combined_categorical_palette() -> list[str]:
    """Concatenate several Matplotlib ListedColormaps into one long palette."""
    names = ['tab20b', 'tab20c', 'Set1', 'Set2', 'Set3', 'Pastel1', 'Pastel2', 'Accent', 'Dark2']
    colors = []
    for name in names:
        try:
            cmap = plt.get_cmap(name)
            listed = getattr(cmap, 'colors', None)
            if listed is None:
                continue  # skip non-listed cmaps just in case
            colors.extend(mcolors.to_hex(c) for c in listed)
        except ValueError:
            # palette not available in this Matplotlib build; skip it
            continue
    # de-duplicate while preserving order
    seen = set()
    unique = []
    for c in colors:
        if c not in seen:
            seen.add(c)
            unique.append(c)
    #randomly shuffle the colors
    import random
    random.shuffle(unique)
    return unique

def plot_umap_clusters(
    umap_df: pd.DataFrame,
    clusters: pd.Series,
    title: str,
    out_base: Optional[str] = None,
) -> Tuple[Optional[str], Optional[str]]:
    """
    Scatter UMAP colored by cluster labels using a concatenated categorical palette.
    Saves both .svg and .png if out_base is given. (No legend, no grid.)
    """
    df = umap_df.copy()
    cl = clusters.astype(str).reindex(df.index)
    df["cluster"] = cl.values
    labels = sorted(df["cluster"].dropna().unique().tolist(), key=str)

    # Build the combined palette and assign colors by cycling if needed
    base_colors = _combined_categorical_palette()
    if not base_colors:
        raise RuntimeError("No categorical palettes found. Check your Matplotlib installation.")
    color_cycle = list(islice(cycle(base_colors), len(labels)))
    pal = dict(zip(labels, color_cycle))

    fig, ax = plt.subplots(figsize=(7.0, 6.0))
    for lab in labels:
        sub = df[df["cluster"] == lab]
        ax.scatter(sub["UMAP1"], sub["UMAP2"], s=8, alpha=0.8, color=pal[lab])

    ax.set_xlabel("UMAP1")
    ax.set_ylabel("UMAP2")
    ax.set_title(title)
    ax.grid(False)  # ensure no grid even if a style sets one

    paths = _save_dual_formats(fig, out_base)
    plt.close(fig)
    return paths


def _write_all_umap_embeddings(all_umap_rows: List[pd.DataFrame], out_prefix: str, out_dir: str) -> Optional[str]:
    """
    Concatenate and write all UMAP rows to one CSV.
    Columns: set, resolution, gene, UMAP1, UMAP2, cluster
    """
    out_path = os.path.join(out_dir, f"{out_prefix}_umap_embeddings.csv")
    if not all_umap_rows:
        pd.DataFrame(columns=["set","resolution","gene","UMAP1","UMAP2","cluster"]).to_csv(out_path, index=False)
        print(f"[INFO] Wrote (empty UMAP table): {out_path}")
        return out_path
    df = pd.concat(all_umap_rows, axis=0, ignore_index=True)
    df.to_csv(out_path, index=False)
    print(f"[INFO] Wrote UMAP embeddings: {out_path}")
    return out_path

# -----------------------------
# Co-expression (precomputed clusters) loader
# -----------------------------
_COEXP_FILE_RE = re.compile(r"^leiden_clusters_([0-9]*\.?[0-9]+)\.csv$", re.IGNORECASE)

def load_coexpression_clusters_from_dir(
    dir_path: str,
    expected_resolutions: Optional[Sequence[float]] = None,
) -> Dict[float, pd.Series]:
    """
    Read precomputed Leiden clusters from CSVs in dir_path.
    Each file must be named like: leiden_clusters_<resolution>.csv
    Example header:
        ,leiden
        ENSDARG00000000001,12
        ...

    Returns:
        Dict[float, pd.Series] mapping resolution -> cluster series (index=genes, values=str cluster IDs)
    """
    if not dir_path or not os.path.isdir(dir_path):
        raise ValueError(f"Co-expression directory not found: {dir_path}")

    # Normalize expected resolutions to a set of float (if provided)
    want: Optional[set] = None
    if expected_resolutions is not None:
        want = set(float(x) for x in expected_resolutions)

    clusters_by_res: Dict[float, pd.Series] = {}
    seen_res = set()

    for fname in os.listdir(dir_path):
        m = _COEXP_FILE_RE.match(fname)
        if not m:
            continue
        res_str = m.group(1)
        try:
            res = float(res_str)
        except Exception:
            continue
        if want is not None and res not in want:
            # Skip files not requested by --resolutions
            continue

        path = os.path.join(dir_path, fname)
        try:
            df = pd.read_csv(path, index_col=0)
        except Exception as e:
            print(f"[WARN] Failed reading {path}: {e}", file=sys.stderr)
            continue

        if df.shape[1] == 0:
            print(f"[WARN] No columns in {path}; skipping.", file=sys.stderr)
            continue

        # Prefer a column explicitly named 'leiden'; else fall back to the first column.
        col = "leiden" if "leiden" in df.columns else df.columns[0]
        cl = df[col].copy()

        # Clean: ensure string clusters, string index; drop NAs.
        cl = cl.dropna()
        cl.index = cl.index.astype(str)
        cl = cl.astype(str)

        # Guard duplicates: if index not unique, keep the first occurrence.
        if not cl.index.is_unique:
            cl = cl[~cl.index.duplicated(keep="first")]

        clusters_by_res[res] = pd.Series(cl.values, index=cl.index, name="cluster")
        seen_res.add(res)

    if want is not None:
        missing = sorted(want.difference(seen_res))
        if missing:
            print(f"[WARN] Missing co-expression cluster files for resolutions: {missing}", file=sys.stderr)

    return clusters_by_res


# -----------------------------
# Orchestration across sets and spaces
# -----------------------------
def run_pipeline(
    embedding_sets: Mapping[str, pd.DataFrame],
    resolutions: Sequence[float],
    go_sources: Sequence[str] = DEFAULT_GO_SOURCES,
    alpha: float = 0.05,
    min_cluster_size: int = 10,
    n_neighbors: int = 20,
    random_state: int = 0,
    align_by_intersection: bool = True,
    out_prefix: str = "enrichment",
    organism: str = "drerio",
    correction: str = "fdr",
    domain_scope: str = "annotated",
    exclude_iea: bool = False,
    # Clustering space controls
    cluster_spaces: Sequence[str] = ("embedding",),
    metric_embedding: str = "cosine",
    metric_pca: str = "euclidean",
    n_pcs: int = 0,
    scale_pca: bool = False,
    pca_svd_solver: str = "arpack",
    # Explained variance plot controls
    make_pca_ev_plot: bool = True,
    ev_max_pcs: int = 100,
    ev_scale_pca: Optional[bool] = None,
    ev_mode: str = "cumulative",
    # Parallel enrichment controls
    parallel_gprofiler: bool = True,
    max_workers: int = 4,
    max_retries: int = 3,
    retry_backoff: float = 1.5,
    out_dir: str = "enrich_results",
    umap_out_dir: str = "umap_results",
    custom_palette: Optional[object] = None,
    coexp_cluster_dir: Optional[str] = None,
    plot_series: Optional[Sequence[str]] = None,  # <-- NEW: specify series to plot
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Full pipeline: clustering (raw and/or PCA) + GO enrichment with g:Profiler.

    cluster_spaces:
      - ("embedding",)        → only raw embedding dims
      - ("pca",)              → only PCA
      - ("embedding","pca")   → both; results labeled for side-by-side comparison

    Returns:
      details_df, summary_df
    """
    os.makedirs(out_dir, exist_ok=True)
    os.makedirs(umap_out_dir, exist_ok=True)
    # Validate sets and make indices unique
    sets = {k: v.copy() for k, v in embedding_sets.items()}
    for k, df in sets.items():
        if not isinstance(df.index, pd.Index):
            raise ValueError(f"Embedding set '{k}' must be a DataFrame with gene IDs as the index.")
        if not df.index.is_unique:
            sets[k] = df.groupby(df.index).mean()

    # Align by gene intersection across sets (fairness)
    if align_by_intersection:
        common = None
        for df in sets.values():
            common = df.index if common is None else common.intersection(df.index)
        for k in list(sets.keys()):
            sets[k] = sets[k].loc[common]
        print(f"[INFO] Using intersection of genes across sets: {len(common)} genes.")

    # PCA explained variance across sets
    if make_pca_ev_plot:
        ev_scale = scale_pca if ev_scale_pca is None else bool(ev_scale_pca)
        ev_df = compute_pca_explained_variance_across_sets(
            embedding_sets=sets,
            max_pcs=ev_max_pcs,
            scale=ev_scale,
            pca_svd_solver=pca_svd_solver,
        )
        ev_csv = os.path.join(out_dir, f"{out_prefix}_pca_explained_variance.csv")
        ev_png = os.path.join(out_dir, f"{out_prefix}_pca_explained_variance.png")
        ev_df.to_csv(ev_csv, index=False)

        vline = n_pcs if ("pca" in set(map(str, cluster_spaces))) else None
        plot_pca_explained_variance(
            ev_df,
            out_png=ev_png,
            mode=ev_mode,
            vline_at=vline,
            title="PCA explained variance across embedding sets",
            custom_palette=custom_palette,
            include_sets=plot_series,  # <-- NEW
        )
        print(f"[INFO] Wrote: {ev_csv}, {ev_png}")

    all_details = []
    all_summary = []
    all_umap_tables: List[pd.DataFrame] = []

    # Run per set × cluster_space
    for set_name, emb_df in sets.items():
        for space in cluster_spaces:
            label = f"{set_name}" if space == "embedding" else f"{set_name} [PCA{n_pcs}]"

            print(f"[INFO] Clustering '{set_name}' in space='{space}' ...")
            clusters_by_res = leiden_clusters_from_embeddings(
                emb_df=emb_df,
                resolutions=resolutions,
                cluster_space=space,
                n_pcs=n_pcs,
                scale_pca=scale_pca,
                n_neighbors=n_neighbors,
                metric=(metric_embedding if space == "embedding" else metric_pca),
                random_state=random_state,
                pca_svd_solver=pca_svd_solver,
            )

            print(f"[INFO] Enrichment (g:Profiler) for '{label}', GO sources: {go_sources} "
                  f"({'parallel' if parallel_gprofiler else 'serial'}, workers={max_workers}) ...")
            details_df, summary_df = gprofiler_enrich_clusters_for_set(
                clusters_by_res=clusters_by_res,
                go_sources=go_sources,
                alpha=alpha,
                min_cluster_size=min_cluster_size,
                organism=organism,
                correction=correction,
                domain_scope=domain_scope,
                exclude_iea=exclude_iea,
                set_name=label,   # keep label with space identifier
                parallel=parallel_gprofiler,
                max_workers=max_workers,
                max_retries=max_retries,
                retry_backoff=retry_backoff,
            )
            all_details.append(details_df)
            all_summary.append(summary_df)

            # --- NEW: Cosine-based UMAP per set×space, with a plot for each resolution ---
            print(f"[INFO] UMAP (cosine) for '{label}' ...")
            umap_coords = compute_umap_embedding(
                emb_df=emb_df,
                cluster_space="embedding",
                n_pcs=0, # use all PCs
                scale_pca=False,
                n_neighbors=n_neighbors,
                metric_umap="cosine",
                random_state=random_state,
                pca_svd_solver=pca_svd_solver,
            )

            # For each resolution: plot same UMAP coords, colored by clusters at that resolution,
            # and add rows (duplicated coords per resolution) to the combined CSV.
            for res, cl_ser in sorted(clusters_by_res.items(), key=lambda x: float(x[0])):
                safe_label = _sanitize_filename(label)
                out_base = os.path.join(umap_out_dir, f"{out_prefix}_{safe_label}_res{float(res):g}_umap")
                plot_umap_clusters(
                    umap_df=umap_coords,
                    clusters=cl_ser,
                    title=f"UMAP (cosine) — {label} — res={float(res):g}",
                    out_base=out_base,
                )
                # gather rows for the combined CSV
                df_rows = umap_coords.copy()
                df_rows["set"] = label
                df_rows["resolution"] = float(res)
                df_rows["gene"] = df_rows.index.astype(str)
                df_rows["cluster"] = cl_ser.reindex(df_rows.index).astype(str).values
                all_umap_tables.append(df_rows.reset_index(drop=True)[["set","resolution","gene","UMAP1","UMAP2","cluster"]])

    # --- NEW: Precomputed co-expression set (clusters only; no embeddings/UMAP/EV) ---
    if coexp_cluster_dir:
        print(f"[INFO] Loading precomputed clusters for 'co-expression' from: {coexp_cluster_dir}")
        try:
            coexp_clusters_by_res = load_coexpression_clusters_from_dir(
                coexp_cluster_dir, expected_resolutions=resolutions
            )
        except Exception as e:
            coexp_clusters_by_res = {}
            print(f"[WARN] Could not load co-expression clusters: {e}", file=sys.stderr)

        if not coexp_clusters_by_res:
            print("[WARN] No valid co-expression cluster files found; skipping co-expression enrichment.", file=sys.stderr)
        else:
            print(f"[INFO] Enrichment (g:Profiler) for 'co-expression', GO sources: {go_sources} "
                f"({'parallel' if parallel_gprofiler else 'serial'}, workers={max_workers}) ...")
            details_df, summary_df = gprofiler_enrich_clusters_for_set(
                clusters_by_res=coexp_clusters_by_res,
                go_sources=go_sources,
                alpha=alpha,
                min_cluster_size=min_cluster_size,
                organism=organism,
                correction=correction,
                domain_scope=domain_scope,
                exclude_iea=exclude_iea,
                set_name="co-expression",
                parallel=parallel_gprofiler,
                max_workers=max_workers,
                max_retries=max_retries,
                retry_backoff=retry_backoff,
            )
            all_details.append(details_df)
            all_summary.append(summary_df)
           


    details = pd.concat(all_details, axis=0, ignore_index=True) if len(all_details) else pd.DataFrame()
    summary = pd.concat(all_summary, axis=0, ignore_index=True) if len(all_summary) else pd.DataFrame()

    # Save combined tables
    details_out = os.path.join(out_dir, f"{out_prefix}_details.csv")
    summary_out = os.path.join(out_dir, f"{out_prefix}_summary.csv")
    details.to_csv(details_out, index=False)
    summary.to_csv(summary_out, index=False)

    # One TSV for all resolutions with requested fields
    _write_all_enrichments(details, out_prefix=out_prefix, out_dir=out_dir)

    # Plots into the folder
    plot_base = os.path.join(out_dir, f"{out_prefix}_summary_plot")
    plot_summary_from_df(summary, out_base=plot_base, title_suffix=out_prefix, custom_palette=custom_palette, include_sets=plot_series)

    # NEW: Save all UMAP embeddings into one file
    _write_all_umap_embeddings(all_umap_tables, out_prefix=out_prefix, out_dir=out_dir)

    return details, summary


# -----------------------------
# CLI
# -----------------------------
def _parse_kv_embedding(arg: str) -> Tuple[str, str]:
    """
    Parse "name=/path/to.csv" into (name, path).
    """
    if "=" not in arg:
        raise argparse.ArgumentTypeError(f"Embedding must be NAME=PATH, got: {arg}")
    name, path = arg.split("=", 1)
    if not name:
        raise argparse.ArgumentTypeError("Embedding NAME cannot be empty.")
    if not os.path.exists(path):
        raise argparse.ArgumentTypeError(f"Embedding path does not exist: {path}")
    return name, path


def _parse_float_list(arg: str) -> List[float]:
    try:
        vals = [float(x.strip()) for x in arg.split(",") if x.strip() != ""]
    except Exception:
        raise argparse.ArgumentTypeError(f"Could not parse float list: {arg}")
    if not vals:
        raise argparse.ArgumentTypeError("At least one resolution required.")
    return vals


def _parse_str_list(arg: str) -> List[str]:
    vals = [x.strip() for x in arg.split(",") if x.strip() != ""]
    if not vals:
        raise argparse.ArgumentTypeError("List must contain at least one item.")
    return vals


def build_embedding_sets_from_args(pairs: List[Tuple[str, str]]) -> Dict[str, pd.DataFrame]:
    embedding_sets: Dict[str, pd.DataFrame] = {}
    for name, path in pairs:
        df = pd.read_csv(path, index_col=0)
        if not isinstance(df.index, pd.Index):
            raise ValueError(f"CSV {path} must have gene IDs as row index.")
        embedding_sets[name] = df
    return embedding_sets


def main():
    parser = argparse.ArgumentParser(
        description="Cluster gene embeddings and run GO enrichment via g:Profiler (parallelized).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--embedding", "-e", action="append", required=False, type=_parse_kv_embedding,
        help="Embedding CSV as NAME=PATH (repeat for multiple). Index must be gene IDs."
    )
    parser.add_argument(
        "--resolutions", "-r", type=_parse_float_list, required=False,
        help="Comma-separated Leiden resolutions, e.g., 0.2,0.4,0.6,0.8,1.0"
    )
    parser.add_argument(
        "--cluster-spaces", type=_parse_str_list, default=["embedding", "pca"],
        help="Comma-separated among: embedding,pca"
    )
    parser.add_argument("--n-pcs", type=int, default=225, help="#PCs for PCA space.")
    parser.add_argument("--scale-pca", action="store_true", help="Standardize features before PCA.")
    parser.add_argument("--no-scale-pca", dest="scale_pca", action="store_false")
    parser.set_defaults(scale_pca=False)
    parser.add_argument("--pca-svd-solver", type=str, default="arpack", choices=["arpack", "randomized"])

    parser.add_argument("--n-neighbors", type=int, default=20)
    parser.add_argument("--metric-embedding", type=str, default="cosine")
    parser.add_argument("--metric-pca", type=str, default="euclidean")
    parser.add_argument("--random-state", type=int, default=0)
    parser.add_argument("--align-by-intersection", action="store_true", default=True)
    parser.add_argument("--no-align-by-intersection", dest="align_by_intersection", action="store_false")

    parser.add_argument("--go-sources", type=_parse_str_list, default=DEFAULT_GO_SOURCES,
                        help="Comma-separated GO libraries, e.g., GO:BP,GO:MF,GO:CC,REAC")
    parser.add_argument("--alpha", type=float, default=0.05)
    parser.add_argument("--min-cluster-size", type=int, default=10)
    parser.add_argument("--organism", type=str, default="drerio")
    parser.add_argument("--correction", type=str, default="fdr", choices=["fdr", "bonferroni", "g_SCS"])
    parser.add_argument("--domain-scope", type=str, default="annotated", choices=["annotated", "known"])
    parser.add_argument("--exclude-iea", action="store_true", default=False,
                        help="Kept for compatibility with original code; mapped to g:Profiler 'no_evidences' flag.")

    parser.add_argument("--out-prefix", type=str, default="enrichment")

    # Explained variance (EV) plot
    parser.add_argument("--make-ev-plot", action="store_true", default=True)
    parser.add_argument("--no-ev-plot", dest="make_ev_plot", action="store_false")
    parser.add_argument("--ev-max-pcs", type=int, default=250)
    parser.add_argument("--ev-mode", type=str, default="cumulative", choices=["cumulative", "per_pc"])
    parser.add_argument("--ev-scale-pca", type=str, default=None,
                        help="Override EV scaling: 'true', 'false', or omit to inherit --scale-pca.")

    # Parallel enrichment controls
    parser.add_argument("--parallel", action="store_true", default=True)
    parser.add_argument("--no-parallel", dest="parallel", action="store_false")
    parser.add_argument("--max-workers", type=int, default=min(4, max(mp.cpu_count(), 1)))
    parser.add_argument("--max-retries", type=int, default=3)
    parser.add_argument("--retry-backoff", type=float, default=1.5)

    # Color palette controls
    parser.add_argument("--palette-name", type=str, default=None,
                        help="Matplotlib palette name (e.g., tab10, Set2, viridis).")
    parser.add_argument("--palette-list", type=_parse_str_list, default=None,
                        help="Comma-separated colors (hex or names), used in series order.")
    parser.add_argument("--palette-map", type=str, default=None,
                        help='Comma-separated NAME=COLOR pairs (match "set" or its base before " [").')
    parser.add_argument("--palette-file", type=str, default=None,
                        help="Path to JSON file containing {\"Series Name\": \"#hex\"} mapping.")

    # co-expression cluster directory
    parser.add_argument(
        "--coexp-dir", type=str, default=None,
        help="Directory containing precomputed Leiden cluster CSVs named like 'leiden_clusters_<res>.csv'. Adds a set named 'co-expression'."
    )
    
    # mock run
    parser.add_argument("--mock-summary-csv", default=None,
                        help="If provided, skip pipeline and just plot the summary from this CSV, then exit.")
    parser.add_argument("--mock-umap-csv", default=None,
                        help="If provided, skip pipeline and just plot the UMAP from this CSV, then exit.")
    parser.add_argument("--mock-title-suffix", default=None,
                        help="Optional second-line title suffix for mock summary plot.")

    # output directory
    parser.add_argument(
        "--out-dir", type=str, default="enrich_results",
        help="Directory to write outputs (created if missing)."
    )
    parser.add_argument(
        "--umap-out-dir", type=str, default="umap_results",
        help="Directory to write UMAP outputs (created if missing)."
    )

    # NEW: specify which series (sets) to include in plots
    parser.add_argument(
        "--plot-series", type=_parse_str_list, default=None,
        help="Comma-separated series names (matching the 'set' labels) to include in EV and summary plots. If omitted, plot all."
    )

    args = parser.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)

    def _parse_color_map(arg: str) -> Dict[str, str]:
        # "A=#1f77b4,B=#ff7f0e" → {"A":"#1f77b4","B":"#ff7f0e"}
        out = {}
        for pair in arg.split(","):
            pair = pair.strip()
            if not pair:
                continue
            if "=" not in pair:
                raise argparse.ArgumentTypeError(f"Bad NAME=COLOR pair: {pair}")
            k, v = pair.split("=", 1)
            out[k.strip()] = v.strip()
        return out

    # Palette resolution (priority: file > map > list > name)
    custom_palette = None
    if args.palette_file:
        import json
        with open(args.palette_file, "r") as fh:
            loaded = json.load(fh)
        if not isinstance(loaded, dict):
            raise ValueError("--palette-file must contain a JSON object mapping names to colors.")
        custom_palette = loaded
    elif args.palette_map:
        custom_palette = _parse_color_map(args.palette_map)
    elif args.palette_list:
        custom_palette = args.palette_list
    elif args.palette_name:
        custom_palette = args.palette_name

    # Interpret ev_scale_pca override
    if args.ev_scale_pca is None:
        ev_scale_pca = None
    else:
        s = str(args.ev_scale_pca).strip().lower()
        if s in ("true", "t", "1", "yes", "y"):
            ev_scale_pca = True
        elif s in ("false", "f", "0", "no", "n"):
            ev_scale_pca = False
        else:
            raise ValueError("--ev-scale-pca must be 'true' or 'false' if provided.")

    # --- Mock mode: just plot the summary from a CSV and exit ---
    if args.mock_summary_csv:
        df = pd.read_csv(args.mock_summary_csv)
        plot_base = os.path.join(args.out_dir, f"{args.out_prefix}_summary_plot")
        plot_summary_from_df(
            df,
            out_base=plot_base,
            title_suffix=(args.mock_title_suffix if args.mock_title_suffix else args.out_prefix),
            custom_palette=custom_palette,
            include_sets=args.plot_series,  # <-- NEW
        )
        print(f"[INFO] Mock run: wrote {plot_base}.svg and {plot_base}.png from {args.mock_summary_csv}")
        
        # NEW: UMAP mock plots from precomputed embeddings
        umap_csv_default = args.mock_umap_csv
        umap_csv_prefixed = os.path.join(args.out_dir, f"{args.out_prefix}_umap_embeddings.csv")
        umap_csv = None
        if os.path.exists(umap_csv_default):
            umap_csv = umap_csv_default
        elif os.path.exists(umap_csv_prefixed):
            umap_csv = umap_csv_prefixed

        if umap_csv is None:
            print(f"[INFO] Mock run: no UMAP embeddings file found at {umap_csv_default} or {umap_csv_prefixed}; skipping UMAP plots.")
            return

        umap_tbl = pd.read_csv(umap_csv)
        required_cols = {"set", "resolution", "gene", "UMAP1", "UMAP2", "cluster"}
        missing = required_cols.difference(umap_tbl.columns)
        if missing:
            print(f"[WARN] Mock run: UMAP embeddings missing columns {sorted(missing)}; skipping UMAP plots.")
            return

        # Clean types and NA
        umap_tbl["resolution"] = pd.to_numeric(umap_tbl["resolution"], errors="coerce")
        umap_tbl = umap_tbl.dropna(subset=["resolution", "UMAP1", "UMAP2", "gene", "set", "cluster"])

        # For each set×resolution, reuse the existing scatter helper
        for (set_label, res), group in umap_tbl.groupby(["set", "resolution"]):
            umap_df = group.set_index("gene")[["UMAP1", "UMAP2"]]
            clusters = group.set_index("gene")["cluster"]
            safe_label = _sanitize_filename(str(set_label))
            out_base = os.path.join(args.umap_out_dir, f"{args.out_prefix}_{safe_label}_res{float(res):g}_umap")
            plot_umap_clusters(
                umap_df=umap_df,
                clusters=clusters,
                title=f"UMAP (mock) — {set_label} — res={float(res):g}",
                out_base=out_base,
            )
        print(f"[INFO] Mock run: wrote UMAP plots from {umap_csv}")
        return

    # If not mock mode, enforce required inputs:
    # - Require at least embeddings or coexp-dir
    # - Require resolutions for both modes
    if (not args.embedding and not args.coexp_dir) or not args.resolutions:
        parser.error("At least one of --embedding or --coexp-dir, and also --resolutions, are required (unless --mock-summary-csv is provided).")

    # Build input sets
    embedding_sets = build_embedding_sets_from_args(args.embedding)
    resolutions = args.resolutions
    cluster_spaces = args.cluster_spaces

    # Ensure safe start method for multiprocessing across platforms
    try:
        mp.set_start_method("spawn", force=False)
    except RuntimeError:
        # start method already set; OK
        pass

    # Run
    details, summary = run_pipeline(
        embedding_sets=embedding_sets,
        resolutions=resolutions,
        go_sources=args.go_sources,
        alpha=args.alpha,
        min_cluster_size=args.min_cluster_size,
        n_neighbors=args.n_neighbors,
        random_state=args.random_state,
        align_by_intersection=args.align_by_intersection,
        out_prefix=args.out_prefix,
        organism=args.organism,
        correction=args.correction,
        domain_scope=args.domain_scope,
        exclude_iea=args.exclude_iea,
        cluster_spaces=cluster_spaces,
        metric_embedding=args.metric_embedding,
        metric_pca=args.metric_pca,
        n_pcs=args.n_pcs,
        scale_pca=args.scale_pca,
        pca_svd_solver=args.pca_svd_solver,
        make_pca_ev_plot=args.make_ev_plot,
        ev_max_pcs=args.ev_max_pcs,
        ev_scale_pca=ev_scale_pca,
        ev_mode=args.ev_mode,
        parallel_gprofiler=args.parallel,
        max_workers=args.max_workers,
        max_retries=args.max_retries,
        retry_backoff=args.retry_backoff,
        out_dir=args.out_dir,  
        umap_out_dir=args.umap_out_dir,
        custom_palette=custom_palette,
        coexp_cluster_dir=args.coexp_dir,
        plot_series=args.plot_series,  # <-- NEW
    )

    print(f"[INFO] Wrote: "
        f"{os.path.join(args.out_dir, f'{args.out_prefix}_details.csv')}, "
        f"{os.path.join(args.out_dir, f'{args.out_prefix}_summary.csv')}, "
        f"{os.path.join(args.out_dir, f'{args.out_prefix}_summary_plot.png')}")
    if args.make_ev_plot:
        print(f"[INFO] Wrote: "
            f"{os.path.join(args.out_dir, f'{args.out_prefix}_pca_explained_variance.csv')}, "
            f"{os.path.join(args.out_dir, f'{args.out_prefix}_pca_explained_variance.png')}")


if __name__ == "__main__":
    main()
