# GO enrichment analysis for Geneformer embeddings

from __future__ import annotations
import warnings, os
from typing import Dict, Mapping, Sequence, Tuple, Optional, List
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scanpy as sc
# g:Profiler official Python client
from gprofiler import GProfiler

# -----------------------------
# Config
# -----------------------------
DEFAULT_GO_SOURCES = ["GO:BP", "GO:MF", "GO:CC"]  # biological process, molecular function, cellular component


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
    mode: str = "cumulative",     # "cumulative" or "per_pc"
    vline_at: Optional[int] = None,  # e.g., n_pcs used for PCA clustering
    title: str = "PCA explained variance across embedding sets",
):
    """
    Makes a single-panel line plot (no seaborn, no explicit colors).
    """
    if ev_df is None or ev_df.empty:
        warnings.warn("No explained-variance data to plot.")
        return

    plt.figure(figsize=(8.0, 5.0))
    for set_name in ev_df["set"].unique():
        sub = ev_df[ev_df["set"] == set_name].sort_values("pc")
        y = sub["cumulative_variance_ratio"] if mode == "cumulative" else sub["variance_ratio"]
        plt.plot(sub["pc"].values, y.values, marker="o", label=set_name)

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
    n_pcs: int = 50,
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
        # Optionally standardize features before PCA (safer if dims vary in scale)
        if scale_pca:
            sc.pp.scale(adata, zero_center=True, max_value=None)

        # Compute PCA scores; cap n_pcs to rank
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
    n_neighbors: int = 15,
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

    # Neighbor graph
    if cluster_space == "pca":
        sc.pp.neighbors(adata, n_neighbors=n_neighbors, use_rep="X_pca", metric=metric)
    else:
        sc.pp.neighbors(adata, n_neighbors=n_neighbors, use_rep='X',metric=metric)

    # Leiden across resolutions
    out = {}
    for res in resolutions:
        key = f"leiden_{res:g}"
        sc.tl.leiden(adata, resolution=float(res), key_added=key, random_state=random_state)
        out[float(res)] = pd.Series(
            adata.obs[key].astype(str).values, index=adata.obs_names, name="cluster"
        )
    return out


    
# -----------------------------
# g:Profiler enrichment per set
# -----------------------------
def gprofiler_enrich_clusters_for_set(
    clusters_by_res: Dict[float, pd.Series],
    go_sources: Sequence[str],
    alpha: float = 0.05,
    min_cluster_size: int = 10,
    organism: str = "drerio",
    correction: str = "fdr",          # {"fdr","bonferroni","g_SCS"}
    domain_scope: str = "annotated",  # {"annotated","known"}
    exclude_iea: bool = False,        # True = exclude IEA annotations
    set_name: str = "set",
    gp: Optional[GProfiler] = None,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Run g:Profiler GO enrichment for each cluster at each resolution.

    Returns:
      details_df: rows = significant terms (p_value <= alpha under chosen correction)
      summary_df: counts per (set, resolution) — total and unique terms
    """
    if gp is None:
        gp = GProfiler(return_dataframe=True)

    records: List[dict] = []
    summary_rows: List[dict] = []

    for res, clusters in clusters_by_res.items():
        terms_seen_this_res: List[str] = []

        for cl in sorted(clusters.unique(), key=str):
            genes = clusters.index[clusters == cl].tolist()
            if len(genes) < min_cluster_size:
                continue

            try:
                enr = gp.profile(
                    organism=organism,
                    query=genes,
                    sources=list(go_sources),
                    user_threshold=alpha,
                    significance_threshold_method=correction,
                    domain_scope=domain_scope,
                    no_evidences=exclude_iea,
                    all_results=True,
                )
            except Exception as e:
                warnings.warn(f"g:Profiler failed for set={set_name}, res={res}, cluster={cl}: {e}")
                continue

            if enr is None or len(enr) == 0:
                continue

            df = enr.copy()
            go_id_col = "native" if "native" in df.columns else ("term_id" if "term_id" in df.columns else None)
            name_col = "name" if "name" in df.columns else ("description" if "description" in df.columns else None)
            p_col = "p_value" if "p_value" in df.columns else None

            if p_col is None or go_id_col is None or name_col is None:
                warnings.warn(f"g:Profiler result missing expected columns (have: {list(df.columns)}); skipping cluster.")
                continue

            if "source" in df.columns:
                df = df[df["source"].isin(go_sources)].copy()

            df = df[np.isfinite(df[p_col]) & (df[p_col] <= alpha)].copy()
            if df.empty:
                continue

            for _, row in df.iterrows():
                go_id = str(row[go_id_col]) if pd.notna(row[go_id_col]) else None
                term_key = go_id if go_id else f"{row.get('source','GO')}::{row[name_col]}"
                terms_seen_this_res.append(term_key)

                records.append({
                    "set": set_name,
                    "resolution": float(res),
                    "cluster": str(cl),
                    "library": str(row.get("source", "GO")),
                    "term": str(row[name_col]),
                    "go_id": go_id,
                    "pvalue_adj": float(row[p_col]),
                    "term_size": int(row.get("term_size")) if pd.notna(row.get("term_size")) else None,
                    "query_size": int(row.get("query_size")) if pd.notna(row.get("query_size")) else len(genes),
                    "intersection_size": int(row.get("intersection_size")) if pd.notna(row.get("intersection_size")) else None,
                })

        summary_rows.append({
            "set": set_name,
            "resolution": float(res),
            "n_sig_terms_total": len(terms_seen_this_res),
            "n_sig_terms_unique": len(set(terms_seen_this_res)),
        })

    details_df = pd.DataFrame.from_records(records).sort_values(
        ["set", "resolution", "cluster", "library", "pvalue_adj"], ignore_index=True
    )
    summary_df = pd.DataFrame.from_records(summary_rows).sort_values(
        ["set", "resolution"], ignore_index=True
    )
    return details_df, summary_df


# -----------------------------
# Plotting
# -----------------------------
def plot_summary(summary_df: pd.DataFrame, out_png: Optional[str] = None):
    """
    Line plot: for each (embedding set [+ clustering space] label),
    number of *unique* significant terms vs resolution.
    """
    if summary_df.empty:
        warnings.warn("Empty summary_df; nothing to plot.")
        return

    df = summary_df.copy()
    if "n_sig_terms_unique" not in df.columns:
        raise ValueError("summary_df must contain 'n_sig_terms_unique'.")

    labels = df["set"].unique().tolist()

    plt.figure(figsize=(8.0, 5.2))
    for label in labels:
        sub = df[df["set"] == label].sort_values("resolution")
        plt.plot(sub["resolution"].values, sub["n_sig_terms_unique"].values, marker="o", label=label)

    plt.xlabel("Leiden resolution")
    plt.ylabel("Unique significant GO terms (q ≤ threshold)")
    plt.title("GO enrichment vs. clustering resolution (g:Profiler; zebrafish)")
    plt.grid(True, linestyle="--", alpha=0.4)
    plt.legend(title="Series", frameon=False, ncol=1)
    plt.tight_layout()
    if out_png:
        plt.savefig(out_png, dpi=150)
    plt.close()

    # -----------------------------
# Orchestration across sets and spaces
# -----------------------------
def run_pipeline(
    embedding_sets: Mapping[str, pd.DataFrame],
    resolutions: Sequence[float],
    go_sources: Sequence[str] = DEFAULT_GO_SOURCES,
    alpha: float = 0.05,
    min_cluster_size: int = 10,
    n_neighbors: int = 15,
    metric: str = "cosine",  # (unused in revised version if you split metrics; keep if you didn't)
    random_state: int = 0,
    align_by_intersection: bool = True,
    out_prefix: str = "enrichment",
    organism: str = "drerio",
    correction: str = "fdr",
    domain_scope: str = "annotated",
    exclude_iea: bool = False,
    # --- clustering space controls you already added ---
    cluster_spaces: Sequence[str] = ("embedding",),
    metric_embedding: str = "cosine",
    metric_pca: str = "euclidean",
    n_pcs: int = 50,
    scale_pca: bool = False,
    pca_svd_solver: str = "arpack",
    # --- NEW: explained variance plot controls ---
    make_pca_ev_plot: bool = True,
    ev_max_pcs: int = 100,          # how far to extend the EV curve
    ev_scale_pca: Optional[bool] = None,  # None → inherit scale_pca
    ev_mode: str = "cumulative",    # "cumulative" or "per_pc"
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


    # -----------------------------
    # (NEW) PCA explained variance across sets
    # -----------------------------
    if make_pca_ev_plot:
        ev_scale = scale_pca if ev_scale_pca is None else bool(ev_scale_pca)
        ev_df = compute_pca_explained_variance_across_sets(
            embedding_sets=sets,
            max_pcs=ev_max_pcs,
            scale=ev_scale,
            pca_svd_solver=pca_svd_solver,
        )
        ev_csv = f"{out_prefix}_pca_explained_variance.csv"
        ev_png = f"{out_prefix}_pca_explained_variance.png"
        ev_df.to_csv(ev_csv, index=False)

        # Draw a vertical line at the n_pcs you use for PCA clustering (if applicable)
        vline = n_pcs if ("pca" in set(map(str, cluster_spaces))) else None
        plot_pca_explained_variance(
            ev_df,
            out_png=ev_png,
            mode=ev_mode,
            vline_at=vline,
            title="PCA explained variance across embedding sets",
        )
        print(f"[INFO] Wrote: {ev_csv}, {ev_png}")
        
        
    gp = GProfiler(return_dataframe=True)

    all_details = []
    all_summary = []

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

            print(f"[INFO] Enrichment (g:Profiler) for '{label}' ...")
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
                gp=gp,
            )
            all_details.append(details_df)
            all_summary.append(summary_df)

    details = pd.concat(all_details, axis=0, ignore_index=True) if len(all_details) else pd.DataFrame()
    summary = pd.concat(all_summary, axis=0, ignore_index=True) if len(all_summary) else pd.DataFrame()

    # Save & plot
    details.to_csv(f"{out_prefix}_details.csv", index=False)
    summary.to_csv(f"{out_prefix}_summary.csv", index=False)
    plot_summary(summary, out_png=f"{out_prefix}_summary_plot.png")

    return details, summary



# -----------------------------
# Example usage 
# -----------------------------
inputdir = "/path/to/Geneformer_RQfork/embeddings/extracted_embeddings"

RQ_last = os.path.join(inputdir, "RNAquarium_last_layer.csv")
RQ_second2last = os.path.join(inputdir, "RNAquarium_second_to_last_layer.csv")

RQ_scrambled_last = os.path.join(inputdir, "RNAquarium_scrambled_last_layer.csv")
RQ_scrambled_second2last = os.path.join(inputdir, "RNAquarium_scrambled_second_to_last_layer.csv")

RQ_randomized_last = os.path.join(inputdir, "RNAquarium_randomized_last_layer.csv")
RQ_randomized_second2last = os.path.join(inputdir, "RNAquarium_randomized_second_to_last_layer.csv")

RQ_zebrahubfinetuned_last = os.path.join(inputdir, "RNAquarium_zebrahub_finetune_last_layer.csv")
RQ_zebrahubfinetuned_second2last = os.path.join(inputdir, "RNAquarium_zebrahub_finetune_second_to_last_layer.csv")

RQ_zebrahubData_last = os.path.join(inputdir, "RNAquarium_zebrahubData_last_layer.csv")
RQ_zebrahubData_second2last = os.path.join(inputdir, "RNAquarium_zebrahubData_second_to_last_layer.csv")

import pandas as pd
# Example CSVs with index as ENSDARG IDs, columns=embedding dims
emb_a = pd.read_csv(RQ_last, index_col=0)
emb_b = pd.read_csv(RQ_randomized_last, index_col=0)
embedding_sets = {"RQ": emb_a, "RQ_randomized": emb_b}
resolutions = [0.2, 0.4, 0.6, 0.8, 1.0]

# # (A) Raw-embedding clustering only
# details, summary = run_pipeline(
#     embedding_sets=embedding_sets,
#     resolutions=resolutions,
#     cluster_spaces=("embedding",),
#     metric_embedding="cosine",
#     out_prefix="enrichment_raw",
# )

# # (B) PCA clustering only
# details, summary = run_pipeline(
#     embedding_sets=embedding_sets,
#     resolutions=resolutions,
#     cluster_spaces=("pca",),
#     n_pcs=80, scale_pca=True, metric_pca="euclidean",
#     out_prefix="enrichment_pca",
# )

# (C) Compare raw vs PCA side-by-side
details, summary = run_pipeline(
    embedding_sets=embedding_sets,
    resolutions=resolutions,
    cluster_spaces=("embedding","pca"),
    n_pcs=225, scale_pca=False,
    out_prefix="enrichment_both",
    make_pca_ev_plot=True,   # <— enable EV plot
    ev_max_pcs=255,          # length of EV curve (per set, capped by rank)
    ev_scale_pca=None,       # inherit scale_pca (set True/False to override)
    ev_mode="cumulative",    # or "per_pc" for a classic scree plot
)
print("Wrote: enrichment_details.csv, enrichment_summary.csv, enrichment_summary_plot.png")

