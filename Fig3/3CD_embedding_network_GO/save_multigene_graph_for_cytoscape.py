import time
import argparse
from pathlib import Path
from typing import Iterable, Optional, Tuple, Union

import numpy as np
import pandas as pd
import networkx as nx
from sklearn.metrics.pairwise import cosine_similarity


def _dedupe_preserve_order(items: Iterable[str]) -> list[str]:
    seen = set()
    out = []
    for x in items:
        if x not in seen:
            seen.add(x)
            out.append(x)
    return out


def _collapse_duplicate_index(df: pd.DataFrame) -> pd.DataFrame:
    if df.index.has_duplicates:
        df = df.groupby(level=0).mean(numeric_only=True)
    return df


def _collapse_duplicate_axes_square(df: pd.DataFrame) -> pd.DataFrame:
    if df.index.has_duplicates:
        df = df.groupby(level=0).mean(numeric_only=True)
    if df.columns.has_duplicates:
        df = df.groupby(level=0, axis=1).mean(numeric_only=True)
    return df


def _parse_genes_arg(tokens: list[str]) -> list[str]:
    out: list[str] = []
    for tok in tokens:
        if tok is None:
            continue
        for part in str(tok).split(","):
            part = part.strip()
            if part:
                out.append(part)
    return out


def _normalize_percentile(p: float, name: str) -> float:
    p = float(p)
    if 0.0 <= p <= 1.0:
        p = p * 100.0
    if not (0.0 <= p <= 100.0):
        raise ValueError(f"{name} must be in [0, 100] (or [0, 1] as a fraction). Got: {p}")
    return p


def _percentile_threshold_from_square(
    M: pd.DataFrame,
    percentile: float,
    *,
    directed: bool,
    keep_self_loops: bool,
) -> float:
    arr = M.values.astype(float, copy=False)

    if directed:
        mask = np.isfinite(arr)
        if not keep_self_loops:
            mask &= ~np.eye(arr.shape[0], dtype=bool)
        vals = arr[mask]
    else:
        k = 0 if keep_self_loops else 1
        iu, ju = np.triu_indices_from(arr, k=k)
        vals = arr[iu, ju]
        vals = vals[np.isfinite(vals)]

    if vals.size == 0:
        return float("inf")

    return float(np.nanpercentile(vals, percentile))


def _force_connect_each_node_to_best_query(
    edges: pd.DataFrame,
    M: pd.DataFrame,
    nodes: list[str],
    query_genes: list[str],
    layer_name: str,
    *,
    directed: bool,
    keep_self_loops: bool,
) -> Tuple[pd.DataFrame, int, int]:
    """
    Ensure every node has an edge (in this layer) to the query gene with the highest weight.
    If the edge is missing due to thresholding, add it ("override threshold").

    Returns (edges_updated, n_forced_added, n_failed).
    """
    if edges is None:
        edges = pd.DataFrame(columns=["source", "target", "type", "weight"])

    pos = {g: i for i, g in enumerate(nodes)}
    q_in_nodes = [q for q in query_genes if q in pos]
    if not q_in_nodes:
        return edges, 0, 0

    # Build a fast lookup for existing edges in this layer.
    # For undirected graphs we store edges in (source,target) orientation consistent with node order.
    existing = set(
        zip(
            edges["source"].astype(str),
            edges["target"].astype(str),
            edges["type"].astype(str),
        )
    )

    forced_rows = []
    n_failed = 0

    for u in nodes:
        # Candidates are query genes; exclude self to avoid self-loop
        candidates = q_in_nodes
        if u in candidates:
            candidates = [x for x in candidates if x != u]

        if not candidates:
            # e.g. only one query gene and u is that gene
            continue

        # Weights from u -> candidates (may include NaNs if missing in this layer)
        try:
            s = M.loc[u, candidates]
        except KeyError:
            n_failed += 1
            continue

        if isinstance(s, pd.DataFrame):
            # shouldn't happen, but guard anyway
            s = s.iloc[0]

        s = s.astype(float)
        s = s[np.isfinite(s.to_numpy())]
        if s.empty:
            n_failed += 1
            continue

        v = str(s.idxmax())
        w = float(s.loc[v])

        if not keep_self_loops and u == v:
            continue

        if directed:
            src, tgt = u, v
        else:
            iu, iv = pos[u], pos[v]
            if iu < iv:
                src, tgt = u, v
            else:
                src, tgt = v, u

        key = (str(src), str(tgt), str(layer_name))
        if key in existing:
            continue

        forced_rows.append({"source": src, "target": tgt, "type": layer_name, "weight": w})
        existing.add(key)

    if forced_rows:
        edges = pd.concat([edges, pd.DataFrame(forced_rows)], ignore_index=True)

    return edges, len(forced_rows), n_failed


def build_multiplex_percentile_threshold_network_for_cytoscape(
    simA: pd.DataFrame,
    simB: pd.DataFrame,
    nodes: Iterable[str],
    *,
    query_genes: list[str],
    nameA: str = "embeddings",
    nameB: str = "coexpr",
    pctlA: float = 0.0,
    pctlB: float = 0.0,
    includeA: bool = True,
    includeB: bool = True,
    directed: bool = False,
    keep_self_loops: bool = False,
    node_table: Optional[pd.DataFrame] = None,  # must include "id"
    output_dir: Optional[Union[str, Path]] = None,
    write_graphml_path: Optional[Union[str, Path]] = "multiplex.graphml",
    write_csv_prefix: Optional[Union[str, Path]] = "multiplex",
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Percentile-thresholded multiplex export, plus:
      - FORCED edges: each node connects to its max-weight query gene (per included layer),
        overriding thresholds as needed.
    """
    nodes = list(nodes)

    outdir = Path(output_dir).resolve() if output_dir else Path.cwd()
    if output_dir:
        outdir.mkdir(parents=True, exist_ok=True)

    graphml_path = None
    if write_graphml_path:
        wp = Path(write_graphml_path)
        graphml_path = wp if wp.is_absolute() else (outdir / wp.name)

    nodes_csv_path = edges_csv_path = None
    if write_csv_prefix:
        base = Path(write_csv_prefix)
        if not base.is_absolute():
            base = outdir / base
        base_str = base.name
        nodes_csv_path = base.with_name(f"{base_str}_nodes.csv")
        edges_csv_path = base.with_name(f"{base_str}_edges.csv")

    # Align matrices
    A = simA.reindex(index=nodes, columns=nodes).astype(float)
    B = simB.reindex(index=nodes, columns=nodes).astype(float)

    if not keep_self_loops:
        np.fill_diagonal(A.values, np.nan)
        np.fill_diagonal(B.values, np.nan)

    # Percentile -> absolute threshold
    pA = _normalize_percentile(pctlA, "emb-thresh")
    pB = _normalize_percentile(pctlB, "coexpr-thresh")

    thrA = None
    thrB = None
    if includeA:
        thrA = _percentile_threshold_from_square(A, pA, directed=directed, keep_self_loops=keep_self_loops)
        print(f"[info] {nameA} threshold: percentile={pA:g} => value={thrA:.6g} (keep weight >= value)")
    else:
        print(f"[info] {nameA} layer excluded (no {nameA} edges will be written).")

    if includeB:
        thrB = _percentile_threshold_from_square(B, pB, directed=directed, keep_self_loops=keep_self_loops)
        print(f"[info] {nameB} threshold: percentile={pB:g} => value={thrB:.6g} (keep weight >= value)")
    else:
        print(f"[info] {nameB} layer excluded (no {nameB} edges will be written).")

    def _edges_from_matrix(M: pd.DataFrame, layer_name: str, thresh: float) -> pd.DataFrame:
        arr = M.values
        ids = M.index.to_numpy()

        if directed:
            m = np.isfinite(arr) & (arr >= thresh)
            if not keep_self_loops:
                m &= ~np.eye(arr.shape[0], dtype=bool)
            i_idx, j_idx = np.where(m)
            weights = arr[i_idx, j_idx].astype(float)
        else:
            iu, ju = np.triu_indices_from(arr, k=(0 if keep_self_loops else 1))
            vals = arr[iu, ju]
            m = np.isfinite(vals) & (vals >= thresh)
            i_idx, j_idx = iu[m], ju[m]
            weights = vals[m].astype(float)

        return pd.DataFrame({"source": ids[i_idx], "target": ids[j_idx], "type": layer_name, "weight": weights})

    edges_parts = []
    if includeA:
        edges_parts.append(_edges_from_matrix(A, nameA, thrA))
    if includeB:
        edges_parts.append(_edges_from_matrix(B, nameB, thrB))

    edges = (
        pd.concat(edges_parts, ignore_index=True)
        if edges_parts
        else pd.DataFrame(columns=["source", "target", "type", "weight"])
    )
    edges = edges.drop_duplicates(subset=["source", "target", "type"]).reset_index(drop=True)

    # ---- FORCED CONNECTIONS (override thresholds) ----
    # Ensure each node has an edge to its best query gene for each included layer.
    if includeA:
        edges, n_forced, n_failed = _force_connect_each_node_to_best_query(
            edges, A, nodes, query_genes, nameA, directed=directed, keep_self_loops=keep_self_loops
        )
        print(f"[info] Forced {n_forced} {nameA} edges to best query gene (failed for {n_failed} nodes).")
    if includeB:
        edges, n_forced, n_failed = _force_connect_each_node_to_best_query(
            edges, B, nodes, query_genes, nameB, directed=directed, keep_self_loops=keep_self_loops
        )
        print(f"[info] Forced {n_forced} {nameB} edges to best query gene (failed for {n_failed} nodes).")

    edges = edges.drop_duplicates(subset=["source", "target", "type"]).reset_index(drop=True)

    # Node table
    nodes_df = pd.DataFrame({"id": nodes})
    nodes_df["label"] = nodes_df["id"]

    # Edge-derived stats
    if edges.empty:
        if directed:
            nodes_df["in_degree_any"] = 0
            nodes_df["out_degree_any"] = 0
            nodes_df[f"in_strength_{nameA}"] = 0.0
            nodes_df[f"out_strength_{nameA}"] = 0.0
            nodes_df[f"in_strength_{nameB}"] = 0.0
            nodes_df[f"out_strength_{nameB}"] = 0.0
        else:
            nodes_df["degree_any"] = 0
            nodes_df[f"strength_{nameA}"] = 0.0
            nodes_df[f"strength_{nameB}"] = 0.0
    else:
        if directed:
            any_edges = edges.groupby(["source", "target"], as_index=False)["weight"].sum()
            indeg = any_edges.groupby("target").size().rename("in_degree_any")
            outdeg = any_edges.groupby("source").size().rename("out_degree_any")

            nodes_df["in_degree_any"] = nodes_df["id"].map(indeg).fillna(0).astype(int)
            nodes_df["out_degree_any"] = nodes_df["id"].map(outdeg).fillna(0).astype(int)

            for nm in (nameA, nameB):
                sub = edges[edges["type"] == nm]
                in_str = sub.groupby("target")["weight"].sum().rename(f"in_strength_{nm}")
                out_str = sub.groupby("source")["weight"].sum().rename(f"out_strength_{nm}")
                nodes_df[f"in_strength_{nm}"] = nodes_df["id"].map(in_str).fillna(0.0)
                nodes_df[f"out_strength_{nm}"] = nodes_df["id"].map(out_str).fillna(0.0)
        else:
            any_edges = edges.groupby(["source", "target"], as_index=False)["weight"].sum()
            deg_counts = pd.concat([any_edges["source"], any_edges["target"]]).value_counts()
            nodes_df["degree_any"] = nodes_df["id"].map(deg_counts).fillna(0).astype(int)

            for nm in (nameA, nameB):
                sub = edges[edges["type"] == nm]
                strn = sub.groupby("source")["weight"].sum().add(
                    sub.groupby("target")["weight"].sum(), fill_value=0.0
                )
                nodes_df[f"strength_{nm}"] = nodes_df["id"].map(strn).fillna(0.0)

    # Merge in extra node metadata
    if node_table is not None and not node_table.empty:
        if "id" not in node_table.columns:
            raise ValueError('node_table must include a column named "id".')
        extra = node_table.copy()
        dup_cols = (set(extra.columns) & set(nodes_df.columns)) - {"id"}
        if dup_cols:
            extra = extra.rename(columns={c: f"meta_{c}" for c in dup_cols})
        nodes_df = nodes_df.merge(extra, on="id", how="left")

    # Write files
    if graphml_path:
        G = nx.MultiDiGraph() if directed else nx.MultiGraph()

        for _, r in nodes_df.iterrows():
            attrs = r.drop(labels=["id"]).to_dict()
            G.add_node(r["id"], **attrs)

        for _, r in edges.iterrows():
            G.add_edge(
                r["source"],
                r["target"],
                key=r["type"],
                type=str(r["type"]),
                weight=float(r["weight"]),
            )

        nx.write_graphml(G, graphml_path)
        print(f"[info] GraphML written to: {graphml_path}")

    if nodes_csv_path and edges_csv_path:
        nodes_df.to_csv(nodes_csv_path, index=False)
        edges.to_csv(edges_csv_path, index=False)
        print(f"[info] Wrote CSVs:\n  {nodes_csv_path}\n  {edges_csv_path}")

    return nodes_df, edges


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Query gene list -> summed embedding score ranking -> take top-N (excluding query genes) -> "
            "export Cytoscape network on (query + top-N). Thresholding is percentile-based per layer. "
            "Additionally, each node is forced to connect to its max-weight query gene (per included layer)."
        )
    )

    parser.add_argument("--genes", required=True, nargs="+")
    parser.add_argument("--top-n", type=int, required=True)

    parser.add_argument(
        "--embedding-score",
        choices=["cosine_similarity", "cosine_distance"],
        default="cosine_similarity",
        help="Edge weights used in embeddings layer and for ranking.",
    )

    parser.add_argument(
        "--layers",
        choices=["both", "embeddings", "coexpr"],
        default="both",
        help='Which edge layers to export: "both" (default), "embeddings" only, or "coexpr" only.',
    )

    parser.add_argument(
        "--emb-thresh",
        type=float,
        default=0.0,
        help="Embedding edge percentile threshold (0-100). Keeps edges with weight >= this percentile. 0 keeps all.",
    )
    parser.add_argument(
        "--coexpr-thresh",
        type=float,
        default=0.0,
        help="Coexpression edge percentile threshold (0-100). Keeps edges with weight >= this percentile. 0 keeps all.",
    )

    parser.add_argument("--rank-csv", default="ranked_by_embedding_sum.csv")
    parser.add_argument("--output-dir", default="for_cytoscape")
    parser.add_argument("--output-prefix", default="multiplex_graph")

    parser.add_argument(
        "--mapping-csv",
        default="/path/to/Geneformer_RQfork/data/RNAquarium/xref/drerio_z12-geneid_to_z11-ENSDARG.csv",
    )
    parser.add_argument(
        "--embeddings-csv",
        default="/path/to/Geneformer_RQfork/embeddings/extracted_embeddings/RNAquarium_second_to_last_layer.csv",
    )
    parser.add_argument(
        "--coexpr-npz",
        default="/path/to/Geneformer_RQfork/embeddings/coexpression/corr_matrix/corr_matrix.npz",
    )

    args = parser.parse_args()

    include_emb = args.layers in ("both", "embeddings")
    include_coexpr = args.layers in ("both", "coexpr")

    outdir = Path(args.output_dir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    # ----- mapping -----
    mapping_df = pd.read_csv(args.mapping_csv)
    mapper = mapping_df.set_index("Ensembl_gene_id")["gene_id"].to_dict()

    # ----- parse query genes -----
    raw_genes = _parse_genes_arg(args.genes)
    query_genes = _dedupe_preserve_order([mapper.get(g.strip(), g.strip()) for g in raw_genes if g.strip()])

    # ----- embeddings (needed for ranking regardless) -----
    t0 = time.time()
    embeddings = pd.read_csv(args.embeddings_csv, index_col=0)
    embeddings = embeddings.rename(index=mapper)
    embeddings = _collapse_duplicate_index(embeddings)
    embeddings_num = embeddings.select_dtypes(include="number").astype(float)
    print(f"[info] Loaded embeddings in {time.time() - t0:.2f}s. shape={embeddings_num.shape}")

    # ----- summed embedding score vs query list -----
    query_in_emb = [g for g in query_genes if g in embeddings_num.index]
    missing_q = [g for g in query_genes if g not in embeddings_num.index]
    if missing_q:
        print(f"[warning] {len(missing_q)} query genes not found in embeddings (ignored for scoring). Example: {missing_q[:10]}")
    if not query_in_emb:
        raise ValueError("None of the query genes were found in embeddings, so ranking cannot be computed.")

    sim_all_q = cosine_similarity(embeddings_num.values, embeddings_num.loc[query_in_emb].values)
    score_mat = (1.0 - sim_all_q) if args.embedding_score == "cosine_distance" else sim_all_q
    sum_score_s = pd.Series(score_mat.sum(axis=1), index=embeddings_num.index, name="embedding_sum_score")
    ranked_all = sum_score_s.sort_values(ascending=False)

    # ----- optionally load coexpr (only if exporting it) -----
    coexpr_corr_df = None
    ok_coexpr = set()
    if include_coexpr:
        t0 = time.time()
        npz = np.load(args.coexpr_npz, allow_pickle=True)
        coexpr_corr = npz["corr"]
        coexpr_genes = npz["genes"]
        try:
            coexpr_genes = coexpr_genes.astype(str)
        except Exception:
            coexpr_genes = np.array([str(x) for x in coexpr_genes])

        coexpr_corr_df = pd.DataFrame(coexpr_corr, index=coexpr_genes, columns=coexpr_genes)
        coexpr_corr_df = coexpr_corr_df.rename(index=mapper).rename(columns=mapper)
        coexpr_corr_df = _collapse_duplicate_axes_square(coexpr_corr_df)
        ok_coexpr = set(coexpr_corr_df.index).intersection(coexpr_corr_df.columns)
        print(f"[info] Loaded coexpr matrix in {time.time() - t0:.2f}s. shape={coexpr_corr_df.shape}")

    ok_emb = set(embeddings_num.index)
    ok_for_output = set()
    if include_emb:
        ok_for_output |= ok_emb
    if include_coexpr:
        ok_for_output |= ok_coexpr

    # ----- choose top genes that are actually usable in the output layers -----
    ranked_nonquery = ranked_all.drop(labels=query_genes, errors="ignore")
    ranked_nonquery = ranked_nonquery[ranked_nonquery.index.isin(ok_for_output)]
    top_genes = ranked_nonquery.head(args.top_n).index.tolist()

    # Write full ranking table
    rank_df = ranked_all.reset_index()
    rank_df.columns = ["gene", "embedding_sum_score"]
    rank_df["is_query"] = rank_df["gene"].isin(query_genes).astype(int)
    (outdir / Path(args.rank_csv).name).write_text(rank_df.to_csv(index=False))
    print(f"[info] Wrote ranking CSV: {outdir / Path(args.rank_csv).name}")

    # ----- nodes = query + top (filtered to output-available genes) -----
    nodes = _dedupe_preserve_order(query_genes + top_genes)
    nodes_kept = [g for g in nodes if g in ok_for_output]
    dropped = [g for g in nodes if g not in nodes_kept]
    if dropped:
        print(f"[warning] Dropped {len(dropped)} nodes not present in requested output layer(s). Example: {dropped[:10]}")
    nodes = nodes_kept

    if not nodes:
        raise ValueError("No nodes remain after filtering to those present in the requested output layer(s).")

    # ----- build layer matrices on node set -----
    if include_emb:
        emb_nodes_present = [g for g in nodes if g in ok_emb]
        X = embeddings_num.loc[emb_nodes_present].values
        S = cosine_similarity(X)
        if args.embedding_score == "cosine_distance":
            S = 1.0 - S
        emb_mat_sub = pd.DataFrame(S, index=emb_nodes_present, columns=emb_nodes_present)
        emb_mat_df = emb_mat_sub.reindex(index=nodes, columns=nodes)
    else:
        emb_mat_df = pd.DataFrame(np.nan, index=nodes, columns=nodes)

    if include_coexpr and coexpr_corr_df is not None:
        coexpr_sub_df = coexpr_corr_df.reindex(index=nodes, columns=nodes)
    else:
        coexpr_sub_df = pd.DataFrame(np.nan, index=nodes, columns=nodes)

    # ----- node metadata -----
    rank_all = ranked_all.rank(method="min", ascending=False).astype(int)
    node_meta = pd.DataFrame({"id": nodes})
    node_meta["is_query"] = node_meta["id"].isin(query_genes).astype(int)
    node_meta["is_top"] = node_meta["id"].isin(top_genes).astype(int)
    node_meta["embedding_sum_score"] = node_meta["id"].map(sum_score_s)
    node_meta["embedding_sum_rank"] = node_meta["id"].map(rank_all)

    # ----- output naming -----
    def _fmt_token(x: Union[str, float, int]) -> str:
        return str(x).replace(" ", "_").replace("/", "_").replace(".", "p")

    suffix = (
        f".layers{args.layers}"
        f".q{len(query_genes)}.top{len(top_genes)}.nodes{len(nodes)}"
        f".{_fmt_token(args.embedding_score)}"
        f".embp{_fmt_token(args.emb_thresh)}"
        f".coexprp{_fmt_token(args.coexpr_thresh)}"
        f".forcedToQuery"
    )
    graphml_name = f"{args.output_prefix}{suffix}.graphml"
    csv_prefix = f"{args.output_prefix}{suffix}"

    build_multiplex_percentile_threshold_network_for_cytoscape(
        emb_mat_df,
        coexpr_sub_df,
        nodes,
        query_genes=query_genes,
        nameA="embeddings",
        nameB="coexpr",
        pctlA=args.emb_thresh,
        pctlB=args.coexpr_thresh,
        includeA=include_emb,
        includeB=include_coexpr,
        directed=False,
        keep_self_loops=False,
        node_table=node_meta,
        output_dir=args.output_dir,
        write_graphml_path=graphml_name,
        write_csv_prefix=csv_prefix,
    )


if __name__ == "__main__":
    main()
