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
    # If multiple Ensembl IDs map to the same gene symbol, collapse by mean.
    if df.index.has_duplicates:
        df = df.groupby(level=0).mean(numeric_only=True)
    return df


def _collapse_duplicate_axes_square(df: pd.DataFrame) -> pd.DataFrame:
    # Collapse duplicates on both axes of a square matrix (mean).
    if df.index.has_duplicates:
        df = df.groupby(level=0).mean(numeric_only=True)
    if df.columns.has_duplicates:
        df = df.groupby(level=0, axis=1).mean(numeric_only=True)  # ok for pandas<=2.x
    return df


def build_multiplex_threshold_network_for_cytoscape(
    simA: pd.DataFrame,
    simB: pd.DataFrame,
    nodes: Iterable[str],
    *,
    nameA: str = "embeddings",
    nameB: str = "coexpr",
    threshA: float = 0.0,
    threshB: float = 0.0,
    directed: bool = False,
    keep_self_loops: bool = False,
    output_dir: Optional[Union[str, Path]] = None,
    write_graphml_path: Optional[Union[str, Path]] = "multiplex.graphml",
    write_csv_prefix: Optional[Union[str, Path]] = "multiplex",
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Build a 2-layer multiplex Cytoscape network from two similarity matrices using thresholds:
      - add an edge in layer A if simA(i,j) > threshA
      - add an edge in layer B if simB(i,j) > threshB

    Exports GraphML + nodes/edges CSVs. Returns (nodes_df, edges_df).
    """
    nodes = list(nodes)

    # ----- output paths -----
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
        base_str = base.name  # preserve dots, etc.
        nodes_csv_path = base.with_name(f"{base_str}_nodes.csv")
        edges_csv_path = base.with_name(f"{base_str}_edges.csv")

    # ----- align matrices to nodes -----
    A = simA.reindex(index=nodes, columns=nodes).astype(float)
    B = simB.reindex(index=nodes, columns=nodes).astype(float)

    if not keep_self_loops:
        np.fill_diagonal(A.values, np.nan)
        np.fill_diagonal(B.values, np.nan)

    # ----- edge extraction -----
    def _edges_from_matrix(M: pd.DataFrame, layer_name: str, thresh: float) -> pd.DataFrame:
        arr = M.values
        ids = M.index.to_numpy()

        if directed:
            i_idx, j_idx = np.where(np.isfinite(arr) & (arr > thresh))
        else:
            iu, ju = np.triu_indices_from(arr, k=(0 if keep_self_loops else 1))
            vals = arr[iu, ju]
            m = np.isfinite(vals) & (vals > thresh)
            i_idx, j_idx = iu[m], ju[m]
            vals = vals[m]

        return pd.DataFrame(
            {
                "source": ids[i_idx],
                "target": ids[j_idx],
                "type": layer_name,
                "weight": arr[i_idx, j_idx].astype(float) if directed else vals.astype(float),
            }
        )

    edgesA = _edges_from_matrix(A, nameA, threshA)
    edgesB = _edges_from_matrix(B, nameB, threshB)
    edges = pd.concat([edgesA, edgesB], ignore_index=True)

    # Safety: avoid accidental duplicates
    edges = edges.drop_duplicates(subset=["source", "target", "type"]).reset_index(drop=True)

    # ----- node table -----
    nodes_df = pd.DataFrame({"id": nodes})
    nodes_df["label"] = nodes_df["id"]

    if edges.empty:
        # no edges at all
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
            # Union degree: count unique neighbors ignoring layer
            any_edges = edges.groupby(["source", "target"], as_index=False)["weight"].sum()
            deg_counts = pd.concat([any_edges["source"], any_edges["target"]]).value_counts()
            nodes_df["degree_any"] = nodes_df["id"].map(deg_counts).fillna(0).astype(int)

            for nm in (nameA, nameB):
                sub = edges[edges["type"] == nm]
                strn = sub.groupby("source")["weight"].sum().add(
                    sub.groupby("target")["weight"].sum(), fill_value=0.0
                )
                nodes_df[f"strength_{nm}"] = nodes_df["id"].map(strn).fillna(0.0)

    # ----- write files -----
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
        description="Build 2-layer (multiplex) gene network for Cytoscape using thresholds on embeddings cosine-similarity and coexpression correlation."
    )
    parser.add_argument("--genes-csv", required=True, help='CSV containing a "gene" column.')
    parser.add_argument("--top-n", type=int, required=True, help='Use the top N rows from the "gene" column.')
    parser.add_argument("--emb-thresh", type=float, required=True, help="Keep embedding edges with cosine similarity > this threshold.")
    parser.add_argument("--coexpr-thresh", type=float, required=True, help="Keep coexpression edges with correlation > this threshold.")
    parser.add_argument("--output-dir", default="for_cytoscape", help="Output directory for GraphML and CSVs.")

    # Keep your existing defaults, but make them overridable:
    parser.add_argument(
        "--mapping-csv",
        default="/path/to/Geneformer_RQfork/data/RNAquarium/xref/drerio_z12-geneid_to_z11-ENSDARG.csv",
        help="CSV mapping Ensembl_gene_id -> gene_id.",
    )
    parser.add_argument(
        "--embeddings-csv",
        default="/path/to/Geneformer_RQfork/embeddings/extracted_embeddings/RNAquarium_second_to_last_layer.csv",
        help="Embeddings CSV (rows=genes, cols=embedding dims).",
    )
    parser.add_argument(
        "--coexpr-npz",
        default="/path/to/Geneformer_RQfork/embeddings/coexpression/corr_matrix/corr_matrix.npz",
        help='NPZ containing arrays "corr" and "genes".',
    )
    parser.add_argument(
        "--output-prefix",
        default="thresholded_graph",
        help="Output prefix for the graphml and csv files.",
    )

    args = parser.parse_args()

    # ----- load mapping -----
    mapping_df = pd.read_csv(args.mapping_csv)
    gene_col = "gene_id"
    ensembl_col = "Ensembl_gene_id"
    mapper = mapping_df.set_index(ensembl_col)[gene_col].to_dict()

    # ----- read top-N genes from input CSV -----
    genes_df = pd.read_csv(args.genes_csv)
    if "gene" not in genes_df.columns:
        raise ValueError(f'Input file {args.genes_csv} must contain a column named "gene".')

    raw_genes = genes_df["gene"].head(args.top_n).astype(str).tolist()

    # Map Ensembl IDs -> gene symbols if needed; otherwise keep as-is
    selected = []
    for g in raw_genes:
        g = g.strip()
        if not g:
            continue
        selected.append(mapper.get(g, g))

    selected = _dedupe_preserve_order(selected)
    if not selected:
        raise ValueError("No valid genes found in the top-N of the input CSV.")

    # ----- load embeddings -----
    t0 = time.time()
    embeddings = pd.read_csv(args.embeddings_csv, index_col=0)
    embeddings = embeddings.rename(index=mapper)
    embeddings = _collapse_duplicate_index(embeddings)
    t1 = time.time()
    print(f"[info] Loaded embeddings in {t1 - t0:.2f}s. shape={embeddings.shape}")

    # ----- load coexpression corr matrix -----
    t0 = time.time()
    npz = np.load(args.coexpr_npz, allow_pickle=True)
    coexpr_corr = npz["corr"]
    coexpr_genes = npz["genes"]
    coexpr_corr_df = pd.DataFrame(coexpr_corr, index=coexpr_genes, columns=coexpr_genes)
    coexpr_corr_df = coexpr_corr_df.rename(index=mapper).rename(columns=mapper)
    coexpr_corr_df = _collapse_duplicate_axes_square(coexpr_corr_df)
    t1 = time.time()
    print(f"[info] Loaded coexpr matrix in {t1 - t0:.2f}s. shape={coexpr_corr_df.shape}")

    # ----- keep only nodes that exist in at least one source -----
    ok_emb = set(embeddings.index)
    ok_coexpr = set(coexpr_corr_df.index).intersection(coexpr_corr_df.columns)
    nodes = [g for g in selected if (g in ok_emb) or (g in ok_coexpr)]
    dropped = [g for g in selected if g not in nodes]
    if dropped:
        print(f"[warning] Dropped {len(dropped)} genes not found in embeddings or coexpr matrix. Example: {dropped[:10]}")

    if not nodes:
        raise ValueError("None of the selected genes are present in embeddings or coexpression matrix.")

    # ----- embeddings cosine similarity on the available subset -----
    emb_nodes = [g for g in nodes if g in ok_emb]
    if emb_nodes:
        X = embeddings.loc[emb_nodes].select_dtypes(include="number").astype(float)
        S = cosine_similarity(X.values)
        emb_sim_sub = pd.DataFrame(S, index=emb_nodes, columns=emb_nodes)
    else:
        emb_sim_sub = pd.DataFrame(index=[], columns=[])

    # Expand to full node set (missing values = NaN)
    emb_sim_df = emb_sim_sub.reindex(index=nodes, columns=nodes)

    # ----- coexpr subset to node set -----
    coexpr_sub_df = coexpr_corr_df.reindex(index=nodes, columns=nodes)

    # ----- output naming -----
    def _fmt(x: float) -> str:
        # Compact + filename-friendly-ish
        return f"{x:g}".replace(".", "p")

    output_prefix = args.output_prefix
    suffix = f".n{len(nodes)}.embgt{_fmt(args.emb_thresh)}.coexprgt{_fmt(args.coexpr_thresh)}"
    graphml_name = f"{output_prefix}{suffix}.graphml"
    csv_prefix = f"{output_prefix}{suffix}"

    # ----- build + export multiplex -----
    build_multiplex_threshold_network_for_cytoscape(
        emb_sim_df,
        coexpr_sub_df,
        nodes,
        nameA="embeddings",
        nameB="coexpr",
        threshA=args.emb_thresh,
        threshB=args.coexpr_thresh,
        directed=False,
        keep_self_loops=False,
        output_dir=args.output_dir,
        write_graphml_path=graphml_name,
        write_csv_prefix=csv_prefix,
    )


if __name__ == "__main__":
    main()
