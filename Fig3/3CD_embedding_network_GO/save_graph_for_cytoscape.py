import time
import pandas as pd
import numpy as np

import argparse

parser = argparse.ArgumentParser(description="Build multiplex network for Cytoscape.")
parser.add_argument(
    "--top-k",
    dest="top_k",
    type=int,
    default=None,
    help="Per-node top-k neighbors to keep in each layer (symmetrized). Omit for no top-k filtering."
)
args = parser.parse_args()
TOP_K = args.top_k

# define paths
mapping_file_path = "/path/to/Geneformer_RQfork/data/RNAquarium/xref/drerio_z12-geneid_to_z11-ENSDARG.csv"
embeddings_path = "/path/to/Geneformer_RQfork/embeddings/extracted_embeddings/RNAquarium_second_to_last_layer.csv"
coexpr_corr_path = "/path/to/Geneformer_RQfork/embeddings/coexpression/corr_matrix/corr_matrix.npz"

# load mapping table
mapping_df = pd.read_csv(mapping_file_path)
print(mapping_df.head())


# create mapping dict
gene_col = "gene_id"
ensembl_col = "Ensembl_gene_id"
mapper = mapping_df.set_index(ensembl_col)[gene_col].to_dict()


# load embeddings
time_start = time.time()
embeddings = pd.read_csv(embeddings_path, index_col=0)
# translate ensembl to gene name
embeddings=embeddings.rename(index=mapper)
time_end = time.time()
print(f"Time taken to load embeddings: {time_end - time_start} seconds")
print(f"embeddings shape: {embeddings.shape}")
embeddings.head()


# compute cosine distance between embeddings
from sklearn.metrics.pairwise import cosine_similarity
# df: rows = observations, columns = features
X = embeddings.select_dtypes(include="number").astype(float)
S = cosine_similarity(X)                   # numpy array (n_obs x n_obs)
emb_sim_df = pd.DataFrame(S, index=X.index, columns=X.index)
# translate ensembl to gene name 
print(f"emb_sim_df shape: {emb_sim_df.shape}")
emb_sim_df.head()


# load co-expression correlation matrix
time_start = time.time()
coexpr_corr = np.load(coexpr_corr_path, allow_pickle=True)["corr"]
coexpr_genes = np.load(coexpr_corr_path, allow_pickle=True)["genes"]
time_end = time.time()
# construct a pandas dataframe 
coexpr_corr_df = pd.DataFrame(coexpr_corr, index=coexpr_genes, columns=coexpr_genes)
# translate ensembl to gene name
coexpr_corr_df=coexpr_corr_df.rename(index=mapper)
coexpr_corr_df=coexpr_corr_df.rename(columns=mapper)
print(f"Time taken to load co-expression correlation matrix: {time_end - time_start} seconds")
print(f"corr_df shape: {coexpr_corr_df.shape}")
coexpr_corr_df.head()


# select genes
IFNs = {
"ifnphi1":"ENSDARG00000025607",
"ifnphi2":"ENSDARG00000069012",
"ifnphi3":"ENSDARG00000070676",
"ifnphi4":"ENSDARG00000100678",
"ifng1":"ENSDARG00000024211",
"ifng1r":"ENSDARG00000045671",  
}

TLRs = {
"tlr1":"ENSDARG00000100649",
"tlr2":"ENSDARG00000037758",
"tlr3":"ENSDARG00000016065",
"tlr4al":"ENSDARG00000075671",
"tlr4ba":"ENSDARG00000019742",
"tlr4bb":"ENSDARG00000022048",
"tlr5a":"ENSDARG00000044415",
"tlr5b":"ENSDARG00000052322",
#"tlr7":"ENSDARG00000068812",
#"tlr8a":"ENSDARG00000090119",
#"CU914164.1":"ENSDARG00000104832",
"tlr8b":"ENSDARG00000073675",
"tlr9":"ENSDARG00000044490",
"tlr18":"ENSDARG00000040249",
"tlr19.1":"ENSDARG00000026663",
#"tlr20.1":"ENSDARG00000115923",
"tlr20.2":"ENSDARG00000088701",
#"tlr20.3":"ENSDARG00000114057",
"tlr20.4":"ENSDARG00000092668",
"tlr21":"ENSDARG00000058045",
"tlr22":"ENSDARG00000104045",
}

proinflammatory = ["caspa", "cxcl8a", "ifng1", "ifnphi1", "ifnphi2", "il1b", "irf3", "irf7", "lta", "sting1", "pycard", "rarres3", "tnfa", "tnfb"]
antiviral = [#"defb2", 
            "foxo3b", "ifit8", "ifit14", "isg15", "mavs", "ifih1", "mxa", "mxb", "mxc", "nod2", "pkz", "prmt3", "rela", 
            #"ddx58", 
            "ripk2", "tbk1", "tlr3", 
            #"tlr7", 
            "tlr8a", "tlr22", "rsad2"]


# check if all gene names are in the mapping df
all_names  = mapping_df["gene_id"].to_list()
all_ENSids = mapping_df["Ensembl_gene_id"].to_list()
def check_in_list(list1, list2):
    '''check if all items in list1 is in list2'''
    for i in list1:
        if not i in list2:
            print(f"{i} not found")

check_in_list(IFNs.keys(), all_names)
check_in_list(IFNs.values(), all_ENSids)
check_in_list(TLRs.keys(), all_names)
check_in_list(TLRs.values(), all_ENSids)

check_in_list(proinflammatory, all_names)
check_in_list(antiviral, all_names)


# construct selected gene list
from collections import defaultdict
selected = []
node_types = defaultdict(list)
for gene in proinflammatory:
    selected.append(gene)
    node_types[gene].append("proinflammatory")
for gene in antiviral:
    selected.append(gene)
    node_types[gene].append("antiviral")
for gene in IFNs.keys():
    selected.append(gene)
    node_types[gene].append("IFN")
for gene in TLRs.keys():
    selected.append(gene)
    node_types[gene].append("TLR")
    


# export graph to cytoscape
import pandas as pd
import numpy as np
import networkx as nx
from pathlib import Path
from typing import Iterable, Optional, Tuple, Union, Mapping, Any
import re

def union_multiplex_edges_for_cytoscape(
    simA: pd.DataFrame,
    simB: pd.DataFrame,
    selected_obs: Iterable,
    *,
    nameA: str = "embA",
    nameB: str = "embB",
    directed: bool = False,
    keep_self_loops: bool = False,
    drop_pairs_where_both_zero: bool = True,
    # Outputs
    output_dir: Optional[Union[str, Path]] = None,
    write_graphml_path: Optional[Union[str, Path]] = "multiplex_union.graphml",
    write_csv_prefix: Optional[Union[str, Path]] = None,
    # NEW: node type support
    node_types: Optional[Mapping[Any, Union[str, Iterable[str]]]] = None,
    type_separator: str = "|",
    make_type_indicator_cols: bool = False,
    # NEW: per-node top-k threshold (applied to each layer independently). None = no change.
    top_k_per_node: Optional[int] = None,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Build a multiplex (two edge types) Cytoscape network from two similarity matrices
    by taking the union of edges and filling missing similarities with 0.
    Adds optional node 'types' (multi-value) attribute.

    If `top_k_per_node` is provided, for each node i in each layer we keep only the
    top-k highest-weight neighbors in row i (self-loops are always excluded when
    keep_self_loops=False). For undirected graphs, the mask is symmetrized so an
    edge (i, j) is kept if it is in the top-k of i OR j.

    Returns:
        nodes_df, edges_df
          nodes_df includes: id, label, degree/strength (union-based), and optional
          type columns.
        edges_df columns: source, target, type, weight
    """
    # ---------- Resolve output paths ----------
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
        base_str = base.name  # preserve dots like ".topk5"
        nodes_csv_path = base.with_name(f"{base_str}_nodes.csv")
        edges_csv_path = base.with_name(f"{base_str}_edges.csv")

    # ---------- Align inputs ----------
    selected = list(selected_obs)
    okA = set(simA.index).intersection(simA.columns)
    okB = set(simB.index).intersection(simB.columns)
    ok_any = okA.union(okB)
    present = [o for o in selected if o in ok_any]
    dropped = sorted(set(selected) - set(present))
    if dropped:
        print(f"[warning] {len(dropped)} selected obs not found in A or B -> dropped: "
              f"{dropped[:8]}{'...' if len(dropped) > 8 else ''}")
    if not present:
        raise ValueError("None of the selected observations were found (as row+column) in simA or simB.")

    A = simA.reindex(index=present, columns=present).astype(float).fillna(0.0)
    B = simB.reindex(index=present, columns=present).astype(float).fillna(0.0)

    if not keep_self_loops:
        np.fill_diagonal(A.values, np.nan)
        np.fill_diagonal(B.values, np.nan)

    # ---------- Optional per-node top-k thresholding ----------
    def _apply_topk(M: pd.DataFrame, k: int) -> pd.DataFrame:
        arr = M.values.copy()
        n = arr.shape[0]
        # treat NaN as -inf so they won't be selected; self-loops already NaN
        arr_nonan = np.where(np.isnan(arr), -np.inf, arr)
        mask = np.zeros_like(arr_nonan, dtype=bool)

        for i in range(n):
            row = arr_nonan[i]
            finite_mask = np.isfinite(row)
            m = int(finite_mask.sum())
            if m == 0:
                continue
            kk = min(k, m)
            # indices of the kk largest values in the row
            top_idx = np.argpartition(row, -kk)[-kk:]
            # ignore positions that were -inf (non-finite)
            top_idx = top_idx[np.isfinite(row[top_idx])]
            mask[i, top_idx] = True

        # for undirected graphs, keep edges that are in top-k of i OR j
        if not directed:
            mask = mask | mask.T

        # zero out everything not selected; keep NaN on diagonal if needed
        out = np.where(mask, np.where(np.isfinite(arr), arr, 0.0), 0.0)
        if not keep_self_loops:
            np.fill_diagonal(out, np.nan)
        return pd.DataFrame(out, index=M.index, columns=M.columns)

    if top_k_per_node is not None:
        if top_k_per_node <= 0:
            raise ValueError("top_k_per_node must be a positive integer or None.")
        A = _apply_topk(A, top_k_per_node)
        B = _apply_topk(B, top_k_per_node)

    # ---------- Union of pairs ----------
    if directed:
        mask_union = ((~np.isnan(A.values)) & (A.values != 0)) | ((~np.isnan(B.values)) & (B.values != 0))
        i_idx, j_idx = np.where(mask_union)
    else:
        k = 0 if keep_self_loops else 1
        iu, ju = np.triu_indices_from(A.values, k=k)
        AU, BU = A.values[iu, ju], B.values[iu, ju]
        m = ((~np.isnan(AU)) & (AU != 0)) | ((~np.isnan(BU)) & (BU != 0))
        i_idx, j_idx = iu[m], ju[m]

    ids = np.array(present)
    edges_list = []
    for i, j in zip(i_idx, j_idx):
        wA = 0.0 if np.isnan(A.values[i, j]) else float(A.values[i, j])
        wB = 0.0 if np.isnan(B.values[i, j]) else float(B.values[i, j])
        s, t = ids[i], ids[j]
        if drop_pairs_where_both_zero and (wA == 0.0 and wB == 0.0):
            continue
        edges_list.append((s, t, nameA, wA))
        edges_list.append((s, t, nameB, wB))

    edges = pd.DataFrame(edges_list, columns=["source", "target", "type", "weight"])

    # ---------- Nodes ----------
    if directed:
        any_edges = edges.groupby(["source", "target"], as_index=False)["weight"].sum()
        indeg = any_edges.groupby("target").size().rename("in_degree_any")
        outdeg = any_edges.groupby("source").size().rename("out_degree_any")
        nodes = pd.DataFrame({"id": present})
        nodes["in_degree_any"] = nodes["id"].map(indeg).fillna(0).astype(int)
        nodes["out_degree_any"] = nodes["id"].map(outdeg).fillna(0).astype(int)
        for nm in (nameA, nameB):
            sub = edges[edges["type"] == nm]
            in_str = sub.groupby("target")["weight"].sum().rename(f"in_strength_{nm}")
            out_str = sub.groupby("source")["weight"].sum().rename(f"out_strength_{nm}")
            nodes[f"in_strength_{nm}"] = nodes["id"].map(in_str).fillna(0.0)
            nodes[f"out_strength_{nm}"] = nodes["id"].map(out_str).fillna(0.0)
    else:
        any_edges = edges.groupby(["source", "target"], as_index=False)["weight"].sum()
        deg_counts = pd.concat([any_edges["source"], any_edges["target"]]).value_counts()
        nodes = pd.DataFrame({"id": present})
        nodes["degree_any"] = nodes["id"].map(deg_counts).fillna(0).astype(int)
        for nm in (nameA, nameB):
            sub = edges[edges["type"] == nm]
            strn = sub.groupby("source")["weight"].sum().add(
                sub.groupby("target")["weight"].sum(), fill_value=0.0
            )
            nodes[f"strength_{nm}"] = nodes["id"].map(strn).fillna(0.0)

    nodes["label"] = nodes["id"]

    # ---------- Add node multi-types ----------
    if node_types is not None:
        def _norm_types(v) -> list:
            if v is None:
                return []
            if isinstance(v, str):
                s = v.strip()
                return [s] if s else []
            try:
                it = list(v)
            except TypeError:
                return []
            out = []
            for x in it:
                if x is None:
                    continue
                s = str(x).strip()
                if s:
                    out.append(s)
            # unique & stable order
            seen, unique = set(), []
            for s in out:
                if s not in seen:
                    seen.add(s); unique.append(s)
            return unique

        tmap = {k: _norm_types(v) for k, v in node_types.items()}

        nodes["types"] = nodes["id"].map(lambda i: type_separator.join(tmap.get(i, [])))
        nodes["types_count"] = nodes["id"].map(lambda i: len(tmap.get(i, []))).fillna(0).astype(int)
        nodes["primary_type"] = nodes["id"].map(lambda i: (tmap.get(i, [""] ) or [""])[0])

        if make_type_indicator_cols:
            # collect all unique types
            all_types = sorted({t for lst in tmap.values() for t in lst})
            def _safe(s: str) -> str:
                s = re.sub(r"\W+", "_", s).strip("_")
                return s or "type"
            for t in all_types:
                col = f"type__{_safe(t)}"
                nodes[col] = nodes["id"].map(lambda i: t in tmap.get(i, []))

    # ---------- Write files ----------
    if graphml_path:
        G = nx.MultiDiGraph() if directed else nx.MultiGraph()
        for _, r in nodes.iterrows():
            attrs = r.drop(labels=["id"]).to_dict()
            G.add_node(r["id"], **attrs)
        for _, r in edges.iterrows():
            G.add_edge(r["source"], r["target"], key=r["type"], type=str(r["type"]), weight=float(r["weight"]))
        nx.write_graphml(G, graphml_path)
        print(f"[info] GraphML written to: {graphml_path}")

    if nodes_csv_path and edges_csv_path:
        nodes.to_csv(nodes_csv_path, index=False)
        edges.to_csv(edges_csv_path, index=False)
        print(f"[info] Wrote CSVs:\n  {nodes_csv_path}\n  {edges_csv_path}")

    return nodes, edges



# name suffix only if top-k was specified
suffix = f".topk{TOP_K}" if TOP_K is not None else ""

nodes_df, edges_df = union_multiplex_edges_for_cytoscape(
    emb_sim_df, coexpr_corr_df, selected,
    nameA="embeddings", nameB="coexpr",
    directed=False,
    keep_self_loops=False,
    drop_pairs_where_both_zero=False,
    output_dir="for_cytoscape",
    write_graphml_path=f"union{suffix}.graphml",         # e.g., union.topk5.graphml
    write_csv_prefix=f"embeddings_union{suffix}",        # e.g., embeddings_union.topk5_{nodes,edges}.csv
    node_types=node_types,
    type_separator="|",
    make_type_indicator_cols=True,
    top_k_per_node=TOP_K                                  # <-- script-level arg plumbed through
)