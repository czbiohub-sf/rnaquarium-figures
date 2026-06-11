import time
import pandas as pd
import numpy as np
import argparse
from pathlib import Path

# -----------------------------
# Helpers
# -----------------------------
def _l2_normalize(vec, eps=1e-12):
    vec = np.asarray(vec, dtype=float)
    n = np.linalg.norm(vec)
    return vec / n if n > eps else vec * 0.0

def _rowwise_l2_normalize(M, eps=1e-12):
    M = np.asarray(M, dtype=float)
    norms = np.linalg.norm(M, axis=1, keepdims=True)
    norms = np.where(norms > eps, norms, np.inf)  # avoid divide by 0
    return M / norms

# -----------------------------
# Data loading (keeps defaults + mapping)
# -----------------------------
def load_and_process_data(mapping_file_path, embeddings_path):
    """Load mapping table and embeddings, keep numeric columns, and rename indices via mapping."""
    # Load mapping
    print("Loading mapping table...")
    mapping_df = pd.read_csv(mapping_file_path)
    print(f"Mapping table shape: {mapping_df.shape}")

    # Create mapping dict
    gene_col = "gene_id"
    ensembl_col = "Ensembl_gene_id"
    mapper = mapping_df.set_index(ensembl_col)[gene_col].to_dict()

    # Load embeddings
    print("\nLoading embeddings...")
    t0 = time.time()
    embeddings = pd.read_csv(embeddings_path, index_col=0)
    # Translate Ensembl -> gene symbol using mapping
    embeddings = embeddings.rename(index=mapper)
    # Keep only numeric columns
    X = embeddings.select_dtypes(include="number").astype(float)
    t1 = time.time()
    print(f"Time taken to load embeddings: {t1 - t0:.2f} s")
    print(f"Embeddings shape (genes x dims): {X.shape}")

    # Drop rows with all-NaN (just in case)
    X = X.dropna(how="all")

    return mapping_df, X

# -----------------------------
# Pairwise cosine for random pairs
# -----------------------------
def compute_pairwise_cosine_random_pairs(m, embeddings_df, seed=None):
    """
    Sample m random unordered pairs of distinct genes (with replacement across iterations),
    compute cosine similarity between the two genes in each pair, and return a dataframe.
    """
    rng = np.random.default_rng(seed)
    pool = list(embeddings_df.index)
    if len(pool) < 2:
        raise SystemExit("Embeddings must contain at least 2 genes to form a pair.")

    # Pre-normalize all rows once for efficient dot products
    Xn = _rowwise_l2_normalize(embeddings_df.to_numpy(dtype=float))
    Xn_df = pd.DataFrame(Xn, index=embeddings_df.index, columns=embeddings_df.columns)

    records = []
    for i in range(1, m + 1):
        g1, g2 = rng.choice(pool, size=2, replace=False).tolist()
        v1 = Xn_df.loc[g1].to_numpy()
        v2 = Xn_df.loc[g2].to_numpy()
        # Cosine similarity = dot of L2-normalized vectors
        cos = float(np.dot(v1, v2))

        if i <= 3 or i == m:
            print(f"[Iter {i}] pair: ({g1}, {g2}) | cosine={cos:.4f}")

        records.append({
            "iteration": i,
            "gene_a": g1,
            "gene_b": g2,
            "cosine_similarity": cos
        })

    return pd.DataFrame(records)

# -----------------------------
# Main
# -----------------------------
def main():
    parser = argparse.ArgumentParser(
        description="Sample random gene pairs (n=2) and compute pairwise cosine similarity."
    )

    # Core sampling params
    parser.add_argument("--n", type=int, default=2,
                        help="Number of genes per sample; must be 2 for pairwise cosine.")
    parser.add_argument("--m", type=int, default=1,
                        help="Number of random pairs to sample.")
    parser.add_argument("--seed", type=int, default=None,
                        help="Random seed for reproducibility.")

    # Paths (defaults preserved)
    parser.add_argument('--mapping', type=str, 
        default="/path/to/Geneformer_RQfork/data/RNAquarium/xref/drerio_z12-geneid_to_z11-ENSDARG.csv",
        help='Path to gene ID mapping file')
    parser.add_argument('--embeddings', type=str,
        default="/path/to/Geneformer_RQfork/embeddings/extracted_embeddings/RNAquarium_second_to_last_layer.csv",
        help='Path to embeddings CSV file')

    # Output
    parser.add_argument('--output-dir', type=str, default='gene_pair_cosines',
        help='Output directory')
    parser.add_argument('--output-prefix', type=str, default='',
        help='Prefix for output filename')

    args = parser.parse_args()

    if args.n != 2:
        raise SystemExit(f"--n must be 2 for pairwise cosine computation; got {args.n}")

    # Load data
    _, embeddings_df = load_and_process_data(args.mapping, args.embeddings)

    # Compute random pairwise cosines
    print("\n" + "=" * 60)
    print(f"Sampling {args.m} random gene pairs (n=2) | seed={args.seed}")
    print("=" * 60)

    pairs_df = compute_pairwise_cosine_random_pairs(m=args.m, embeddings_df=embeddings_df, seed=args.seed)

    # Save outputs
    outdir = Path(args.output_dir)
    outdir.mkdir(parents=True, exist_ok=True)

    prefix = f"{args.output_prefix}_" if args.output_prefix else ""
    base = f"{prefix}pairwise_cosine_n2_m{args.m}"
    out_path = outdir / f"{base}.csv"

    pairs_df.to_csv(out_path, index=False)

    print("\n" + "=" * 60)
    print("Output")
    print("=" * 60)
    print(f"Pairwise cosine results: {out_path}")
    print("\nDone.")

if __name__ == "__main__":
    main()
