import time
import pandas as pd
import numpy as np
from sklearn.metrics.pairwise import cosine_similarity
import argparse
from pathlib import Path
import sys
import random

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

def _parse_gene_list_arg(list_arg, file_arg):
    """
    Accept genes via --genes G1 G2 ... and/or --genes-file path (txt/csv with one per line or comma/space-separated).
    Returns a list of unique, order-preserved gene names (as strings).
    """
    out = []
    if list_arg:
        for tok in list_arg:
            parts = [p for p in str(tok).replace(",", " ").split() if p]
            out.extend(parts)
    if file_arg:
        p = Path(file_arg)
        if not p.exists():
            raise FileNotFoundError(f"Gene list file not found: {file_arg}")
        txt = p.read_text()
        parts = [p for p in txt.replace(",", " ").replace("\t", " ").split() if p]
        out.extend(parts)
    # Deduplicate preserving order
    seen = set()
    uniq = []
    for g in out:
        if g not in seen:
            seen.add(g)
            uniq.append(g)
    return uniq

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
# Mean vector builders
# -----------------------------
def build_mean_vector_all_at_once(genes, embeddings_df):
    """Arithmetic mean across selected genes, then L2-normalize once."""
    sub = embeddings_df.loc[genes].to_numpy(dtype=float)
    mean_vec = sub.mean(axis=0)
    return _l2_normalize(mean_vec)

def build_mean_vector_sequential(genes, embeddings_df):
    """
    Sequential mean with L2-normalization at each step:
      v_0 = 0
      for i, g_i:
         v_i = L2norm( ( (i-1)/i ) * v_{i-1} + (1/i) * emb(g_i) )
    This differs from all-at-once because v_{i-1} is unit-length at each step.
    """
    d = embeddings_df.shape[1]
    v = np.zeros(d, dtype=float)
    for i, g in enumerate(genes, start=1):
        e = embeddings_df.loc[g].to_numpy(dtype=float)
        v = ((i - 1) / i) * v + (1.0 / i) * e
        v = _l2_normalize(v)
    return v

def build_mean_vector(genes, embeddings_df, strategy):
    if strategy == "all":
        return build_mean_vector_all_at_once(genes, embeddings_df)
    elif strategy == "sequential":
        return build_mean_vector_sequential(genes, embeddings_df)
    else:
        raise ValueError(f"Unknown combine strategy: {strategy}")

# -----------------------------
# Similarity computation
# -----------------------------
def cosine_to_targets(qvec, embeddings_df, target_genes):
    present = [g for g in target_genes if g in embeddings_df.index]
    missing = [g for g in target_genes if g not in embeddings_df.index]
    if missing:
        print(f"Warning: {len(missing)} target gene(s) not found in embeddings and will be NaN: {', '.join(missing[:10])}{' ...' if len(missing)>10 else ''}")

    if present:
        T = embeddings_df.loc[present].to_numpy(dtype=float)
        Tn = _rowwise_l2_normalize(T)
        qn = _l2_normalize(qvec)
        # Cosine between each target row and query vector
        sims_present = Tn.dot(qn)
        ser_present = pd.Series(sims_present, index=present)
    else:
        ser_present = pd.Series(dtype=float)

    # Combine present + NaNs for missing, preserving original order
    sims = []
    for g in target_genes:
        sims.append(ser_present[g] if g in ser_present.index else np.nan)
    return pd.Series(sims, index=target_genes)

# -----------------------------
# Modes
# -----------------------------
def run_mode_random(n, m, strategy, embeddings_df, target_genes, seed=None):
    rng = np.random.default_rng(seed) if seed is not None else np.random.default_rng()
    pool = list(embeddings_df.index)
    if n < 1:
        raise SystemExit("--n must be >= 1")
    if n > len(pool):
        raise SystemExit(f"--n ({n}) exceeds number of available genes in embeddings ({len(pool)})")

    mean_vectors = []
    sims_rows = []
    picked_records = []

    for i in range(1, m + 1):
        picked = rng.choice(pool, size=n, replace=False).tolist()
        qvec = build_mean_vector(picked, embeddings_df, strategy=strategy)
        sims = cosine_to_targets(qvec, embeddings_df, target_genes)

        mean_vectors.append(qvec)
        sims_rows.append(sims)
        picked_records.append({"iteration": i, "genes": ",".join(picked)})

        if i <= 3 or i == m:
            print(f"[Iter {i}] first 5 picked: {picked[:5]} ... | first 5 sims: {sims.head(5).round(4).to_dict()}")

    mean_vectors_df = pd.DataFrame(mean_vectors, index=[f"iter_{i}" for i in range(1, m + 1)], columns=embeddings_df.columns)
    sims_df = pd.DataFrame(sims_rows, index=[f"iter_{i}" for i in range(1, m + 1)])
    picked_df = pd.DataFrame(picked_records)

    return mean_vectors_df, sims_df, picked_df

def run_mode_specified(genes, strategy, embeddings_df, target_genes):
    present = [g for g in genes if g in embeddings_df.index]
    missing = [g for g in genes if g not in embeddings_df.index]
    if missing:
        print(f"Warning: {len(missing)} specified gene(s) not found and will be ignored: {', '.join(missing[:10])}{' ...' if len(missing)>10 else ''}")
    if len(present) == 0:
        raise SystemExit("None of the specified genes were found in the embeddings index.")
    qvec = build_mean_vector(present, embeddings_df, strategy=strategy)
    sims = cosine_to_targets(qvec, embeddings_df, target_genes)

    mean_vectors_df = pd.DataFrame([qvec], index=["specified_set"], columns=embeddings_df.columns)
    sims_df = pd.DataFrame([sims], index=["specified_set"])
    picked_df = pd.DataFrame([{"genes": ",".join(present)}])
    return mean_vectors_df, sims_df, picked_df

# -----------------------------
# Main
# -----------------------------
def main():
    parser = argparse.ArgumentParser(
        description="Compute mean embedding vectors (random or specified gene sets) and cosine similarity to a target gene list."
    )

    # Mode & behavior
    parser.add_argument("--mode", choices=["random", "specified"], required=True,
                        help="Mode 1 = 'random': sample n genes m times; Mode 2 = 'specified': use provided n genes once.")
    parser.add_argument("--combine", choices=["sequential", "all"], default="all",
                        help="How to build the mean vector: 'sequential' = sequential mean with per-step L2; 'all' = average all then L2 once.")
    parser.add_argument("--n", type=int, required=True, help="Number of genes to sample/use.")
    parser.add_argument("--m", type=int, default=1, help="[random mode only] Number of repetitions.")
    parser.add_argument("--seed", type=int, default=None, help="[random mode] Random seed for reproducibility.")

    # Specified genes (mode 2)
    parser.add_argument("--genes", nargs="*", help="[specified mode] Gene symbols, space/comma-separated.")
    parser.add_argument("--genes-file", type=str, default=None, help="[specified mode] File with gene symbols.")

    # Targets to score against
    parser.add_argument("--targets", nargs="+", required=True, help="Target gene symbols to compute cosine similarity to (space/comma-separated).")
    parser.add_argument("--targets-file", type=str, default=None, help="Optional file with target gene symbols.")

    # Paths (defaults preserved from your original script)
    parser.add_argument('--mapping', type=str, 
        default="/path/to/Geneformer_RQfork/data/RNAquarium/xref/drerio_z12-geneid_to_z11-ENSDARG.csv",
        help='Path to gene ID mapping file')
    parser.add_argument('--embeddings', type=str,
        default="/path/to/Geneformer_RQfork/embeddings/extracted_embeddings/RNAquarium_second_to_last_layer.csv",
        help='Path to embeddings CSV file')

    # Output
    parser.add_argument('--output-dir', type=str, default='gene_rankings',
        help='Output directory')
    parser.add_argument('--output-prefix', type=str, default='',
        help='Prefix for output filenames')

    args = parser.parse_args()

    # Load data (kept preprocessing & default paths)
    mapping_df, embeddings_df = load_and_process_data(args.mapping, args.embeddings)

    # Parse targets
    target_genes = _parse_gene_list_arg(args.targets, args.targets_file)
    if len(target_genes) == 0:
        raise SystemExit("Please provide at least one target gene via --targets and/or --targets-file.")
    print(f"\nNumber of target genes provided: {len(target_genes)}")

    # Run modes
    print("\n" + "=" * 60)
    print(f"Running mode: {args.mode} | combine={args.combine} | n={args.n} | m={args.m if args.mode=='random' else 1}")
    print("=" * 60)

    if args.mode == "random":
        mean_vectors_df, sims_df, picked_df = run_mode_random(
            n=args.n, m=args.m, strategy=args.combine,
            embeddings_df=embeddings_df,
            target_genes=target_genes, seed=args.seed
        )
        mode_tag = "model_random"
    else:
        specified_genes = _parse_gene_list_arg(args.genes, args.genes_file)
        if len(specified_genes) == 0:
            raise SystemExit("[specified mode] Provide genes via --genes and/or --genes-file.")
        if args.n is not None and args.n > 0 and len(specified_genes) != args.n:
            print(f"Notice: --n={args.n} but {len(specified_genes)} genes were provided; proceeding with the provided list.")
        mean_vectors_df, sims_df, picked_df = run_mode_specified(
            genes=specified_genes, strategy=args.combine,
            embeddings_df=embeddings_df, target_genes=target_genes
        )
        mode_tag = "mode2_specified"

    # Save outputs
    outdir = Path(args.output_dir)
    outdir.mkdir(parents=True, exist_ok=True)

    prefix = f"{args.output_prefix}_" if args.output_prefix else ""
    base = f"{prefix}{mode_tag}_{args.combine}_n{args.n}"
    if args.mode == "random":
        base += f"_m{args.m}"

    mean_vec_path = outdir / f"{base}_mean_vectors.csv"
    sims_path = outdir / f"{base}_cosine_to_targets.csv"
    picked_path = outdir / f"{base}_picked_genes.csv"

    mean_vectors_df.to_csv(mean_vec_path)
    sims_df.to_csv(sims_path)
    picked_df.to_csv(picked_path, index=False)

    print("\n" + "=" * 60)
    print("Outputs")
    print("=" * 60)
    print(f"Mean vector(s): {mean_vec_path}")
    print(f"Cosine to targets: {sims_path}")
    print(f"Picked genes: {picked_path}")
    print("\nDone.")

if __name__ == "__main__":
    main()
