import time
import pandas as pd
import numpy as np
from sklearn.metrics.pairwise import cosine_similarity
import argparse
from pathlib import Path

# ---- helper: L2-normalize a vector ----
def _l2_normalize(vec, eps=1e-12):
    vec = np.asarray(vec, dtype=float)
    n = np.linalg.norm(vec)
    return vec / n if n > eps else vec * 0.0


def load_and_process_data(mapping_file_path, embeddings_path, coexpr_corr_path, skip_embeddings=False, skip_coexpr=False):
    """Load and process all necessary data files."""
    
    # Load mapping table
    print("Loading mapping table...")
    mapping_df = pd.read_csv(mapping_file_path)
    print(f"Mapping table shape: {mapping_df.shape}")
    
    # Create mapping dict
    gene_col = "gene_id"
    ensembl_col = "Ensembl_gene_id"
    mapper = mapping_df.set_index(ensembl_col)[gene_col].to_dict()
    
    # Load embeddings
    if not skip_embeddings:
        print("\nLoading embeddings...")
        time_start = time.time()
        embeddings = pd.read_csv(embeddings_path, index_col=0)
        # Translate ensembl to gene name
        embeddings = embeddings.rename(index=mapper)
        # Keep only numeric columns
        X = embeddings.select_dtypes(include="number").astype(float)
        time_end = time.time()
        print(f"Time taken to load embeddings: {time_end - time_start:.2f} seconds")
        print(f"Embeddings shape (genes x dims): {X.shape}")
        
        # Compute cosine similarity between embeddings (gene x gene)
        print("\nComputing cosine similarity for embeddings...")
        S = cosine_similarity(X)
        emb_sim_df = pd.DataFrame(S, index=X.index, columns=X.index)
        print(f"Embedding similarity matrix shape: {emb_sim_df.shape}")
    else:
        print("\nSkipping embeddings per --skip-embeddings flag")
        emb_sim_df = pd.DataFrame()
        X = pd.DataFrame()
    
    # Load co-expression correlation matrix
    if not skip_coexpr:
        print("\nLoading co-expression correlation matrix from path:")
        print(f"{coexpr_corr_path}")
        time_start = time.time()
        coexpr_corr = np.load(coexpr_corr_path, allow_pickle=True)["corr"]
        coexpr_genes = np.load(coexpr_corr_path, allow_pickle=True)["genes"]
        time_end = time.time()
        
        # Construct a pandas dataframe
        coexpr_corr_df = pd.DataFrame(coexpr_corr, index=coexpr_genes, columns=coexpr_genes)
        # Translate ensembl to gene name
        coexpr_corr_df = coexpr_corr_df.rename(index=mapper)
        coexpr_corr_df = coexpr_corr_df.rename(columns=mapper)
        print(f"Time taken to load co-expression correlation matrix: {time_end - time_start:.2f} seconds")
        print(f"Co-expression correlation matrix shape: {coexpr_corr_df.shape}")
    else:
        print("\nSkipping co-expression per --skip-coexpr flag")
        coexpr_corr_df = pd.DataFrame()
    
    # Return both the pairwise cosine sim matrix and the raw embedding table X
    return emb_sim_df, coexpr_corr_df, mapping_df, X


def _build_query_vector(genes, embeddings_df, mode):
    """
    Build a query vector from one or two genes.

    L2-normalizes each component BEFORE combining for 'sum'/'mean_l2' to
    avoid norm dominance. Single-gene returns the raw embedding (original behavior).

    genes: list[str] of length 1 or 2
    embeddings_df: DataFrame with gene index and embedding columns
    mode: 'single' | 'sum' | 'mean_l2'
    """
    if len(genes) == 1:
        g = genes[0]
        if g not in embeddings_df.index:
            raise ValueError(f"Gene '{g}' not found in embeddings index")
        return embeddings_df.loc[g].to_numpy()

    # Two-gene composite
    g1, g2 = genes
    missing = [g for g in (g1, g2) if g not in embeddings_df.index]
    if missing:
        raise ValueError(f"Gene(s) not found in embeddings index: {', '.join(missing)}")

    v1 = embeddings_df.loc[g1].to_numpy(dtype=float)
    v2 = embeddings_df.loc[g2].to_numpy(dtype=float)

    # L2-normalize components BEFORE combining
    v1 = _l2_normalize(v1)
    v2 = _l2_normalize(v2)

    if mode == "sum":
        return v1 + v2
    elif mode == "mean_l2":
        q = (v1 + v2) / 2.0
        return _l2_normalize(q)
    else:
        # 'single' with two genes doesn't make sense; fallback to 'sum' for convenience
        print("Warning: --embed-combine 'single' provided with two genes. Defaulting to 'sum'.")
        return v1 + v2


def _cosine_similarities_to_query_vec(embeddings_df, qvec):
    """
    Compute cosine similarity of every row in embeddings_df to the query vector qvec.
    Returns a pandas Series indexed by gene.
    """
    sims = cosine_similarity(embeddings_df.to_numpy(), qvec.reshape(1, -1)).ravel()
    return pd.Series(sims, index=embeddings_df.index)


def _load_external_vector(vec_path, embeddings_df):
    """
    Load an external embedding vector from CSV or NPY and align to embeddings_df columns if possible.
    Returns a 1D numpy array.
    """
    if vec_path is None:
        return None

    p = Path(vec_path)
    if not p.exists():
        raise FileNotFoundError(f"External vector file not found: {vec_path}")

    dims = embeddings_df.shape[1]
    cols = list(embeddings_df.columns)

    if p.suffix.lower() in [".npy", ".npz"]:
        arr = np.load(p)
        if isinstance(arr, np.lib.npyio.NpzFile):
            for key in ("vector", "emb", "q", "data"):
                if key in arr:
                    arr = arr[key]
                    break
            else:
                first_key = list(arr.keys())[0]
                arr = arr[first_key]
        arr = np.asarray(arr).reshape(-1)
        if arr.size != dims:
            raise ValueError(f"External vector dims ({arr.size}) != embeddings dims ({dims})")
        return arr

    # CSV / TSV
    df = pd.read_csv(p)
    # Case A: columns match embedding columns names (order them)
    if set(cols).issubset(set(df.columns)):
        row = df[cols].iloc[0].to_numpy(dtype=float)
        return row
    # Case B: single row with exactly dims values
    if df.shape[0] == 1 and df.shape[1] == dims:
        return df.iloc[0].to_numpy(dtype=float)
    # Case C: single column of length dims
    if df.shape[1] == 1 and df.shape[0] == dims:
        return df.iloc[:, 0].to_numpy(dtype=float)

    raise ValueError(
        f"Could not interpret external vector file '{vec_path}'. "
        f"Provide a CSV/NPY with {dims} values or a CSV with columns matching the embedding dimensions."
    )


def rank_genes_by_similarity(query_genes, emb_sim_df, coexpr_corr_df, embeddings_df,
                             embed_combine_mode="single", output_dir="gene_rankings", output_prefix="",
                             external_qvec=None, external_label=None):
    """
    Rank all genes by their similarity to a query:
      • single gene
      • two-gene composite
      • gene + external vector
      • external vector only (no gene)

    Parameters
    ----------
    query_genes : list[str]
        0, 1, or 2 gene names
    emb_sim_df : pd.DataFrame
        Embedding cosine similarity matrix (gene x gene). Optional if using composite queries.
    coexpr_corr_df : pd.DataFrame
        Co-expression correlation matrix
    embeddings_df : pd.DataFrame
        Raw embeddings (gene x dims), used for composite queries or fallback
    embed_combine_mode : str
        'single' (one gene), 'sum', or 'mean_l2' (two genes or gene+external vector)
    output_dir : str
        Directory to save output files
    output_prefix : str
        Prefix for output filenames (applies to co-expression output only)
    external_qvec : np.ndarray | None
        Optional external embedding vector to combine with a single gene or to use alone
    external_label : str | None
        Name/label for the external embedding vector (used in filenames)
    
    Returns
    -------
    tuple : (embedding_rankings_df or None, coexpr_rankings_df or None)
    """
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    num_genes = len(query_genes)
    two_gene = (num_genes == 2)
    using_external = (external_qvec is not None)
    external_only = using_external and (num_genes == 0)
    qlabel = None

    # ----- Embedding rankings -----
    emb_rankings = None
    if not embeddings_df.empty:
        if two_gene or using_external:
            mode = embed_combine_mode if embed_combine_mode in ("sum", "mean_l2") else "sum"

            if two_gene:
                qvec = _build_query_vector(query_genes, embeddings_df, mode)
                qlabel = f"{query_genes[0]}+{query_genes[1]}_{mode}"
            elif external_only:
                qvec = np.asarray(external_qvec, dtype=float)
                ext_name = external_label if external_label else "external"
                qlabel = f"{ext_name}"
            else:
                # single gene + external vector
                qg = query_genes[0]
                if qg not in embeddings_df.index:
                    raise ValueError(f"Gene '{qg}' not found in embeddings index")
                v_gene = embeddings_df.loc[qg].to_numpy(dtype=float)
                evec = np.asarray(external_qvec, dtype=float)

                # L2-normalize components BEFORE combining
                v_gene = _l2_normalize(v_gene)
                evec = _l2_normalize(evec)

                if mode == "sum":
                    qvec = v_gene + evec
                elif mode == "mean_l2":
                    q = (v_gene + evec) / 2.0
                    qvec = _l2_normalize(q)
                else:
                    qvec = v_gene + evec
                ext_name = external_label if external_label else "external"
                qlabel = f"{qg}+{ext_name}_{mode}"

            # Print & save the query vector (even if it came from file, for consistency)
            print(f"\nComposite/query vector [{qlabel}] (length {qvec.size}):")
            print(np.array2string(qvec, precision=6, suppress_small=False, threshold=qvec.size))

            qvec_df = pd.DataFrame(qvec.reshape(1, -1), columns=embeddings_df.columns)
            qvec_file = output_path / f"{output_prefix}_query_vector.csv"
            qvec_df.to_csv(qvec_file, index=False)
            print(f"Query vector saved to: {qvec_file}")

            # Rank by cosine similarity
            sims = _cosine_similarities_to_query_vec(embeddings_df, qvec).sort_values(ascending=False)
            emb_rankings = pd.DataFrame({
                'gene': sims.index,
                'embedding_cosine_similarity': sims.values,
                'rank': range(1, len(sims) + 1)
            })
            emb_output_file = output_path / f"{output_prefix}_embedding_similarity_ranked.csv"
            emb_rankings.to_csv(emb_output_file, index=False)
            print(f"\nEmbedding similarity rankings saved to: {emb_output_file}")
            print(f"Top 10 genes by embedding similarity for [{qlabel}]:")
            print(emb_rankings.head(10).to_string(index=False))

        else:
            # Single gene only
            qg = query_genes[0]
            if (not emb_sim_df.empty) and (qg in emb_sim_df.index):
                sims = emb_sim_df[qg].copy().sort_values(ascending=False)
            else:
                if qg not in embeddings_df.index:
                    print(f"Warning: {qg} not found in embeddings; skipping embedding rankings.")
                    sims = None
                else:
                    qvec = embeddings_df.loc[qg].to_numpy()
                    sims = _cosine_similarities_to_query_vec(embeddings_df, qvec).sort_values(ascending=False)

            if sims is not None:
                qlabel = qg

                # --- save the single-gene embedding vector ---
                if qg in embeddings_df.index:
                    single_vec = embeddings_df.loc[qg].to_numpy()
                    single_vec_df = pd.DataFrame(single_vec.reshape(1, -1), columns=embeddings_df.columns)
                    single_vec_file = output_path / f"{qlabel}_query_vector.csv"
                    single_vec_df.to_csv(single_vec_file, index=False)
                    print(f"Single-gene query vector saved to: {single_vec_file}")
                # --- End ---

                emb_rankings = pd.DataFrame({
                    'gene': sims.index,
                    'embedding_cosine_similarity': sims.values,
                    'rank': range(1, len(sims) + 1)
                })
                emb_output_file = output_path / f"{qlabel}_embedding_similarity_ranked.csv"
                emb_rankings.to_csv(emb_output_file, index=False)
                print(f"\nEmbedding similarity rankings saved to: {emb_output_file}")
                print(f"Top 10 genes by embedding similarity:")
                print(emb_rankings.head(10).to_string(index=False))
    else:
        print("\nSkipping embedding similarity rankings per --skip-embeddings flag")

    # ----- Co-expression rankings -----
    coexpr_rankings = None
    if not coexpr_corr_df.empty:
        if num_genes >= 1:
            q_for_coexpr = query_genes[0]
            if two_gene or using_external:
                print(f"\nNote: co-expression ranking uses the first gene only: {q_for_coexpr}")
            if q_for_coexpr in coexpr_corr_df.index:
                coexpr_correlations = coexpr_corr_df[q_for_coexpr].copy().sort_values(ascending=False)
                coexpr_rankings = pd.DataFrame({
                    'gene': coexpr_correlations.index,
                    'coexpression_correlation': coexpr_correlations.values,
                    'rank': range(1, len(coexpr_correlations) + 1)
                })
                coexpr_label = q_for_coexpr if output_prefix == "" else f"{q_for_coexpr}_{output_prefix}"
                coexpr_output_file = output_path / f"{coexpr_label}_coexpression_correlation_ranked.csv"
                coexpr_rankings.to_csv(coexpr_output_file, index=False)
                print(f"\nCo-expression correlation rankings saved to: {coexpr_output_file}")
                print(f"Top 10 genes by co-expression correlation:")
                print(coexpr_rankings.head(10).to_string(index=False))
            else:
                print(f"\nSkipping co-expression correlation rankings: {q_for_coexpr} not found in co-expression matrix")
        else:
            print("\nSkipping co-expression correlation rankings: no query gene provided")
    else:
        print("\nSkipping co-expression correlation rankings per --skip-coexpr flag")
    
    return emb_rankings, coexpr_rankings


def main():
    """Main function to run the gene similarity ranking."""
    
    parser = argparse.ArgumentParser(description='Rank genes by similarity to a query gene, two-gene composite, gene+external vector, or external vector alone.')
    # Backward-compat single-gene flag (deprecated)
    parser.add_argument('--gene', type=str, help='[Deprecated] Use --genes instead.')
    # One or two genes (optional)
    parser.add_argument('--genes', nargs='+', type=str, required=False,
                       help='Zero, one, or two gene names (e.g., --genes ifnphi1 or --genes ifnphi1 ifnphi2)')
    parser.add_argument('--embed-combine', type=str, choices=['single', 'sum', 'mean_l2'], default='single',
                       help="How to combine two gene embeddings (or gene+external vector). Ignored for external-only.")
    parser.add_argument('--mapping', type=str, 
                       default="/path/to/Geneformer_RQfork/data/RNAquarium/xref/drerio_z12-geneid_to_z11-ENSDARG.csv",
                       help='Path to gene ID mapping file')
    parser.add_argument('--embeddings', type=str,
                       default="/path/to/Geneformer_RQfork/embeddings/extracted_embeddings/RNAquarium_second_to_last_layer.csv",
                       help='Path to embeddings CSV file')
    parser.add_argument('--coexpr', type=str,
                       default="/path/to/Geneformer_RQfork/embeddings/coexpression/corr_matrix/corr_matrix.npz",
                       help='Path to co-expression correlation matrix')
    parser.add_argument('--output-dir', type=str, default='gene_rankings',
                       help='Output directory for ranking files')
    parser.add_argument('--output-prefix', type=str, default='',
                       help='Output prefix for ranking files (applied to co-expression output only)')
    parser.add_argument('--skip-embeddings', action='store_true',
                       help='Skip loading embeddings and embedding similarity computation')
    parser.add_argument('--skip-coexpr', action='store_true',
                       help='Skip loading co-expression correlation matrix')

    # External vector support (CSV or NPY). Can be used alone or with one gene.
    parser.add_argument('--query-vector-path', type=str, default=None,
                        help='Path to a CSV/NPY external embedding vector (can be used without genes)')
    parser.add_argument('--query-vector-name', type=str, default=None,
                        help='Name/label for the external vector (used in output filenames). Defaults to file stem.')

    args = parser.parse_args()

    # Backward compat: map --gene to --genes
    if args.genes is None and args.gene is not None:
        print("Notice: --gene is deprecated; please use --genes. Proceeding with single-gene query.")
        args.genes = [args.gene]
    if args.genes is None:
        args.genes = []
    if len(args.genes) > 2:
        raise SystemExit("Please provide at most two genes for --genes.")

    using_external = args.query_vector_path is not None

    # Require at least one of genes or external vector
    if (not args.genes) and (not using_external):
        raise SystemExit("Provide --genes <gene1> [gene2] and/or --query-vector-path PATH.")

    # If external vector is provided, allow at most one gene
    if using_external and len(args.genes) > 1:
        raise SystemExit("When using --query-vector-path, provide at most one gene with --genes.")

    # Defaults for combine mode
    if len(args.genes) == 2 and args.embed_combine == 'single':
        print("Two genes provided with --embed-combine 'single'; defaulting to 'sum'.")
        args.embed_combine = 'sum'
    if using_external and len(args.genes) == 1 and args.embed_combine == 'single':
        print("Gene + external vector with --embed-combine 'single'; defaulting to 'sum'.")
        args.embed_combine = 'sum'

    # Derive external label default from file stem if not provided
    if using_external and (args.query_vector_name is None or args.query_vector_name.strip() == ""):
        args.query_vector_name = Path(args.query_vector_path).stem

    print(f"Query genes: {args.genes if args.genes else '[]'}")
    print(f"Embedding combine mode: {args.embed_combine}")
    print(f"Output directory: {args.output_dir}")
    print(f"Output prefix: {args.output_prefix}")
    if using_external:
        print(f"External vector path: {args.query_vector_path}")
        print(f"External vector name: {args.query_vector_name}")
    print("=" * 60)
    
    # Load and process data
    emb_sim_df, coexpr_corr_df, mapping_df, embeddings_df = load_and_process_data(
        args.mapping, 
        args.embeddings, 
        args.coexpr,
        skip_embeddings=args.skip_embeddings,
        skip_coexpr=args.skip_coexpr
    )

    # Load external vector if provided
    external_qvec = None
    if using_external:
        if embeddings_df.empty:
            raise SystemExit("Embeddings are required when using --query-vector-path (remove --skip-embeddings).")
        external_qvec = _load_external_vector(args.query_vector_path, embeddings_df)
    
    print("\n" + "=" * 60)
    print("Generating similarity rankings...")
    print("=" * 60)
    
    # Rank genes by similarity
    emb_rankings, coexpr_rankings = rank_genes_by_similarity(
        args.genes, 
        emb_sim_df, 
        coexpr_corr_df,
        embeddings_df,
        embed_combine_mode=args.embed_combine,
        output_dir=args.output_dir,
        output_prefix=args.output_prefix,
        external_qvec=external_qvec,
        external_label=args.query_vector_name
    )
    
    # Print summary statistics
    print("\n" + "=" * 60)
    print("Summary Statistics")
    print("=" * 60)
    
    if emb_rankings is not None:
        if (len(args.genes) == 1 and not using_external):
            q_label = args.genes[0]
            print(f"\nEmbedding similarity to {q_label}:")
            print(f"  Self-similarity: {emb_rankings.iloc[0]['embedding_cosine_similarity']:.4f}")
        else:
            if using_external and len(args.genes) == 0:
                q_label = f"{args.query_vector_name}"
                prefix_text = "external vector"
            elif using_external and len(args.genes) == 1:
                q_label = f"{args.genes[0]}+{args.query_vector_name}_{args.embed_combine}"
                prefix_text = "gene + external vector"
            else:
                q_label = f"{args.genes[0]}+{args.genes[1]}_{args.embed_combine}"
                prefix_text = "composite query"
            top_gene = emb_rankings.iloc[0]['gene']
            top_sim = emb_rankings.iloc[0]['embedding_cosine_similarity']
            print(f"\nEmbedding similarity for {prefix_text} [{q_label}]:")
            print(f"  Top match: {top_gene} (cosine={top_sim:.4f})")
        print(f"  Mean similarity: {emb_rankings['embedding_cosine_similarity'].mean():.4f}")
        print(f"  Median similarity: {emb_rankings['embedding_cosine_similarity'].median():.4f}")
        print(f"  Min similarity: {emb_rankings['embedding_cosine_similarity'].min():.4f}")
        print(f"  Max similarity: {emb_rankings['embedding_cosine_similarity'].max():.4f}")
    
    if coexpr_rankings is not None:
        q_ce = args.genes[0]
        print(f"\nCo-expression correlation with {q_ce}:")
        print(f"  Self-correlation: {coexpr_rankings.iloc[0]['coexpression_correlation']:.4f}")
        print(f"  Mean correlation: {coexpr_rankings['coexpression_correlation'].mean():.4f}")
        print(f"  Median correlation: {coexpr_rankings['coexpression_correlation'].median():.4f}")
        print(f"  Min correlation: {coexpr_rankings['coexpression_correlation'].min():.4f}")
        print(f"  Max correlation: {coexpr_rankings['coexpression_correlation'].max():.4f}")
    
    print("\nProcess complete!")

if __name__ == "__main__":
    main()
