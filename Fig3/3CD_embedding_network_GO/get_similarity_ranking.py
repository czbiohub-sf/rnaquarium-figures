import time
import pandas as pd
import numpy as np
from sklearn.metrics.pairwise import cosine_similarity
import argparse
from pathlib import Path

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

        # TODO save a copy of the embeddings with gene names as index
        embeddings.to_csv(f"{embeddings_path.removesuffix('.csv')}_with_genenames.csv")
        
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

    v1 = embeddings_df.loc[g1].to_numpy()
    v2 = embeddings_df.loc[g2].to_numpy()

    if mode == "sum":
        return v1 + v2
    elif mode == "mean_l2":
        q = (v1 + v2) / 2.0
        norm = np.linalg.norm(q)
        return q / norm if norm > 0 else q * 0.0
    else:
        # 'single' with two genes doesn't make sense; fallback to 'sum' for convenience
        print("Warning: --embed-combine 'single' provided with two genes. Defaulting to 'sum'.")
        return v1 + v2


def _cosine_similarities_to_query_vec(embeddings_df, qvec):
    """
    Compute cosine similarity of every row in embeddings_df to the query vector qvec.
    Returns a pandas Series indexed by gene.
    """
    # cosine_similarity handles normalization internally
    sims = cosine_similarity(embeddings_df.to_numpy(), qvec.reshape(1, -1)).ravel()
    return pd.Series(sims, index=embeddings_df.index)


def rank_genes_by_similarity(query_genes, emb_sim_df, coexpr_corr_df, embeddings_df,
                             embed_combine_mode="single", output_dir="gene_rankings", output_prefix=""):
    """
    Rank all genes by their similarity to a query (single gene or two-gene composite) and save results.
    
    Parameters
    ----------
    query_genes : list[str]
        One or two gene names
    emb_sim_df : pd.DataFrame
        Embedding cosine similarity matrix (gene x gene). Optional if using composite queries.
    coexpr_corr_df : pd.DataFrame
        Co-expression correlation matrix
    embeddings_df : pd.DataFrame
        Raw embeddings (gene x dims), used for composite queries or fallback
    embed_combine_mode : str
        'single' (one gene), 'sum', or 'mean_l2' (two genes)
    output_dir : str
        Directory to save output files
    output_prefix : str
        Prefix for output filenames
    
    Returns
    -------
    tuple : (embedding_rankings_df or None, coexpr_rankings_df or None)
    """
    # Create output directory if it doesn't exist
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    two_gene = len(query_genes) == 2
    qlabel = None  # label for filenames

    # ----- Embedding rankings -----
    emb_rankings = None
    if not embeddings_df.empty:
        if two_gene:
            # Build composite query vector and compute cosine to all genes
            mode = embed_combine_mode if embed_combine_mode in ("sum", "mean_l2") else "sum"
            qvec = _build_query_vector(query_genes, embeddings_df, mode)
            sims = _cosine_similarities_to_query_vec(embeddings_df, qvec)
            sims = sims.sort_values(ascending=False)

            qlabel = f"{query_genes[0]}+{query_genes[1]}_{mode}"
            emb_rankings = pd.DataFrame({
                'gene': sims.index,
                'embedding_cosine_similarity': sims.values,
                'rank': range(1, len(sims) + 1)
            })
            emb_output_file = output_path / f"{qlabel}_embedding_similarity_ranked.csv"
            emb_rankings.to_csv(emb_output_file, index=False)
            print(f"\nEmbedding similarity rankings (composite) saved to: {emb_output_file}")
            print(f"Top 10 genes by embedding similarity for composite [{qlabel}]:")
            print(emb_rankings.head(10).to_string(index=False))
        else:
            # Single gene query: prefer fast lookup via emb_sim_df if available
            qg = query_genes[0]
            if (not emb_sim_df.empty) and (qg in emb_sim_df.index):
                sims = emb_sim_df[qg].copy().sort_values(ascending=False)
            else:
                # Fallback: compute cosine to the single gene's embedding vector
                if qg not in embeddings_df.index:
                    print(f"Warning: {qg} not found in embeddings; skipping embedding rankings.")
                    sims = None
                else:
                    qvec = embeddings_df.loc[qg].to_numpy()
                    sims = _cosine_similarities_to_query_vec(embeddings_df, qvec).sort_values(ascending=False)

            if sims is not None:
                qlabel = qg
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
        # Use the first gene if two provided
        q_for_coexpr = query_genes[0]
        if two_gene:
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
        print("\nSkipping co-expression correlation rankings per --skip-coexpr flag")
    
    return emb_rankings, coexpr_rankings


def main():
    """Main function to run the gene similarity ranking."""
    
    # Set up argument parser
    parser = argparse.ArgumentParser(description='Rank genes by similarity to a query gene or two-gene composite (embedding mode supports two-gene queries).')
    # Backward-compat single-gene flag (deprecated)
    parser.add_argument('--gene', type=str, help='[Deprecated] Use --genes instead.')
    # New: one or two genes
    parser.add_argument('--genes', nargs='+', type=str, required=False,
                       help='One or two query gene names (e.g., --genes ifnphi1 or --genes ifnphi1 ifnphi2)')
    parser.add_argument('--embed-combine', type=str, choices=['single', 'sum', 'mean_l2'], default='single',
                       help="How to combine two gene embeddings for the query vector (use with two genes).")
    parser.add_argument('--mapping', type=str, 
                       default="/path/to/Geneformer_RQfork/data/RNAquarium/xref/drerio_z12-geneid_to_z11-ENSDARG.csv",
                       help='Path to gene ID mapping file')
    parser.add_argument('--embeddings', type=str,
                       default="/path/to/Geneformer_RQfork/embeddings/extracted_embeddings/RNAquarium_second_to_last_layer.csv",
                       help='Path to embeddings CSV file')
    parser.add_argument('--coexpr', type=str,
                       default="/path/to/Geneformer_RQfork/embeddings/coexpression/corr_matrix_corr_matrix_robust_winsorize_0.001/corr_matrix.npz",
                       help='Path to co-expression correlation matrix')
    parser.add_argument('--output-dir', type=str, default='gene_rankings',
                       help='Output directory for ranking files')
    parser.add_argument('--output-prefix', type=str, default='',
                       help='Output prefix for ranking files (applied to co-expression output only)')
    parser.add_argument('--skip-embeddings', action='store_true',
                       help='Skip loading embeddings and embedding similarity computation')
    parser.add_argument('--skip-coexpr', action='store_true',
                       help='Skip loading co-expression correlation matrix')

    # NEW: only compute & write full embeddings cosine similarity matrix
    parser.add_argument('--emb-sim-only', action='store_true',
                        help='Only compute and save the full gene x gene cosine similarity matrix from embeddings, then exit.')
    parser.add_argument('--emb-sim-output-file', type=str, default='embedding_cosine_similarity.csv',
                        help='Output filename for the full embedding cosine similarity matrix (used with --emb-sim-only).')

    args = parser.parse_args()

    # Backward compat: map --gene to --genes
    if args.genes is None and args.gene is not None:
        print("Notice: --gene is deprecated; please use --genes. Proceeding with single-gene query.")
        args.genes = [args.gene]
    if args.genes is None and not args.emb_sim_only:
        raise SystemExit("Please provide --genes <gene1> [gene2].")

    if args.genes is not None and len(args.genes) > 2:
        raise SystemExit("Please provide at most two genes for --genes.")

    # If two genes but user left combine mode at 'single', default to 'sum'
    if args.genes is not None and len(args.genes) == 2 and args.embed_combine == 'single':
        print("Two genes provided with --embed-combine 'single'; defaulting to 'sum'.")
        args.embed_combine = 'sum'
    
    print(f"Query genes: {args.genes}")
    print(f"Embedding combine mode: {args.embed_combine}")
    print(f"Output directory: {args.output_dir}")
    print(f"Output prefix: {args.output_prefix}")
    print("=" * 60)
    
    # Load and process data
    emb_sim_df, coexpr_corr_df, mapping_df, embeddings_df = load_and_process_data(
        args.mapping, 
        args.embeddings, 
        args.coexpr,
        skip_embeddings=args.skip_embeddings,
        skip_coexpr=args.skip_coexpr
    )

    # NEW: emb-sim-only mode -> write matrix and exit
    if args.emb_sim_only:
        if emb_sim_df.empty:
            raise SystemExit("Embedding similarity matrix is empty. Make sure embeddings are loaded (do not use --skip-embeddings).")
        # Resolve output path
        output_dir_path = Path(args.output_dir)
        output_dir_path.mkdir(parents=True, exist_ok=True)
        emb_sim_path = Path(args.emb_sim_output_file)
        if not emb_sim_path.is_absolute():
            emb_sim_path = output_dir_path / emb_sim_path
        print(f"\nWriting full embedding cosine similarity matrix to: {emb_sim_path}")
        emb_sim_df.to_csv(emb_sim_path)
        print("Done. Exiting per --emb-sim-only.")
        return
    
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
        output_prefix=args.output_prefix
    )
    
    # Print summary statistics
    print("\n" + "=" * 60)
    print("Summary Statistics")
    print("=" * 60)
    
    if emb_rankings is not None:
        if len(args.genes) == 1:
            q_label = args.genes[0]
            print(f"\nEmbedding similarity to {q_label}:")
            # Self-similarity exists only for single-gene queries
            print(f"  Self-similarity: {emb_rankings.iloc[0]['embedding_cosine_similarity']:.4f}")
        else:
            q_label = f"{args.genes[0]}+{args.genes[1]}_{args.embed_combine}"
            print(f"\nEmbedding similarity for composite query [{q_label}]:")
            top_gene = emb_rankings.iloc[0]['gene']
            top_sim = emb_rankings.iloc[0]['embedding_cosine_similarity']
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
