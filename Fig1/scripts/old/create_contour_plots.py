#!/usr/bin/env python3
"""
Create contour plots for RNAquarium Figure 1:
1. Mapped reads mate 1 vs mate 2 (contour per bioproject)
2. Mapped reads vs gene sparsity (separate contour by mate, per bioproject)

Uses seaborn kdeplot for contour visualization.
"""

import polars as pl
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from scipy.cluster import hierarchy
from scipy.spatial.distance import pdist, squareform
from scipy.stats import wasserstein_distance

# Data paths
SEQDETECTIVE_METRICS = Path("data/75k_unstable/seqdetective_metrics.parquet")
ROSETTA_PATH = Path("data/metadata/zf_rosetta.tsv")
OUTPUT_DIR = Path("figures")


def load_metrics():
    """Load seq-detective metrics with bioproject mapping."""
    if not SEQDETECTIVE_METRICS.exists():
        raise FileNotFoundError(f"Run extract_seqdetective_metrics.py first to create {SEQDETECTIVE_METRICS}")

    df = pl.read_parquet(SEQDETECTIVE_METRICS)
    print(f"Loaded {df.height} runs with seq-detective metrics")
    return df


def filter_top_bioprojects(df: pl.DataFrame, n_top: int = 10) -> pl.DataFrame:
    """Filter to top N bioprojects by run count."""
    if "bioproject" not in df.columns:
        print("Warning: no bioproject column found")
        return df

    top_projects = (
        df.filter(pl.col("bioproject").is_not_null())
        .group_by("bioproject")
        .agg(pl.len().alias("count"))
        .sort("count", descending=True)
        .head(n_top)["bioproject"]
        .to_list()
    )
    print(f"Top {n_top} bioprojects: {top_projects}")
    return df.filter(pl.col("bioproject").is_in(top_projects))


def plot_mate1_vs_mate2_mapping(df: pl.DataFrame, output_path: Path):
    """
    Create contour plot of mate1 vs mate2 mapping rates.
    One contour per bioproject for top bioprojects.
    """
    # Filter to paired-end samples only
    pe_df = df.filter(pl.col("n_mates") == 2)
    print(f"Paired-end samples: {pe_df.height}")

    if pe_df.height == 0:
        print("No paired-end samples found!")
        return

    # Get top bioprojects
    pe_top = filter_top_bioprojects(pe_df, n_top=10)

    fig, ax = plt.subplots(figsize=(10, 10))

    # Convert to pandas for seaborn
    plot_df = pe_top.select([
        "mapping_rate_m1",
        "mapping_rate_m2",
        "bioproject"
    ]).drop_nulls().to_pandas()

    if len(plot_df) == 0:
        print("No data after dropping nulls!")
        return

    # Create contour plot with hue by bioproject
    sns.kdeplot(
        data=plot_df,
        x="mapping_rate_m1",
        y="mapping_rate_m2",
        hue="bioproject",
        levels=5,
        alpha=0.7,
        ax=ax
    )

    # Add diagonal reference line
    ax.plot([0, 1], [0, 1], 'k--', alpha=0.3, label='y=x')

    ax.set_xlabel("Mate 1 Mapping Rate", fontsize=12)
    ax.set_ylabel("Mate 2 Mapping Rate", fontsize=12)
    ax.set_title("Paired-End Mapping Rates by BioProject\n(Seq-Detective QC Metrics)", fontsize=14)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.legend(title="BioProject", bbox_to_anchor=(1.05, 1), loc='upper left')

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"Saved {output_path}")

    # Also save PNG
    plt.savefig(output_path.with_suffix('.png'), dpi=300, bbox_inches='tight')

    plt.close()


def plot_mapping_vs_sparsity(df: pl.DataFrame, output_path: Path):
    """
    Create contour plot of mapping rate vs gene sparsity.
    Separate panels for mate 1 and mate 2, with contours per bioproject.
    """
    # Get top bioprojects
    df_top = filter_top_bioprojects(df, n_top=10)

    fig, axes = plt.subplots(1, 2, figsize=(16, 7))

    # Panel 1: Single-end and Mate 1
    # Combine SE samples (no _m1 suffix) with PE mate 1
    mate1_data = []

    # Single-end samples (columns without _m1/_m2 suffix)
    se_df = df_top.filter(pl.col("n_mates") == 1)
    if se_df.height > 0 and "mapping_rate" in se_df.columns:
        se_plot = se_df.select([
            pl.col("mapping_rate").alias("mapping_rate"),
            pl.col("sparsity").alias("sparsity"),
            "bioproject"
        ]).drop_nulls().with_columns(pl.lit("SE").alias("mate"))
        mate1_data.append(se_plot)

    # Paired-end mate 1
    pe_df = df_top.filter(pl.col("n_mates") == 2)
    if pe_df.height > 0 and "mapping_rate_m1" in pe_df.columns:
        m1_plot = pe_df.select([
            pl.col("mapping_rate_m1").alias("mapping_rate"),
            pl.col("sparsity_m1").alias("sparsity"),
            "bioproject"
        ]).drop_nulls().with_columns(pl.lit("M1").alias("mate"))
        mate1_data.append(m1_plot)

    if mate1_data:
        combined_m1 = pl.concat(mate1_data).to_pandas()
        if len(combined_m1) > 0:
            sns.kdeplot(
                data=combined_m1,
                x="mapping_rate",
                y="sparsity",
                hue="bioproject",
                levels=5,
                alpha=0.7,
                ax=axes[0]
            )
    axes[0].set_xlabel("Mapping Rate", fontsize=12)
    axes[0].set_ylabel("Gene Sparsity", fontsize=12)
    axes[0].set_title("Single-End / Mate 1", fontsize=14)
    axes[0].set_xlim(0, 1)
    axes[0].set_ylim(0, 1)

    # Panel 2: Mate 2 only
    if pe_df.height > 0 and "mapping_rate_m2" in pe_df.columns:
        m2_plot = pe_df.select([
            pl.col("mapping_rate_m2").alias("mapping_rate"),
            pl.col("sparsity_m2").alias("sparsity"),
            "bioproject"
        ]).drop_nulls().to_pandas()

        if len(m2_plot) > 0:
            sns.kdeplot(
                data=m2_plot,
                x="mapping_rate",
                y="sparsity",
                hue="bioproject",
                levels=5,
                alpha=0.7,
                ax=axes[1]
            )
    axes[1].set_xlabel("Mapping Rate", fontsize=12)
    axes[1].set_ylabel("Gene Sparsity", fontsize=12)
    axes[1].set_title("Mate 2", fontsize=14)
    axes[1].set_xlim(0, 1)
    axes[1].set_ylim(0, 1)

    # Shared legend
    handles, labels = axes[0].get_legend_handles_labels()
    if handles:
        axes[0].get_legend().remove()
        axes[1].get_legend().remove() if axes[1].get_legend() else None
        fig.legend(handles, labels, title="BioProject",
                   bbox_to_anchor=(1.02, 0.5), loc='center left')

    fig.suptitle("Mapping Rate vs Gene Sparsity by BioProject\n(Seq-Detective QC Metrics)",
                 fontsize=16, y=1.02)

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"Saved {output_path}")

    plt.savefig(output_path.with_suffix('.png'), dpi=300, bbox_inches='tight')
    plt.close()


def compute_2d_wasserstein_distance(samples1: np.ndarray, samples2: np.ndarray, n_bins: int = 50) -> float:
    """
    Compute approximate 2D Wasserstein distance between two point clouds using histograms.

    Args:
        samples1: Array of shape (n1, 2) with (mate1, mate2) values
        samples2: Array of shape (n2, 2) with (mate1, mate2) values
        n_bins: Number of bins for histogram representation

    Returns:
        Wasserstein distance between the two distributions
    """
    # Create 2D histograms
    bins = np.linspace(0, 1, n_bins + 1)

    hist1, xedges, yedges = np.histogram2d(
        samples1[:, 0], samples1[:, 1],
        bins=[bins, bins],
        density=True
    )
    hist2, _, _ = np.histogram2d(
        samples2[:, 0], samples2[:, 1],
        bins=[bins, bins],
        density=True
    )

    # Normalize to probability distributions
    hist1 = hist1.flatten()
    hist2 = hist2.flatten()
    hist1 = hist1 / (hist1.sum() + 1e-10)
    hist2 = hist2 / (hist2.sum() + 1e-10)

    # Compute bin centers for distance computation
    x_centers = (xedges[:-1] + xedges[1:]) / 2
    y_centers = (yedges[:-1] + yedges[1:]) / 2
    xx, yy = np.meshgrid(x_centers, y_centers)
    positions = np.column_stack([xx.ravel(), yy.ravel()])

    # Compute pairwise Euclidean distances between bins
    dist_matrix = np.sqrt(
        np.sum((positions[:, None, :] - positions[None, :, :]) ** 2, axis=2)
    )

    # Use linear programming to solve optimal transport
    # Simple approximation: weighted sum of distances
    # For exact solution, we'd use ot.emd2(), but this approximation is faster
    distance = 0.0
    for i in range(len(hist1)):
        for j in range(len(hist2)):
            distance += hist1[i] * hist2[j] * dist_matrix[i, j]

    return distance


def compute_sliced_wasserstein_distance(samples1: np.ndarray, samples2: np.ndarray, n_projections: int = 50) -> float:
    """
    Compute sliced Wasserstein distance between two 2D point clouds.
    This is much faster than exact Wasserstein distance and works well for 2D distributions.

    The sliced Wasserstein distance projects the 2D distributions onto multiple 1D lines
    and averages the 1D Wasserstein distances.

    Args:
        samples1: Array of shape (n1, 2) with (mate1, mate2) values
        samples2: Array of shape (n2, 2) with (mate1, mate2) values
        n_projections: Number of random projection directions

    Returns:
        Sliced Wasserstein distance between the two distributions
    """
    distances = []

    for _ in range(n_projections):
        # Random projection direction
        theta = np.random.uniform(0, 2 * np.pi)
        direction = np.array([np.cos(theta), np.sin(theta)])

        # Project samples onto this direction
        proj1 = samples1 @ direction
        proj2 = samples2 @ direction

        # Compute 1D Wasserstein distance
        dist = wasserstein_distance(proj1, proj2)
        distances.append(dist)

    return np.mean(distances)


def extract_bioproject_distributions(
    df: pl.DataFrame,
    min_samples: int = 5,
    max_samples_per_bioproject: int = 1000
) -> tuple[dict[str, np.ndarray], list[str]]:
    """
    Extract mate1-mate2 sample distributions for each bioproject.

    Args:
        df: DataFrame with seq-detective metrics
        min_samples: Minimum samples per bioproject
        max_samples_per_bioproject: Downsample to this many samples for efficiency

    Returns:
        distributions: Dict mapping bioproject -> array of shape (n_samples, 2)
        bioproject_names: Sorted list of bioproject IDs
    """
    # Filter to paired-end samples
    pe_df = df.filter(pl.col("n_mates") == 2)

    # Filter to bioprojects with enough samples
    bioproject_counts = (
        pe_df.filter(pl.col("bioproject").is_not_null())
        .group_by("bioproject")
        .agg(pl.len().alias("n_samples"))
        .filter(pl.col("n_samples") >= min_samples)
    )

    valid_bioprojects = bioproject_counts["bioproject"].to_list()
    pe_filtered = pe_df.filter(
        pl.col("bioproject").is_in(valid_bioprojects) &
        pl.col("mapping_rate_m1").is_not_null() &
        pl.col("mapping_rate_m2").is_not_null()
    )

    print(f"\nExtracting distributions for {len(valid_bioprojects)} bioprojects")
    print(f"Total samples: {pe_filtered.height}")

    # Extract distributions per bioproject
    distributions = {}
    for bioproject in sorted(valid_bioprojects):
        bp_data = pe_filtered.filter(pl.col("bioproject") == bioproject)
        samples = bp_data.select([
            "mapping_rate_m1",
            "mapping_rate_m2"
        ]).to_numpy()

        # Downsample if too many samples (for computational efficiency)
        if len(samples) > max_samples_per_bioproject:
            indices = np.random.choice(len(samples), max_samples_per_bioproject, replace=False)
            samples = samples[indices]

        distributions[bioproject] = samples

    bioproject_names = sorted(distributions.keys())
    return distributions, bioproject_names


def compute_wasserstein_distance_matrix(
    distributions: dict[str, np.ndarray],
    bioproject_names: list[str],
    method: str = "sliced",
    n_projections: int = 50
) -> np.ndarray:
    """
    Compute pairwise Wasserstein distances between bioproject distributions.

    Args:
        distributions: Dict mapping bioproject -> samples array
        bioproject_names: List of bioproject IDs
        method: "sliced" (fast) or "histogram" (slower but more accurate)
        n_projections: Number of projections for sliced Wasserstein

    Returns:
        distance_matrix: Symmetric matrix of pairwise distances
    """
    n = len(bioproject_names)
    distance_matrix = np.zeros((n, n))

    print(f"\nComputing pairwise Wasserstein distances ({method} method)...")
    print(f"Total comparisons: {n * (n - 1) // 2}")

    for i in range(n):
        if i % 100 == 0:
            print(f"  Progress: {i}/{n} bioprojects")

        for j in range(i + 1, n):
            bp1, bp2 = bioproject_names[i], bioproject_names[j]
            samples1 = distributions[bp1]
            samples2 = distributions[bp2]

            if method == "sliced":
                dist = compute_sliced_wasserstein_distance(samples1, samples2, n_projections)
            elif method == "histogram":
                dist = compute_2d_wasserstein_distance(samples1, samples2)
            else:
                raise ValueError(f"Unknown method: {method}")

            distance_matrix[i, j] = dist
            distance_matrix[j, i] = dist

    print(f"  Done! Distance matrix shape: {distance_matrix.shape}")
    return distance_matrix


def extract_bioproject_features(df: pl.DataFrame, min_samples: int = 5) -> tuple[np.ndarray, list[str], pl.DataFrame]:
    """
    Extract features for each bioproject's mate1-mate2 distribution.

    Args:
        df: DataFrame with seq-detective metrics
        min_samples: Minimum number of samples per bioproject to include

    Returns:
        features: 2D array of shape (n_bioprojects, n_features)
        bioproject_names: list of bioproject IDs
        features_df: DataFrame with features and sample counts
    """
    # Filter to paired-end samples only
    pe_df = df.filter(pl.col("n_mates") == 2)

    # Filter to bioprojects with enough samples
    bioproject_counts = (
        pe_df.filter(pl.col("bioproject").is_not_null())
        .group_by("bioproject")
        .agg(pl.len().alias("n_samples"))
        .filter(pl.col("n_samples") >= min_samples)
    )

    valid_bioprojects = bioproject_counts["bioproject"].to_list()
    pe_filtered = pe_df.filter(pl.col("bioproject").is_in(valid_bioprojects))

    print(f"\nFiltering to bioprojects with >={min_samples} samples:")
    print(f"  Total bioprojects: {len(valid_bioprojects)}")
    print(f"  Total samples: {pe_filtered.height}")

    # Compute features per bioproject
    features_df = (
        pe_filtered
        .filter(
            pl.col("mapping_rate_m1").is_not_null() &
            pl.col("mapping_rate_m2").is_not_null()
        )
        .group_by("bioproject")
        .agg([
            pl.len().alias("n_samples"),

            # Mate 1 statistics
            pl.col("mapping_rate_m1").mean().alias("m1_mean"),
            pl.col("mapping_rate_m1").std().alias("m1_std"),
            pl.col("mapping_rate_m1").median().alias("m1_median"),
            pl.col("mapping_rate_m1").quantile(0.25).alias("m1_q25"),
            pl.col("mapping_rate_m1").quantile(0.75).alias("m1_q75"),

            # Mate 2 statistics
            pl.col("mapping_rate_m2").mean().alias("m2_mean"),
            pl.col("mapping_rate_m2").std().alias("m2_std"),
            pl.col("mapping_rate_m2").median().alias("m2_median"),
            pl.col("mapping_rate_m2").quantile(0.25).alias("m2_q25"),
            pl.col("mapping_rate_m2").quantile(0.75).alias("m2_q75"),

            # Cross-mate features
            pl.corr("mapping_rate_m1", "mapping_rate_m2").alias("correlation"),
            (pl.col("mapping_rate_m1") - pl.col("mapping_rate_m2")).abs().mean().alias("mean_abs_diff"),
        ])
        .sort("bioproject")
    )

    bioproject_names = features_df["bioproject"].to_list()
    feature_matrix = features_df.select(pl.exclude("bioproject")).to_numpy()

    # Handle NaN and infinite values
    # Replace NaN std with 0 (no variance)
    # Replace other NaN with column median
    # Replace inf with column max
    print(f"\nChecking for NaN/inf values...")
    nan_mask = ~np.isfinite(feature_matrix)
    if nan_mask.any():
        print(f"  Found {nan_mask.sum()} non-finite values, cleaning...")
        for col_idx in range(feature_matrix.shape[1]):
            col = feature_matrix[:, col_idx]
            finite_values = col[np.isfinite(col)]
            if len(finite_values) > 0:
                # Replace inf with max finite value
                col[np.isinf(col)] = finite_values.max()
                # Replace NaN with median (or 0 for std columns)
                if col_idx in [2, 7]:  # m1_std, m2_std
                    col[np.isnan(col)] = 0
                else:
                    col[np.isnan(col)] = np.median(finite_values)
            else:
                # If entire column is non-finite, set to 0
                col[:] = 0

    return feature_matrix, bioproject_names, features_df


def assign_clusters(linkage_matrix: np.ndarray, n_clusters: int) -> np.ndarray:
    """
    Cut dendrogram to get cluster assignments.

    Args:
        linkage_matrix: Linkage matrix from hierarchical clustering
        n_clusters: Number of clusters to create

    Returns:
        cluster_labels: Array of cluster IDs (0 to n_clusters-1)
    """
    return hierarchy.fcluster(linkage_matrix, n_clusters, criterion='maxclust') - 1


def print_cluster_summary(
    cluster_labels: np.ndarray,
    bioproject_names: list[str],
    features_df: pl.DataFrame
):
    """
    Print summary statistics for each cluster.

    Args:
        cluster_labels: Array of cluster assignments
        bioproject_names: List of bioproject IDs
        features_df: DataFrame with features and metadata
    """
    # Create cluster assignment dataframe
    cluster_df = pl.DataFrame({
        "bioproject": bioproject_names,
        "cluster": cluster_labels
    })

    # Join with features
    cluster_features = cluster_df.join(features_df, on="bioproject")

    print("\n" + "="*80)
    print("CLUSTER SUMMARY")
    print("="*80)

    for cluster_id in sorted(cluster_df["cluster"].unique()):
        cluster_projects = cluster_features.filter(pl.col("cluster") == cluster_id)
        n_projects = cluster_projects.height
        total_samples = cluster_projects["n_samples"].sum()

        print(f"\n--- Cluster {cluster_id} ({n_projects} bioprojects, {total_samples} samples) ---")
        print(f"BioProjects: {', '.join(cluster_projects['bioproject'].to_list())}")
        print(f"\nMapping Rate Statistics:")

        # Handle std() returning None for single-bioproject clusters
        m1_std = cluster_projects['m1_mean'].std()
        m2_std = cluster_projects['m2_mean'].std()
        corr_std = cluster_projects['correlation'].std()
        diff_std = cluster_projects['mean_abs_diff'].std()

        print(f"  Mate1 mean: {cluster_projects['m1_mean'].mean():.3f} ± {m1_std if m1_std is not None else 0:.3f}")
        print(f"  Mate2 mean: {cluster_projects['m2_mean'].mean():.3f} ± {m2_std if m2_std is not None else 0:.3f}")
        print(f"  Correlation: {cluster_projects['correlation'].mean():.3f} ± {corr_std if corr_std is not None else 0:.3f}")
        print(f"  Mean abs diff: {cluster_projects['mean_abs_diff'].mean():.3f} ± {diff_std if diff_std is not None else 0:.3f}")

    print("\n" + "="*80)

    # Save detailed cluster assignments
    output_file = OUTPUT_DIR / "cluster_assignments.tsv"
    cluster_features.write_csv(output_file, separator="\t")
    print(f"\nDetailed cluster assignments saved to: {output_file}")


def cluster_with_wasserstein_and_plot(
    df: pl.DataFrame,
    output_path: Path,
    method: str = "average",
    n_clusters: int = 8,
    min_samples: int = 5,
    wasserstein_method: str = "sliced",
    n_projections: int = 50,
    max_samples_per_bioproject: int = 1000
):
    """
    Perform hierarchical clustering using Wasserstein distances between
    mate1-mate2 distributions and plot dendrogram.

    Args:
        df: DataFrame with seq-detective metrics
        output_path: Where to save dendrogram
        method: Linkage method ('average', 'complete', 'single')
                Note: 'ward' requires Euclidean distances, use 'average' for Wasserstein
        n_clusters: Number of clusters to identify
        min_samples: Minimum samples per bioproject to include
        wasserstein_method: "sliced" (fast) or "histogram" (slower)
        n_projections: Number of projections for sliced Wasserstein
        max_samples_per_bioproject: Downsample to this many samples for speed
    """
    # Extract distributions
    distributions, bioproject_names = extract_bioproject_distributions(
        df, min_samples, max_samples_per_bioproject
    )

    if len(bioproject_names) < 2:
        print("Need at least 2 bioprojects for clustering!")
        return None

    # Also extract feature dataframe for summary statistics
    _, _, features_df = extract_bioproject_features(df, min_samples)

    print(f"\nClustering {len(bioproject_names)} bioprojects using Wasserstein distances")

    # Compute Wasserstein distance matrix
    distance_matrix = compute_wasserstein_distance_matrix(
        distributions, bioproject_names, wasserstein_method, n_projections
    )

    # Convert to condensed distance matrix for scipy
    from scipy.spatial.distance import squareform
    condensed_distances = squareform(distance_matrix)

    # Perform hierarchical clustering
    print(f"\nPerforming hierarchical clustering (method={method})...")
    linkage_matrix = hierarchy.linkage(condensed_distances, method=method)

    # Assign clusters
    cluster_labels = assign_clusters(linkage_matrix, n_clusters)

    # Plot dendrogram
    fig, ax = plt.subplots(figsize=(16, 10))

    # Color threshold to match n_clusters
    max_d = linkage_matrix[-n_clusters+1, 2] if n_clusters > 1 else linkage_matrix[-1, 2]
    dendro = hierarchy.dendrogram(
        linkage_matrix,
        labels=bioproject_names,
        ax=ax,
        leaf_rotation=90,
        leaf_font_size=8,
        color_threshold=max_d
    )

    ax.axhline(y=max_d, color='r', linestyle='--', linewidth=1, alpha=0.5)
    ax.text(0.02, max_d, f'{n_clusters} clusters', va='bottom', fontsize=10, color='r')

    ax.set_xlabel("BioProject", fontsize=12)
    ax.set_ylabel("Wasserstein Distance", fontsize=12)
    ax.set_title(
        f"Hierarchical Clustering of Mate1-Mate2 Distributions (Wasserstein Distance)\n"
        f"(Method: {method}, Clusters: {n_clusters}, Wasserstein: {wasserstein_method})",
        fontsize=14
    )

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"\nSaved dendrogram to {output_path}")

    plt.savefig(output_path.with_suffix('.png'), dpi=300, bbox_inches='tight')
    plt.close()

    # Print cluster summary
    print_cluster_summary(cluster_labels, bioproject_names, features_df)

    # Return results for further analysis
    return {
        'linkage_matrix': linkage_matrix,
        'bioproject_names': bioproject_names,
        'distance_matrix': distance_matrix,
        'features_df': features_df,
        'cluster_labels': cluster_labels,
        'distributions': distributions
    }


def cluster_and_plot_dendrogram(
    df: pl.DataFrame,
    output_path: Path,
    method: str = "ward",
    metric: str = "euclidean",
    n_clusters: int = 5,
    min_samples: int = 5
):
    """
    Perform hierarchical clustering on bioproject mate1-mate2 distributions
    and plot dendrogram.

    Args:
        df: DataFrame with seq-detective metrics
        output_path: Where to save dendrogram
        method: Linkage method ('ward', 'average', 'complete', 'single')
        metric: Distance metric ('euclidean', 'cosine', 'correlation')
        n_clusters: Number of clusters to identify
        min_samples: Minimum samples per bioproject to include
    """
    # Extract features
    features, bioproject_names, features_df = extract_bioproject_features(df, min_samples=min_samples)

    if len(bioproject_names) < 2:
        print("Need at least 2 bioprojects for clustering!")
        return None

    print(f"\nExtracted features for {len(bioproject_names)} bioprojects")
    print(f"Feature matrix shape: {features.shape}")

    # Perform hierarchical clustering
    print(f"\nPerforming hierarchical clustering (method={method}, metric={metric})...")
    linkage_matrix = hierarchy.linkage(features, method=method, metric=metric)

    # Assign clusters
    cluster_labels = assign_clusters(linkage_matrix, n_clusters)

    # Plot dendrogram
    fig, ax = plt.subplots(figsize=(16, 10))

    # Color threshold to match n_clusters
    max_d = linkage_matrix[-n_clusters+1, 2]
    dendro = hierarchy.dendrogram(
        linkage_matrix,
        labels=bioproject_names,
        ax=ax,
        leaf_rotation=90,
        leaf_font_size=8,
        color_threshold=max_d
    )

    ax.axhline(y=max_d, color='r', linestyle='--', linewidth=1, alpha=0.5)
    ax.text(0.02, max_d, f'{n_clusters} clusters', va='bottom', fontsize=10, color='r')

    ax.set_xlabel("BioProject", fontsize=12)
    ax.set_ylabel("Distance", fontsize=12)
    ax.set_title(
        f"Hierarchical Clustering of Mate1-Mate2 Mapping Rate Distributions\n"
        f"(Method: {method}, Metric: {metric}, Clusters: {n_clusters})",
        fontsize=14
    )

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"\nSaved dendrogram to {output_path}")

    plt.savefig(output_path.with_suffix('.png'), dpi=300, bbox_inches='tight')
    plt.close()

    # Print cluster summary
    print_cluster_summary(cluster_labels, bioproject_names, features_df)

    # Return results for further analysis
    return {
        'linkage_matrix': linkage_matrix,
        'bioproject_names': bioproject_names,
        'features': features,
        'features_df': features_df,
        'cluster_labels': cluster_labels
    }


def plot_cluster_hexbins(
    df: pl.DataFrame,
    cluster_labels: np.ndarray,
    bioproject_names: list[str],
    output_path: Path,
    gridsize: int = 30
):
    """
    Create hexbin heatmaps for each cluster showing mate1 vs mate2 distributions.

    Args:
        df: DataFrame with seq-detective metrics
        cluster_labels: Array of cluster assignments
        bioproject_names: List of bioproject IDs
        output_path: Where to save the plot
        gridsize: Number of hexagons in x direction
    """
    # Create cluster assignment mapping
    cluster_map = {bp: cluster for bp, cluster in zip(bioproject_names, cluster_labels)}

    # Get unique clusters, sorted
    unique_clusters = sorted(np.unique(cluster_labels))
    n_clusters = len(unique_clusters)

    # Calculate grid layout (prefer wider layouts)
    ncols = min(3, n_clusters)
    nrows = int(np.ceil(n_clusters / ncols))

    fig, axes = plt.subplots(nrows, ncols, figsize=(6*ncols, 5*nrows))
    if n_clusters == 1:
        axes = np.array([axes])
    axes = axes.flatten()

    # Filter to paired-end samples
    pe_df = df.filter(pl.col("n_mates") == 2).filter(
        pl.col("mapping_rate_m1").is_not_null() &
        pl.col("mapping_rate_m2").is_not_null() &
        pl.col("bioproject").is_in(bioproject_names)
    )

    # Add cluster assignment to dataframe
    cluster_assignment_df = pl.DataFrame({
        "bioproject": bioproject_names,
        "cluster": cluster_labels
    })

    pe_with_clusters = pe_df.join(cluster_assignment_df, on="bioproject")

    # Plot each cluster
    for idx, cluster_id in enumerate(unique_clusters):
        ax = axes[idx]

        # Get data for this cluster
        cluster_data = pe_with_clusters.filter(pl.col("cluster") == cluster_id)

        if cluster_data.height == 0:
            ax.text(0.5, 0.5, "No data", ha='center', va='center', transform=ax.transAxes)
            ax.set_title(f"Cluster {cluster_id}")
            continue

        # Extract mate1 and mate2 values
        m1 = cluster_data["mapping_rate_m1"].to_numpy()
        m2 = cluster_data["mapping_rate_m2"].to_numpy()

        # Create hexbin
        hb = ax.hexbin(
            m1, m2,
            gridsize=gridsize,
            cmap='viridis',
            mincnt=1,
            extent=(0, 1, 0, 1),
            linewidths=0.1,
            edgecolors='face'
        )

        # Add colorbar
        cbar = plt.colorbar(hb, ax=ax)
        cbar.set_label('Sample Count', fontsize=10)

        # Add diagonal reference line
        ax.plot([0, 1], [0, 1], 'r--', alpha=0.5, linewidth=1, label='y=x')

        # Get cluster statistics
        n_bioprojects = len(cluster_data["bioproject"].unique())
        n_samples = cluster_data.height
        mean_m1 = m1.mean()
        mean_m2 = m2.mean()

        ax.set_xlabel("Mate 1 Mapping Rate", fontsize=11)
        ax.set_ylabel("Mate 2 Mapping Rate", fontsize=11)
        ax.set_title(
            f"Cluster {cluster_id}\n"
            f"{n_bioprojects} bioprojects, {n_samples:,} samples\n"
            f"Mean: M1={mean_m1:.2f}, M2={mean_m2:.2f}",
            fontsize=12
        )
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.set_aspect('equal')
        ax.legend(loc='upper left', fontsize=9)

    # Hide unused subplots
    for idx in range(n_clusters, len(axes)):
        axes[idx].axis('off')

    fig.suptitle(
        "Mate1-Mate2 Mapping Rate Distributions by Cluster\n(Hexbin Heatmaps)",
        fontsize=16,
        y=0.995
    )

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"\nSaved cluster hexbin plots to {output_path}")

    plt.savefig(output_path.with_suffix('.png'), dpi=300, bbox_inches='tight')
    plt.close()


def extract_sparsity_distributions(
    df: pl.DataFrame,
    min_samples: int = 5,
    max_samples_per_bioproject: int = 1000
) -> tuple[dict[str, np.ndarray], list[str]]:
    """
    Extract mapping rate vs sparsity distributions for each bioproject (using mate 1).

    Args:
        df: DataFrame with seq-detective metrics
        min_samples: Minimum samples per bioproject
        max_samples_per_bioproject: Downsample to this many samples for efficiency

    Returns:
        distributions: Dict mapping bioproject -> array of shape (n_samples, 2)
                       where columns are [mapping_rate_m1, sparsity_m1]
        bioproject_names: Sorted list of bioproject IDs
    """
    # Filter to paired-end samples with mate 1 data
    pe_df = df.filter(pl.col("n_mates") == 2)

    # Filter to bioprojects with enough samples
    bioproject_counts = (
        pe_df.filter(pl.col("bioproject").is_not_null())
        .group_by("bioproject")
        .agg(pl.len().alias("n_samples"))
        .filter(pl.col("n_samples") >= min_samples)
    )

    valid_bioprojects = bioproject_counts["bioproject"].to_list()
    pe_filtered = pe_df.filter(
        pl.col("bioproject").is_in(valid_bioprojects) &
        pl.col("mapping_rate_m1").is_not_null() &
        pl.col("sparsity_m1").is_not_null()
    )

    print(f"\nExtracting mapping vs sparsity distributions for {len(valid_bioprojects)} bioprojects")
    print(f"Total samples: {pe_filtered.height}")

    # Extract distributions per bioproject
    distributions = {}
    for bioproject in sorted(valid_bioprojects):
        bp_data = pe_filtered.filter(pl.col("bioproject") == bioproject)
        samples = bp_data.select([
            "mapping_rate_m1",
            "sparsity_m1"
        ]).to_numpy()

        # Downsample if too many samples
        if len(samples) > max_samples_per_bioproject:
            indices = np.random.choice(len(samples), max_samples_per_bioproject, replace=False)
            samples = samples[indices]

        distributions[bioproject] = samples

    bioproject_names = sorted(distributions.keys())
    return distributions, bioproject_names


def cluster_sparsity_with_wasserstein(
    df: pl.DataFrame,
    output_path: Path,
    method: str = "average",
    n_clusters: int = 8,
    min_samples: int = 50,
    wasserstein_method: str = "sliced",
    n_projections: int = 50,
    max_samples_per_bioproject: int = 500
):
    """
    Perform hierarchical clustering using Wasserstein distances between
    mapping rate vs sparsity distributions.

    Args:
        df: DataFrame with seq-detective metrics
        output_path: Where to save dendrogram
        method: Linkage method
        n_clusters: Number of clusters to identify
        min_samples: Minimum samples per bioproject to include
        wasserstein_method: "sliced" (fast) or "histogram" (slower)
        n_projections: Number of projections for sliced Wasserstein
        max_samples_per_bioproject: Downsample to this many samples for speed
    """
    # Extract sparsity distributions
    distributions, bioproject_names = extract_sparsity_distributions(
        df, min_samples, max_samples_per_bioproject
    )

    if len(bioproject_names) < 2:
        print("Need at least 2 bioprojects for clustering!")
        return None

    print(f"\nClustering {len(bioproject_names)} bioprojects using Wasserstein distances")

    # Compute Wasserstein distance matrix
    distance_matrix = compute_wasserstein_distance_matrix(
        distributions, bioproject_names, wasserstein_method, n_projections
    )

    # Convert to condensed distance matrix for scipy
    condensed_distances = squareform(distance_matrix)

    # Perform hierarchical clustering
    print(f"\nPerforming hierarchical clustering (method={method})...")
    linkage_matrix = hierarchy.linkage(condensed_distances, method=method)

    # Assign clusters
    cluster_labels = assign_clusters(linkage_matrix, n_clusters)

    # Plot dendrogram
    fig, ax = plt.subplots(figsize=(16, 10))

    # Color threshold to match n_clusters
    max_d = linkage_matrix[-n_clusters+1, 2] if n_clusters > 1 else linkage_matrix[-1, 2]
    dendro = hierarchy.dendrogram(
        linkage_matrix,
        labels=bioproject_names,
        ax=ax,
        leaf_rotation=90,
        leaf_font_size=8,
        color_threshold=max_d
    )

    ax.axhline(y=max_d, color='r', linestyle='--', linewidth=1, alpha=0.5)
    ax.text(0.02, max_d, f'{n_clusters} clusters', va='bottom', fontsize=10, color='r')

    ax.set_xlabel("BioProject", fontsize=12)
    ax.set_ylabel("Wasserstein Distance", fontsize=12)
    ax.set_title(
        f"Hierarchical Clustering of Mapping Rate vs Sparsity Distributions\n"
        f"(Method: {method}, Clusters: {n_clusters}, Wasserstein: {wasserstein_method})",
        fontsize=14
    )

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"\nSaved sparsity dendrogram to {output_path}")

    plt.savefig(output_path.with_suffix('.png'), dpi=300, bbox_inches='tight')
    plt.close()

    # Compute summary statistics for clusters
    print("\n" + "="*80)
    print("SPARSITY CLUSTER SUMMARY")
    print("="*80)

    # Extract features for each bioproject
    pe_df = df.filter(pl.col("n_mates") == 2).filter(
        pl.col("bioproject").is_in(bioproject_names) &
        pl.col("mapping_rate_m1").is_not_null() &
        pl.col("sparsity_m1").is_not_null()
    )

    features_df = (
        pe_df
        .group_by("bioproject")
        .agg([
            pl.len().alias("n_samples"),
            pl.col("mapping_rate_m1").mean().alias("mapping_mean"),
            pl.col("mapping_rate_m1").std().alias("mapping_std"),
            pl.col("sparsity_m1").mean().alias("sparsity_mean"),
            pl.col("sparsity_m1").std().alias("sparsity_std"),
            pl.corr("mapping_rate_m1", "sparsity_m1").alias("correlation"),
        ])
        .sort("bioproject")
    )

    # Create cluster assignment dataframe
    cluster_df = pl.DataFrame({
        "bioproject": bioproject_names,
        "cluster": cluster_labels
    })

    cluster_features = cluster_df.join(features_df, on="bioproject")

    for cluster_id in sorted(cluster_df["cluster"].unique()):
        cluster_projects = cluster_features.filter(pl.col("cluster") == cluster_id)
        n_projects = cluster_projects.height
        total_samples = cluster_projects["n_samples"].sum()

        print(f"\n--- Cluster {cluster_id} ({n_projects} bioprojects, {total_samples} samples) ---")
        print(f"BioProjects: {', '.join(cluster_projects['bioproject'].to_list()[:10])}" +
              (f"... (+{n_projects-10} more)" if n_projects > 10 else ""))
        print(f"\nStatistics:")

        # Handle std() returning None for single-bioproject clusters
        mapping_std = cluster_projects['mapping_mean'].std()
        sparsity_std = cluster_projects['sparsity_mean'].std()
        corr_std = cluster_projects['correlation'].std()

        print(f"  Mapping rate: {cluster_projects['mapping_mean'].mean():.3f} ± {mapping_std if mapping_std is not None else 0:.3f}")
        print(f"  Sparsity: {cluster_projects['sparsity_mean'].mean():.3f} ± {sparsity_std if sparsity_std is not None else 0:.3f}")
        print(f"  Correlation: {cluster_projects['correlation'].mean():.3f} ± {corr_std if corr_std is not None else 0:.3f}")

    print("\n" + "="*80)

    # Save cluster assignments
    output_file = OUTPUT_DIR / "cluster_assignments_sparsity.tsv"
    cluster_features.write_csv(output_file, separator="\t")
    print(f"\nDetailed cluster assignments saved to: {output_file}")

    return {
        'linkage_matrix': linkage_matrix,
        'bioproject_names': bioproject_names,
        'distance_matrix': distance_matrix,
        'cluster_labels': cluster_labels,
        'distributions': distributions,
        'features_df': cluster_features
    }


def plot_sparsity_cluster_hexbins(
    df: pl.DataFrame,
    cluster_labels: np.ndarray,
    bioproject_names: list[str],
    output_path: Path,
    gridsize: int = 30
):
    """
    Create hexbin heatmaps for each cluster showing mapping rate vs sparsity.

    Args:
        df: DataFrame with seq-detective metrics
        cluster_labels: Array of cluster assignments
        bioproject_names: List of bioproject IDs
        output_path: Where to save the plot
        gridsize: Number of hexagons in x direction
    """
    # Get unique clusters, sorted
    unique_clusters = sorted(np.unique(cluster_labels))
    n_clusters = len(unique_clusters)

    # Calculate grid layout
    ncols = min(3, n_clusters)
    nrows = int(np.ceil(n_clusters / ncols))

    fig, axes = plt.subplots(nrows, ncols, figsize=(6*ncols, 5*nrows))
    if n_clusters == 1:
        axes = np.array([axes])
    axes = axes.flatten()

    # Filter to paired-end samples
    pe_df = df.filter(pl.col("n_mates") == 2).filter(
        pl.col("mapping_rate_m1").is_not_null() &
        pl.col("sparsity_m1").is_not_null() &
        pl.col("bioproject").is_in(bioproject_names)
    )

    # Add cluster assignment to dataframe
    cluster_assignment_df = pl.DataFrame({
        "bioproject": bioproject_names,
        "cluster": cluster_labels
    })

    pe_with_clusters = pe_df.join(cluster_assignment_df, on="bioproject")

    # Plot each cluster
    for idx, cluster_id in enumerate(unique_clusters):
        ax = axes[idx]

        # Get data for this cluster
        cluster_data = pe_with_clusters.filter(pl.col("cluster") == cluster_id)

        if cluster_data.height == 0:
            ax.text(0.5, 0.5, "No data", ha='center', va='center', transform=ax.transAxes)
            ax.set_title(f"Cluster {cluster_id}")
            continue

        # Extract mapping rate and sparsity values
        mapping = cluster_data["mapping_rate_m1"].to_numpy()
        sparsity = cluster_data["sparsity_m1"].to_numpy()

        # Create hexbin
        hb = ax.hexbin(
            mapping, sparsity,
            gridsize=gridsize,
            cmap='viridis',
            mincnt=1,
            extent=(0, 1, 0, 1),
            linewidths=0.1,
            edgecolors='face'
        )

        # Add colorbar
        cbar = plt.colorbar(hb, ax=ax)
        cbar.set_label('Sample Count', fontsize=10)

        # Get cluster statistics
        n_bioprojects = len(cluster_data["bioproject"].unique())
        n_samples = cluster_data.height
        mean_mapping = mapping.mean()
        mean_sparsity = sparsity.mean()

        ax.set_xlabel("Mapping Rate (Mate 1)", fontsize=11)
        ax.set_ylabel("Gene Sparsity (Mate 1)", fontsize=11)
        ax.set_title(
            f"Cluster {cluster_id}\n"
            f"{n_bioprojects} bioprojects, {n_samples:,} samples\n"
            f"Mean: Mapping={mean_mapping:.2f}, Sparsity={mean_sparsity:.2f}",
            fontsize=12
        )
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.set_aspect('equal')

    # Hide unused subplots
    for idx in range(n_clusters, len(axes)):
        axes[idx].axis('off')

    fig.suptitle(
        "Mapping Rate vs Sparsity Distributions by Cluster\n(Hexbin Heatmaps, Mate 1)",
        fontsize=16,
        y=0.995
    )

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"\nSaved sparsity cluster hexbin plots to {output_path}")

    plt.savefig(output_path.with_suffix('.png'), dpi=300, bbox_inches='tight')
    plt.close()


def main():
    OUTPUT_DIR.mkdir(exist_ok=True)

    print("Loading metrics...")
    df = load_metrics()

    # Print summary
    print(f"\nDataset summary:")
    print(f"  Total runs: {df.height}")
    print(f"  Single-end: {df.filter(pl.col('n_mates') == 1).height}")
    print(f"  Paired-end: {df.filter(pl.col('n_mates') == 2).height}")

    if "bioproject" in df.columns:
        n_projects = df.filter(pl.col("bioproject").is_not_null()).select("bioproject").n_unique()
        print(f"  Bioprojects: {n_projects}")

    print("\nCreating mate 1 vs mate 2 contour plot...")
    plot_mate1_vs_mate2_mapping(df, OUTPUT_DIR / "contour_mate1_vs_mate2.svg")

    print("\nCreating mapping vs sparsity contour plot...")
    plot_mapping_vs_sparsity(df, OUTPUT_DIR / "contour_mapping_vs_sparsity.svg")

    print("\nPerforming hierarchical clustering of mate1-mate2 distributions...")
    print("Using Wasserstein distance to compare full distributions...")
    cluster_result = cluster_with_wasserstein_and_plot(
        df,
        OUTPUT_DIR / "dendrogram_mate1_mate2_wasserstein.svg",
        n_clusters=8,
        min_samples=50,  # Increased from 5 to reduce number of bioprojects
        wasserstein_method="sliced",
        n_projections=50,  # Reduced from 100 for speed
        max_samples_per_bioproject=500  # Reduced for speed
    )

    if cluster_result is not None:
        print(f"\nMate1-Mate2 clustering complete!")
        print(f"Results include:")
        print(f"  - Dendrogram: {OUTPUT_DIR / 'dendrogram_mate1_mate2_wasserstein.svg'}")
        print(f"  - Cluster assignments: {OUTPUT_DIR / 'cluster_assignments.tsv'}")

        # Plot hexbins for each cluster
        print("\nCreating hexbin heatmaps for mate1-mate2 clusters...")
        plot_cluster_hexbins(
            df,
            cluster_result['cluster_labels'],
            cluster_result['bioproject_names'],
            OUTPUT_DIR / "cluster_hexbins_wasserstein.svg"
        )

    # Clustering for mapping rate vs sparsity
    print("\n" + "="*80)
    print("MAPPING RATE vs SPARSITY CLUSTERING")
    print("="*80)
    print("\nPerforming hierarchical clustering of mapping vs sparsity distributions...")
    print("Using Wasserstein distance to compare full distributions...")

    sparsity_cluster_result = cluster_sparsity_with_wasserstein(
        df,
        OUTPUT_DIR / "dendrogram_mapping_sparsity_wasserstein.svg",
        n_clusters=8,
        min_samples=50,
        wasserstein_method="sliced",
        n_projections=50,
        max_samples_per_bioproject=500
    )

    if sparsity_cluster_result is not None:
        print(f"\nMapping vs sparsity clustering complete!")
        print(f"Results include:")
        print(f"  - Dendrogram: {OUTPUT_DIR / 'dendrogram_mapping_sparsity_wasserstein.svg'}")
        print(f"  - Cluster assignments: {OUTPUT_DIR / 'cluster_assignments_sparsity.tsv'}")

        # Plot hexbins for sparsity clusters
        print("\nCreating hexbin heatmaps for sparsity clusters...")
        plot_sparsity_cluster_hexbins(
            df,
            sparsity_cluster_result['cluster_labels'],
            sparsity_cluster_result['bioproject_names'],
            OUTPUT_DIR / "cluster_hexbins_sparsity_wasserstein.svg"
        )

    print("\nDone!")


if __name__ == "__main__":
    main()
