#!/usr/bin/env python3
"""
Figure 1 Panel B.3: Read/mate distributions (mapping)

Creates hexbin/contour plots of mapped reads (mate1 vs mate2) from seq-detective metrics.
Clusters by bioproject distributions, draws all ~77k samples with log scale density.
"""

import polars as pl
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from scipy.cluster import hierarchy
from scipy.spatial.distance import pdist
from scipy.stats import wasserstein_distance

# =============================================================================
# Configuration
# =============================================================================

SEQDETECTIVE_METRICS = Path("data/75k_unstable/seqdetective_metrics.parquet")
OUTPUT_DIR = Path("figures")
OUTPUT_DIR.mkdir(exist_ok=True)

N_CLUSTERS = 8  # Expected number of meaningful clusters
MIN_BIOPROJECT_SIZE = 50  # Minimum runs per bioproject for clustering

# =============================================================================
# Load data
# =============================================================================

df = pl.read_parquet(SEQDETECTIVE_METRICS)
print(f"Loaded {df.height} runs with seq-detective metrics")

# Filter to paired-end samples only
pe_df = df.filter(pl.col("n_mates") == 2)
print(f"Paired-end samples: {pe_df.height}")

# =============================================================================
# Cluster bioprojects by mapping rate distributions
# =============================================================================

def cluster_bioprojects(data: pl.DataFrame, min_size: int = MIN_BIOPROJECT_SIZE):
    """
    Cluster bioprojects by their mate1 vs mate2 mapping rate distributions.
    Uses Wasserstein distance between bioproject distributions.
    """
    # Filter to bioprojects with enough samples
    bioproject_sizes = (
        data.filter(pl.col("bioproject").is_not_null())
        .group_by("bioproject")
        .agg(pl.len().alias("count"))
    )
    large_bioprojects = (
        bioproject_sizes
        .filter(pl.col("count") >= min_size)
        .sort("count", descending=True)
    )

    print(f"\nFound {len(large_bioprojects)} bioprojects with >={min_size} samples")
    print(large_bioprojects.head(20))

    # Extract mapping rate distributions for each bioproject
    bioprojects = large_bioprojects["bioproject"].to_list()
    distributions = {}

    for bp in bioprojects:
        bp_data = (
            data.filter(pl.col("bioproject") == bp)
            .select(["mapping_rate_m1", "mapping_rate_m2"])
            .drop_nulls()
        )
        if bp_data.height >= min_size:
            distributions[bp] = bp_data.to_pandas()

    print(f"Computing distances for {len(distributions)} bioprojects...")

    # Compute pairwise Wasserstein distances
    bp_list = list(distributions.keys())
    n = len(bp_list)
    dist_matrix = np.zeros((n, n))

    for i in range(n):
        for j in range(i+1, n):
            bp1, bp2 = bp_list[i], bp_list[j]
            d1, d2 = distributions[bp1], distributions[bp2]

            # Wasserstein distance on mate1 mapping rates
            w1 = wasserstein_distance(
                d1["mapping_rate_m1"].dropna(),
                d2["mapping_rate_m1"].dropna()
            )
            # Wasserstein distance on mate2 mapping rates
            w2 = wasserstein_distance(
                d1["mapping_rate_m2"].dropna(),
                d2["mapping_rate_m2"].dropna()
            )

            dist = (w1 + w2) / 2
            dist_matrix[i, j] = dist
            dist_matrix[j, i] = dist

    # Hierarchical clustering
    condensed_dist = dist_matrix[np.triu_indices(n, k=1)]
    linkage = hierarchy.linkage(condensed_dist, method='ward')

    # Cut tree to get cluster assignments
    cluster_labels = hierarchy.fcluster(linkage, t=N_CLUSTERS, criterion='maxclust')

    # Map bioprojects to clusters
    bioproject_clusters = dict(zip(bp_list, cluster_labels))

    return bioproject_clusters, linkage, bp_list, dist_matrix

# Cluster bioprojects
bioproject_clusters, linkage, bp_list, dist_matrix = cluster_bioprojects(pe_df)

# Add cluster labels to data
pe_with_clusters = pe_df.with_columns(
    cluster=pl.col("bioproject").replace_strict(
        bioproject_clusters,
        default=None,
        return_dtype=pl.Int32
    )
)

print(f"\nCluster distribution:")
cluster_counts = (
    pe_with_clusters
    .filter(pl.col("cluster").is_not_null())
    .group_by("cluster")
    .agg(pl.len().alias("count"))
    .sort("cluster")
)
print(cluster_counts)

# =============================================================================
# Create main hexbin plot with all samples
# =============================================================================

def create_main_hexbin_plot(data: pl.DataFrame, output_path: Path):
    """Create main hexbin plot showing all samples with cluster overlays."""

    # Convert to pandas for plotting
    plot_data = (
        data.select(["mapping_rate_m1", "mapping_rate_m2", "cluster"])
        .drop_nulls()
        .to_pandas()
    )

    fig, ax = plt.subplots(figsize=(10, 10))

    # Hexbin of all samples (log scale density)
    hexbin = ax.hexbin(
        plot_data["mapping_rate_m1"],
        plot_data["mapping_rate_m2"],
        gridsize=50,
        cmap='Blues',
        mincnt=1,
        bins='log',
        alpha=0.6
    )

    # Overlay cluster contours
    for cluster_id in sorted(plot_data["cluster"].dropna().unique()):
        cluster_data = plot_data[plot_data["cluster"] == cluster_id]

        if len(cluster_data) > 10:
            sns.kdeplot(
                data=cluster_data,
                x="mapping_rate_m1",
                y="mapping_rate_m2",
                levels=1,
                linewidths=2,
                alpha=0.8,
                ax=ax,
                label=f"Cluster {int(cluster_id)}"
            )

    # Add diagonal reference line
    ax.plot([0, 1], [0, 1], 'k--', alpha=0.3, linewidth=2)

    # Formatting
    ax.set_xlabel("Mate 1 Mapping Rate", fontsize=14, fontweight='bold')
    ax.set_ylabel("Mate 2 Mapping Rate", fontsize=14, fontweight='bold')
    ax.set_title(
        "Paired-End Mapping Rates Distribution\n(Hexbin density + Cluster contours)",
        fontsize=16,
        fontweight='bold',
        pad=20
    )
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.legend(loc='upper left', fontsize=10)

    # Colorbar
    cbar = plt.colorbar(hexbin, ax=ax)
    cbar.set_label('Sample count (log scale)', fontsize=12)

    fig.tight_layout(pad=0)
    plt.savefig(output_path, dpi=300)
    print(f"\nSaved main plot to {output_path}")
    plt.close()

create_main_hexbin_plot(pe_with_clusters, OUTPUT_DIR / "Fig1_B_3_mapping_mate1_vs_mate2_hexbin.svg")

# =============================================================================
# Create dendrogram
# =============================================================================

def create_dendrogram(linkage, labels, output_path: Path):
    """Create dendrogram showing bioproject clustering."""

    fig, ax = plt.subplots(figsize=(12, 8))

    hierarchy.dendrogram(
        linkage,
        labels=labels,
        ax=ax,
        leaf_font_size=8,
        orientation='right'
    )

    ax.set_xlabel("Distance", fontsize=14, fontweight='bold')
    ax.set_ylabel("BioProject", fontsize=14, fontweight='bold')
    ax.set_title(
        "Hierarchical Clustering of BioProjects\n(by Mate1 vs Mate2 Mapping Distribution)",
        fontsize=16,
        fontweight='bold',
        pad=20
    )

    fig.tight_layout(pad=0)
    plt.savefig(output_path, dpi=300)
    print(f"Saved dendrogram to {output_path}")
    plt.close()

create_dendrogram(linkage, bp_list, OUTPUT_DIR / "Fig1_B_3_mapping_dendrogram.svg")

# =============================================================================
# Create separate plots for each cluster (supplemental)
# =============================================================================

def create_cluster_panels(data: pl.DataFrame, output_path: Path):
    """Create multi-panel plot showing each cluster separately."""

    # Get cluster data
    cluster_data = (
        data.filter(pl.col("cluster").is_not_null())
        .select(["mapping_rate_m1", "mapping_rate_m2", "cluster"])
        .drop_nulls()
        .to_pandas()
    )

    clusters = sorted(cluster_data["cluster"].unique())
    n_clusters = len(clusters)

    # Create grid layout
    n_cols = 3
    n_rows = int(np.ceil(n_clusters / n_cols))

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(15, 5 * n_rows))
    axes = axes.flatten() if n_clusters > 1 else [axes]

    for i, cluster_id in enumerate(clusters):
        ax = axes[i]
        cluster_subset = cluster_data[cluster_data["cluster"] == cluster_id]

        # Hexbin for this cluster
        ax.hexbin(
            cluster_subset["mapping_rate_m1"],
            cluster_subset["mapping_rate_m2"],
            gridsize=30,
            cmap='Blues',
            mincnt=1,
            bins='log'
        )

        # Diagonal line
        ax.plot([0, 1], [0, 1], 'k--', alpha=0.3)

        ax.set_xlabel("Mate 1 Mapping Rate", fontsize=10)
        ax.set_ylabel("Mate 2 Mapping Rate", fontsize=10)
        ax.set_title(f"Cluster {int(cluster_id)} (n={len(cluster_subset):,})", fontsize=12, fontweight='bold')
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)

    # Remove empty subplots
    for j in range(i+1, len(axes)):
        fig.delaxes(axes[j])

    plt.suptitle("Mapping Rate Distributions by Cluster", fontsize=16, fontweight='bold', y=1.00)
    fig.tight_layout(pad=0)
    plt.savefig(output_path, dpi=300)
    print(f"Saved cluster panels to {output_path}")
    plt.close()

create_cluster_panels(pe_with_clusters, OUTPUT_DIR / "Fig1_B_3_mapping_cluster_panels.svg")

print("\nDone!")
