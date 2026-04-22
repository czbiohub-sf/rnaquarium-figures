#!/usr/bin/env python3
"""
Figure 1 Panel B.3: Read/mate distributions (mapping) - Seq-Detective filtering view

Creates hexbin plots of mapped reads (mate1 vs mate2) colored by seq-detective's
per-mate filtering decision. Shows why seq-detective is needed before RNA-seq processing.

Categories: B-B, B-T, T-B, T-T (where B=Biological, T=Technical)
"""

import polars as pl
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# =============================================================================
# Configuration
# =============================================================================

SEQDETECTIVE_METRICS = Path("data/75k_unstable/seqdetective_metrics.parquet")
SEQDETECTIVE_JUDGEMENT = Path("data/75k_unstable/seq-detective-judgement-summary-augmented.txt")
OUTPUT_DIR = Path("figures")
OUTPUT_DIR.mkdir(exist_ok=True)

# Custom categorical colors
BATLOW_COLORS = {
    "B-B": "#99882c",  # olive/yellow-green  (both biological, most common)
    "B-T": "#0b2c5c",  # dark navy           (M1 bio, M2 tech, minor issue)
    "T-B": "#426f52",  # forest green        (M1 tech, M2 bio, common issue)
    "T-T": "#f29d6c",  # coral/peach         (both technical, problem)
}

# Figure parameters
CM = 1 / 2.54
FIG_SIZE_MAIN = (3.75 * CM, 3.75 * CM)  # panel B3: 3.25 × 3.25 cm
FIG_SIZE_GRID = (12, 12)                  # supplemental 2×2 grid
FIG_SIZE_PANEL = (16, 10)                 # supplemental 2×3 panel
HEXBIN_GRIDSIZE = 50
DPI = 300
MARKER_SIZE = 1.0                         # ~0.2cm marker diameter
MARKER_ALPHA = 0.5                        # marker opacity
AXIS_COLOR = '#585858'                    # axes and ticks color

# Font settings (Arial, 5–7pt per layout spec)
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial', 'Helvetica Neue', 'DejaVu Sans']
plt.rcParams['svg.fonttype'] = 'none'   # SVG text as text, not paths
FONT_TITLE = 4
FONT_LABEL = 5
FONT_LEGEND = 4
FONT_TICK = 4.5

# =============================================================================
# Load data
# =============================================================================

# Load mapping metrics
metrics_df = pl.read_parquet(SEQDETECTIVE_METRICS)
print(f"Loaded {metrics_df.height} runs with seq-detective metrics")

# Load seq-detective filtering decisions
# Format for PE: ID \t mate1_file \t mate2_file \t grade1 \t grade2 \t reason
# Format for SE: ID \t file \t \t grade \t \t reason
judgement_df = pl.read_csv(
    SEQDETECTIVE_JUDGEMENT,
    separator="\t",
    has_header=False,
    new_columns=["id", "file1", "file2", "grade1", "grade2", "reason"]
)

# Strip "_subsample" suffix from IDs to match metrics
judgement_df = judgement_df.with_columns(
    id=pl.col("id").str.replace("_subsample", "")
)

print(f"Loaded {judgement_df.height} seq-detective judgements")

# Filter to PE samples (have non-null file2 and grade2)
pe_judgement = judgement_df.filter(pl.col("file2").is_not_null()).select([
    "id", "grade1", "grade2", "reason"
])

print(f"\nPE judgements: {pe_judgement.height}")

# Show mate grade combinations
grade_combos = (
    pe_judgement
    .group_by(["grade1", "grade2"])
    .agg(pl.len().alias("count"))
    .sort("count", descending=True)
)
print("\nMate grade combinations:")
print(grade_combos)

# Merge metrics with judgements
merged_df = metrics_df.join(pe_judgement, on="id", how="left")

# Filter to PE samples with judgements
pe_df = merged_df.filter(
    (pl.col("n_mates") == 2) &
    (pl.col("grade1").is_not_null()) &
    (pl.col("grade2").is_not_null())
)
print(f"\nPE samples with metrics and judgements: {pe_df.height}")

# Create combined category column
pe_df = pe_df.with_columns(
    category=pl.concat_str([pl.col("grade1"), pl.lit("-"), pl.col("grade2")])
)

category_counts = (
    pe_df
    .group_by("category")
    .agg(pl.len().alias("count"))
    .sort("count", descending=True)
)
print("\nCategory distribution:")
print(category_counts)

# =============================================================================
# Create scatter plot colored by mate categories
# =============================================================================

def create_category_scatter_plot(data: pl.DataFrame, output_path: Path):
    """Create scatter plot showing mate1 vs mate2 mapping rates, colored by mate categories."""

    plot_data = (
        data.select(["mapping_rate_m1", "mapping_rate_m2", "category"])
        .drop_nulls(subset=["mapping_rate_m1", "mapping_rate_m2"])
        .to_pandas()
    )

    fig, ax = plt.subplots(figsize=FIG_SIZE_MAIN)

    # Define category order (plot in order so T-T is on top)
    categories = [
        ("B-B", 1, 0.7),
        ("B-T", 2, 0.7),
        ("T-B", 3, 0.7),
        ("T-T", 4, 0.7),
    ]

    # Plot each category
    for category, zorder, alpha in categories:
        if category not in plot_data["category"].values:
            continue
        subset = plot_data[plot_data["category"] == category]

        ax.scatter(
            subset["mapping_rate_m1"],
            subset["mapping_rate_m2"],
            c=BATLOW_COLORS[category],
            s=MARKER_SIZE,
            alpha=MARKER_ALPHA,
            linewidths=0.05,
            #edgecolors="black",
            zorder=zorder,
            rasterized=True,
            #edgecolors='none'
        )

    # Add diagonal reference line
    ax.plot([0, 1], [0, 1], 'k--', alpha=0.3, linewidth=1, zorder=0)

    # Formatting - Nature style
    ax.set_xlabel("Mate 1 Mapping Rate", fontsize=FONT_LABEL)
    ax.set_ylabel("Mate 2 Mapping Rate", fontsize=FONT_LABEL)
    ax.set_title(
        "Seq-Detective Mate Pair Filtering Outcomes",
        fontsize=FONT_TITLE
    )
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.tick_params(labelsize=FONT_TICK, direction='out', colors=AXIS_COLOR)

    # Grid: only at 0.5 / 0.5
    ax.set_xticks([0, 0.5, 1.0])
    ax.set_yticks([0, 0.5, 1.0])
    ax.grid(True, alpha=0.3, linewidth=0.5, which='minor')
    ax.set_axisbelow(True)

    # Spines: remove top and right, style others
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_linewidth(0.5)
    ax.spines['left'].set_color(AXIS_COLOR)
    ax.spines['bottom'].set_linewidth(0.5)
    ax.spines['bottom'].set_color(AXIS_COLOR)

    fig.tight_layout(pad=0)
    plt.savefig(output_path, dpi=DPI)
    print(f"\nSaved category scatter plot to {output_path}")
    plt.close()

create_category_scatter_plot(pe_df, OUTPUT_DIR / "Fig1_B_3_mapping_seqdetective_categories.svg")

# =============================================================================
# Create hexbin density plots for each category (2x2 grid)
# =============================================================================

def create_category_hexbin_grid(data: pl.DataFrame, output_path: Path):
    """Create 2x2 grid of hexbin density plots for each mate category."""

    plot_data = (
        data.select(["mapping_rate_m1", "mapping_rate_m2", "category"])
        .drop_nulls(subset=["mapping_rate_m1", "mapping_rate_m2"])
        .to_pandas()
    )

    fig, axes = plt.subplots(2, 2, figsize=FIG_SIZE_GRID)

    categories = [
        ("B-B", "both biological"),
        ("B-T", "M1 biological, M2 technical"),
        ("T-B", "M1 technical, M2 biological"),
        ("T-T", "both technical")
    ]

    for idx, (category, description) in enumerate(categories):
        ax = axes[idx // 2, idx % 2]
        subset = plot_data[plot_data["category"] == category]

        if len(subset) > 0:
            # Hexbin with log scale - consistent gridsize
            hexbin = ax.hexbin(
                subset["mapping_rate_m1"],
                subset["mapping_rate_m2"],
                gridsize=HEXBIN_GRIDSIZE,
                cmap='viridis',
                mincnt=1,
                bins='log',
                alpha=0.9
            )

            # Add diagonal line
            ax.plot([0, 1], [0, 1], 'k--', alpha=0.3, linewidth=1)

            # Formatting - Nature style
            ax.set_xlabel("Mate 1 Mapping Rate", fontsize=FONT_LABEL)
            ax.set_ylabel("Mate 2 Mapping Rate", fontsize=FONT_LABEL)
            ax.set_title(
                f"{category}: {description} (n={len(subset):,})",
                fontsize=FONT_TITLE,
                color=BATLOW_COLORS[category]
            )
            ax.set_xlim(0, 1)
            ax.set_ylim(0, 1)
            ax.set_aspect('equal')  # Ensure square aspect ratio
            ax.tick_params(labelsize=FONT_TICK, direction='out', colors=AXIS_COLOR)

            # Grid: only at 0.5 / 0.5
            ax.set_xticks([0, 0.5, 1.0])
            ax.set_yticks([0, 0.5, 1.0])
            ax.grid(True, alpha=0.3, linewidth=0.5, which='minor')
            ax.set_axisbelow(True)

            # Spines: remove top and right, style others
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            ax.spines['left'].set_linewidth(0.5)
            ax.spines['left'].set_color(AXIS_COLOR)
            ax.spines['bottom'].set_linewidth(0.5)
            ax.spines['bottom'].set_color(AXIS_COLOR)

            # Colorbar
            cbar = plt.colorbar(hexbin, ax=ax)
            cbar.set_label('Sample count (log scale)', fontsize=FONT_LEGEND)
            cbar.ax.tick_params(labelsize=FONT_TICK)

    fig.tight_layout(pad=0)
    plt.savefig(output_path, dpi=DPI)
    print(f"Saved category hexbin grid to {output_path}")
    plt.close()

create_category_hexbin_grid(pe_df, OUTPUT_DIR / "Fig1_B_3_mapping_seqdetective_density_grid.svg")

# =============================================================================
# Create reason breakdown for technical classifications (supplemental)
# =============================================================================

def create_technical_reason_breakdown(data: pl.DataFrame, output_path: Path):
    """Show distributions by filtering reason for samples with technical mates."""

    # Filter to samples with at least one technical mate
    technical_data = (
        data.filter(
            (pl.col("grade1") == "T") | (pl.col("grade2") == "T")
        )
        .select(["mapping_rate_m1", "mapping_rate_m2", "reason", "category"])
        .drop_nulls(subset=["mapping_rate_m1", "mapping_rate_m2"])
        .to_pandas()
    )

    # Get top reasons
    reason_counts = technical_data["reason"].value_counts().head(6)
    top_reasons = reason_counts.index.tolist()

    fig, axes = plt.subplots(2, 3, figsize=FIG_SIZE_PANEL)
    axes = axes.flatten()

    for i, reason in enumerate(top_reasons):
        ax = axes[i]
        subset = technical_data[technical_data["reason"] == reason]

        if len(subset) > 0:
            # Plot by category in consistent order
            for category in ["T-B", "B-T", "T-T"]:
                cat_subset = subset[subset["category"] == category]
                if len(cat_subset) > 0:
                    # Short label format: "B T (n=...)"
                    label = f"{category.replace('-', ' ')} (n={len(cat_subset):,})"
                    ax.scatter(
                        cat_subset["mapping_rate_m1"],
                        cat_subset["mapping_rate_m2"],
                        c=BATLOW_COLORS[category],
                        s=MARKER_SIZE,
                        alpha=MARKER_ALPHA,
                        label=label,
                        rasterized=True,
                        edgecolors='none'
                    )

            # Diagonal line
            ax.plot([0, 1], [0, 1], 'k--', alpha=0.3, linewidth=1)

            # Formatting - Nature style
            ax.set_xlabel("Mate 1 Mapping Rate", fontsize=FONT_LABEL)
            ax.set_ylabel("Mate 2 Mapping Rate", fontsize=FONT_LABEL)
            ax.set_title(
                f"{reason} (n={len(subset):,})",
                fontsize=FONT_TITLE
            )
            ax.set_xlim(0, 1)
            ax.set_ylim(0, 1)
            ax.tick_params(labelsize=FONT_TICK, direction='out', colors=AXIS_COLOR)

            # Grid: only at 0.5 / 0.5
            ax.set_xticks([0, 0.5, 1.0])
            ax.set_yticks([0, 0.5, 1.0])
            ax.grid(True, alpha=0.3, linewidth=0.5, which='minor')
            ax.set_axisbelow(True)

            # Spines: remove top and right, style others
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            ax.spines['left'].set_linewidth(0.5)
            ax.spines['left'].set_color(AXIS_COLOR)
            ax.spines['bottom'].set_linewidth(0.5)
            ax.spines['bottom'].set_color(AXIS_COLOR)

            # Legend with larger markers and semitransparent background
            legend = ax.legend(
                fontsize=FONT_LEGEND,
                loc='upper left',
                frameon=True,
                framealpha=0.7,
                edgecolor='black',
                fancybox=False,
                markerscale=3.0
            )
            legend.get_frame().set_linewidth(0.5)

    fig.tight_layout(pad=0)
    plt.savefig(output_path, dpi=DPI)
    print(f"Saved technical reason breakdown to {output_path}")
    plt.close()

create_technical_reason_breakdown(pe_df, OUTPUT_DIR / "Fig1_B_3_mapping_seqdetective_technical_reasons.svg")

print("\nDone!")
