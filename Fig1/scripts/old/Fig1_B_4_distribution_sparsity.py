#!/usr/bin/env python3
"""
ABANDONED (archived for reference): sparsity vs mapping-ratio scatter

Original Panel B.4 candidate, marked NOT USING in design.md — the sparsity
view did not yield a compelling narrative next to the Seq-Detective
classification scatter. Kept here only so the plotting recipe (including the
min-sparsity computation across mates and the technology-colored variant)
stays on record.
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
TECH_ANNOTATIONS = Path("data/zf-core-v2-74K_problematic_with_bulk.kmers.allcols.csv")
OUTPUT_DIR = Path("figures")
OUTPUT_DIR.mkdir(exist_ok=True)

# Custom categorical colors for seq-detective categories (matching Fig1_B_3)
SEQDETECTIVE_COLORS = {
    "B-B": "#99882c",  # Olive/yellow-green (both biological, most common, good)
    "B-T": "#0b2c5c",  # Dark navy (M1 bio, M2 tech, minor issue)
    "T-B": "#426f52",  # Forest green/teal (M1 tech, M2 bio, common issue)
    "T-T": "#f29d6c",  # Coral/peach (both technical, problem)
}

# Load batlow categorical palette for technology colors
def load_batlow_categorical_palette():
    """Load batlow categorical color palette from .txt file."""
    palette_file = Path("palette/batlow/CategoricalPalettes/batlowS.txt")
    colors = []
    with open(palette_file) as f:
        for line in f:
            line = line.strip()
            if line:
                # Parse RGB values (0-1 range)
                try:
                    r, g, b = map(float, line.split())
                    # Convert to hex
                    hex_color = f"#{int(r*255):02x}{int(g*255):02x}{int(b*255):02x}"
                    colors.append(hex_color)
                except ValueError:
                    continue
    return colors

BATLOW_CATEGORICAL = load_batlow_categorical_palette()
print(f"Loaded {len(BATLOW_CATEGORICAL)} categorical colors from batlow palette")

# Figure parameters
CM = 1 / 2.54
FIG_SIZE_MAIN = (5.5 * CM, 6.178 * CM)  # panel B4: 5.5 × 6.178 cm
HEXBIN_GRIDSIZE = 50
DPI = 300

# Font settings (Arial, 5–7pt per layout spec)
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial', 'Helvetica Neue', 'DejaVu Sans']
FONT_TITLE = 7
FONT_LABEL = 7
FONT_LEGEND = 6
FONT_TICK = 5

# =============================================================================
# Load data
# =============================================================================

df = pl.read_parquet(SEQDETECTIVE_METRICS)
print(f"Loaded {df.height} runs with seq-detective metrics")

# Load seq-detective filtering decisions
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

# Merge with metrics
df = df.join(judgement_df.select(["id", "grade1", "grade2"]), on="id", how="left")

# Load manual technology annotations
# Handle Excel error values and only load needed columns
tech_df = pl.read_csv(
    TECH_ANNOTATIONS,
    null_values=["#DIV/0!", ""],
    ignore_errors=True
)
print(f"Loaded {tech_df.height} technology annotations")

# Select relevant columns and rename for joining
# Use the first 'technology' column (column 6)
tech_df = (
    tech_df.select([
        pl.col("accession").alias("id"),
        pl.col("technology").alias("technology_raw")
    ])
    .with_columns(
        # Clean up technology values
        pl.col("technology_raw")
        .str.strip_chars()
        .str.to_uppercase()
        .alias("technology")
    )
    .select(["id", "technology"])
)

# Merge technology annotations
df = df.join(tech_df, on="id", how="left")
print(f"Merged technology annotations: {df.filter(pl.col('technology').is_not_null()).height} runs with technology labels")

# =============================================================================
# Process sparsity data
# =============================================================================

# For paired-end: calculate mapping ratio and min sparsity
pe_df = (
    df.filter(pl.col("n_mates") == 2)
    .with_columns([
        (pl.col("mapping_rate_m1") / pl.col("mapping_rate_m2")).alias("mapping_ratio"),
        pl.min_horizontal("sparsity_m1", "sparsity_m2").alias("min_sparsity"),
        pl.concat_str([
            pl.col("grade1").fill_null("U"),
            pl.lit("-"),
            pl.col("grade2").fill_null("U")
        ]).alias("category")
    ])
)

# For single-end: use mapping rate directly
se_df = (
    df.filter(pl.col("n_mates") == 1)
    .with_columns([
        pl.col("mapping_rate").alias("mapping_ratio"),
        pl.col("sparsity").alias("min_sparsity"),
        pl.concat_str([pl.col("grade1").fill_null("U"), pl.lit("-SE")]).alias("category")
    ])
)

# Combine
combined = pl.concat([pe_df, se_df])
print(f"Total samples: {combined.height} (PE: {pe_df.height}, SE: {se_df.height})")

# Show category distribution
category_counts = (
    combined
    .group_by("category")
    .agg(pl.len().alias("count"))
    .sort("count", descending=True)
)
print("\nCategory distribution:")
print(category_counts.head(10))

# =============================================================================
# Seq-detective category scatter plot
# =============================================================================

def create_seqdetective_category_scatter(data: pl.DataFrame, output_path: Path):
    """Create scatter plot colored by seq-detective mate categories."""

    # Filter to PE samples with valid categories
    plot_data = (
        data.filter(
            (pl.col("n_mates") == 2) &
            (pl.col("category").is_in(["B-B", "B-T", "T-B", "T-T"]))
        )
        .select(["mapping_ratio", "min_sparsity", "category"])
        .drop_nulls()
        .filter(
            (pl.col("mapping_ratio") > 0.01) &
            (pl.col("mapping_ratio") < 100) &
            (pl.col("min_sparsity") > 0) &
            (pl.col("min_sparsity") <= 1)
        )
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

        label = f"{category.replace('-', ' ')} (n={len(subset):,})"

        ax.scatter(
            subset["mapping_ratio"],
            subset["min_sparsity"],
            c=SEQDETECTIVE_COLORS[category],
            s=10,
            alpha=alpha,
            label=label,
            zorder=zorder,
            rasterized=True,
            edgecolors='none'
        )

    # Add reference line at ratio = 1
    ax.axvline(1, color='k', linestyle='--', alpha=0.3, linewidth=1, zorder=0)

    # Formatting - Nature style with log x-axis, linear y-axis
    ax.set_xlabel("Mapping Ratio (Mate1/Mate2, log scale)", fontsize=FONT_LABEL)
    ax.set_ylabel("Minimum Sparsity", fontsize=FONT_LABEL)
    ax.set_title(
        "Gene Sparsity by Seq-Detective Filtering Decision",
        fontsize=FONT_TITLE
    )
    ax.set_xlim(0.01, 100)
    ax.set_ylim(0.5, 1)
    ax.set_xscale('log')
    ax.tick_params(labelsize=FONT_TICK)

    # Legend with larger markers and semitransparent background
    legend = ax.legend(
        loc='lower left',
        fontsize=FONT_LEGEND,
        frameon=True,
        framealpha=0.7,
        edgecolor='black',
        fancybox=False,
        markerscale=2.5
    )
    legend.get_frame().set_linewidth(0.5)

    # Clean grid
    ax.grid(True, alpha=0.15, linewidth=0.5)
    ax.set_axisbelow(True)

    # Spines
    for spine in ax.spines.values():
        spine.set_linewidth(0.5)

    fig.tight_layout(pad=0)
    plt.savefig(output_path, dpi=DPI)
    print(f"\nSaved seq-detective category scatter to {output_path}")
    plt.close()

create_seqdetective_category_scatter(combined, OUTPUT_DIR / "Fig1_B_4_sparsity_seqdetective_categories.svg")

# =============================================================================
# Technology category scatter plot
# =============================================================================

def create_technology_scatter(data: pl.DataFrame, output_path: Path):
    """Create scatter plot colored by technology category."""

    # Filter to samples with valid technology and metrics
    plot_data = (
        data.select(["mapping_ratio", "min_sparsity", "technology"])
        .filter(pl.col("technology").is_not_null())
        .filter(
            (pl.col("mapping_ratio") > 0.01) &
            (pl.col("mapping_ratio") < 100) &
            (pl.col("min_sparsity") > 0) &
            (pl.col("min_sparsity") <= 1)
        )
    )

    if plot_data.height == 0:
        print("\nWarning: No samples with technology annotations. Skipping technology scatter plot.")
        return

    # Group rare technologies and "unknown" as "Other"
    tech_counts = (
        plot_data
        .group_by("technology")
        .agg(pl.len().alias("count"))
        .sort("count", descending=True)
    )

    # Keep top technologies, group rest as Other
    # Also treat "unknown" and empty strings as "Other"
    top_technologies = [
        t for t in tech_counts.head(8)["technology"].to_list()
        if t not in ["unknown", "", "Unknown"]
    ]

    plot_data = plot_data.with_columns(
        pl.when(
            pl.col("technology").is_in(top_technologies) &
            ~pl.col("technology").is_in(["unknown", "", "Unknown"])
        )
        .then(pl.col("technology"))
        .otherwise(pl.lit("Other"))
        .alias("technology_grouped")
    ).to_pandas()

    print(f"\nTechnology distribution:")
    print(plot_data["technology_grouped"].value_counts())

    fig, ax = plt.subplots(figsize=FIG_SIZE_MAIN)

    # Get unique technologies and sort by frequency
    tech_order = plot_data["technology_grouped"].value_counts().index.tolist()

    # Assign colors from batlow categorical palette
    tech_colors = {}
    for i, tech in enumerate(tech_order):
        tech_colors[tech] = BATLOW_CATEGORICAL[i % len(BATLOW_CATEGORICAL)]

    # Plot each technology (least frequent first so most frequent is on top)
    for i, tech in enumerate(reversed(tech_order)):
        subset = plot_data[plot_data["technology_grouped"] == tech]
        color = tech_colors[tech]

        label = f"{tech} (n={len(subset):,})"

        ax.scatter(
            subset["mapping_ratio"],
            subset["min_sparsity"],
            c=color,
            s=10,
            alpha=0.6,
            label=label,
            zorder=i+1,
            rasterized=True,
            edgecolors='none'
        )

    # Add reference line at ratio = 1
    ax.axvline(1, color='k', linestyle='--', alpha=0.3, linewidth=1, zorder=0)

    # Formatting - Nature style with log x-axis, linear y-axis
    ax.set_xlabel("Mapping Ratio (Mate1/Mate2, log scale)", fontsize=FONT_LABEL)
    ax.set_ylabel("Minimum Sparsity", fontsize=FONT_LABEL)
    ax.set_title(
        "Gene Sparsity by Sequencing Technology",
        fontsize=FONT_TITLE
    )
    ax.set_xlim(0.01, 100)
    ax.set_ylim(0.5, 1)
    ax.set_xscale('log')
    ax.tick_params(labelsize=FONT_TICK)

    # Legend with larger markers and semitransparent background
    legend = ax.legend(
        loc='lower left',
        fontsize=FONT_LEGEND,
        frameon=True,
        framealpha=0.7,
        edgecolor='black',
        fancybox=False,
        markerscale=2.5,
        ncol=1
    )
    legend.get_frame().set_linewidth(0.5)

    # Clean grid
    ax.grid(True, alpha=0.15, linewidth=0.5)
    ax.set_axisbelow(True)

    # Spines
    for spine in ax.spines.values():
        spine.set_linewidth(0.5)

    fig.tight_layout(pad=0)
    plt.savefig(output_path, dpi=DPI)
    print(f"\nSaved technology scatter to {output_path}")
    plt.close()

create_technology_scatter(combined, OUTPUT_DIR / "Fig1_B_4_sparsity_technology_categories.svg")

print("\nDone!")
