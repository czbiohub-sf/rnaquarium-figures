#!/usr/bin/env python3
"""
Figure 1 Panel B.1: Technology and Seq-Detective outcome pies

Two polar bar plots for Panel B, drawn at different physical sizes:
  1. `Fig1_B_1_tech_pie.svg` — simple technology-category pie covering all
     accessions (including "Other"); no outer ring.
  2. `Fig1_B_1_seqdetective_pie.svg` — multi-level pie on the named-tech
     subset. Inner ring is the same technology category; outer ring is the
     Seq-Detective per-accession outcome (BB / TB / BT / TT / B / T).

Both charts share the technology-color mapping (computed once from the full
batlow palette across all tech categories, count-sorted). Outcome-ring colors
stay in sync with `Fig1_B_2_seqdetective_scatter.py`.
"""

import polars as pl
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# =============================================================================
# Configuration
# =============================================================================

CM = 1 / 2.54          # inches per cm

# Both figures are square; the polar circle fills the canvas.
FIG_SIMPLE = (2.34 * CM, 1.72 * CM)   # simple chart: 1.72 cm diameter
FIG_MULTI  = (3.55 * CM, 3.55 * CM)   # multilevel chart: 3.55 cm diameter

# Outer ring calibrated for exactly 0.5 cm physical height in a 3.55 cm figure.
#   physical radius  = 3.55/2 = 1.775 cm
#   target outer ring = 0.5 cm
#   → outer_h / (inner_top + gap + outer_h) = 0.5 / 1.775
#   → outer_h = 0.431,  ylim = 1.1 + 0.431 = 1.531
OUTER_BOTTOM = 1.1
OUTER_HEIGHT = 0.431
YLIM_MULTI   = OUTER_BOTTOM + OUTER_HEIGHT   # 1.531
YLIM_SIMPLE  = 1.02    # tiny overshoot so top bar edge isn't clipped

plt.rcParams["svg.fonttype"] = "none"   # SVG text as text, not paths

ACCESSIONS_FILE = Path("data/75k_unstable/ZF_SraEsearch-2025-06-22.csv")
SEQDETECTIVE_FILE = Path("data/75k_unstable/seq-detective-judgement-summary-augmented.txt")
ANNOTATIONS_FILE = Path("data/zf-core-v2-74K_problematic_with_bulk.kmers.allcols.csv")
BATLOW_PALETTE = Path("palette/batlow/DiscretePalettes/batlow100.txt")
OUTPUT_DIR = Path("figures")
OUTPUT_DIR.mkdir(exist_ok=True)

# =============================================================================
# Load batlow palette
# =============================================================================

def load_batlow_colors(n_colors=10):
    """Load colors from batlow palette file."""
    if BATLOW_PALETTE.exists():
        colors_raw = BATLOW_PALETTE.read_text().strip().split('\n')
        # Convert RGB values to hex
        colors = []
        for line in colors_raw:
            if line.strip() and not line.startswith('#') and not line.startswith('S'):
                parts = line.split()
                # Format: R G B batlow-X #RRGGBB
                # Use the hex value directly if available
                if '#' in line:
                    hex_color = [p for p in parts if p.startswith('#')][0]
                    colors.append(hex_color)
                elif len(parts) >= 3:
                    # Parse RGB (0-255 range)
                    try:
                        rgb = [int(parts[i]) for i in range(3)]
                        hex_color = '#{:02x}{:02x}{:02x}'.format(rgb[0], rgb[1], rgb[2])
                        colors.append(hex_color)
                    except (ValueError, IndexError):
                        continue
        # Sample evenly if we need fewer colors
        if len(colors) > 0:
            step = max(1, len(colors) // n_colors)
            return [colors[i * step] for i in range(min(n_colors, len(colors)))]

    # Fallback
    return [
        '#03185a', '#0B2F8A', '#1E4B99', '#3568A6', '#4C84B1',
        '#63A0BC', '#7ABCC7', '#92D8D2', '#AAF4DD', '#be9036'
    ]

# =============================================================================
# Load and process data
# =============================================================================

# Load accessions
accessions = pl.read_csv(ACCESSIONS_FILE)
print(f"Loaded {len(accessions)} accessions")

# Load seq-detective judgements
seqdetective = pl.read_csv(
    SEQDETECTIVE_FILE,
    separator="\t",
    has_header=False,
    new_columns=['id', 'file1', 'file2', 'judgement1', 'judgement2', 'reason']
).with_columns(
    accession=pl.col("id").str.strip_suffix("_subsample"),
    outcome=pl.concat_str(
        pl.col("judgement1"),
        pl.col("judgement2"),
        ignore_nulls=True
    )
)
print(f"Loaded {len(seqdetective)} seq-detective judgements")

# Load manual technology annotations
annotations = pl.read_csv(
    ANNOTATIONS_FILE,
    columns=['accession', 'resolution', 'group', 'technology', 'chemistry_ver', 'library_var']
)
print(f"Loaded {len(annotations)} manual annotations")

# Join all data
merged = (
    accessions
    .rename({"Run": "accession"})
    .join(seqdetective.select("accession", "outcome"), on="accession", how="left")
    .join(annotations, on="accession", how="left")
)
print(f"Merged dataset: {len(merged)} records")

# =============================================================================
# Categorize technologies with hierarchical fallback
# =============================================================================

# Simplify technology categories for visualization
def map_technology_category(df):
    """Map technology to simplified categories."""
    return (
        df.with_columns(
            tech_category=pl.when(pl.col("technology").str.contains("Smart-seq|SMART-seq"))
                .then(pl.lit("Smart-seq"))
                .when(pl.col("technology").str.contains("DeTCT"))
                .then(pl.lit("DeTCT"))
                .when(pl.col("technology").str.contains("10x"))
                .then(pl.lit("10x Genomics"))
                .when(pl.col("technology").str.contains("CEL-Seq"))
                .then(pl.lit("CEL-Seq"))
                .when(pl.col("technology").str.contains("CORALL"))
                .then(pl.lit("CORALL"))
                .when(pl.col("technology").str.contains("sci-RNA-seq"))
                .then(pl.lit("sci-RNA-seq"))
                .when(pl.col("technology").str.contains("inDrop"))
                .then(pl.lit("inDrop"))
                .when(pl.col("technology").str.contains("CLIP|RIP|ChIP"))
                .then(pl.lit("Protein-RNA"))
                .otherwise(pl.lit("Other"))
        )
    )

merged = map_technology_category(merged)

# Count technologies
tech_counts = (
    merged.filter(pl.col("outcome").is_not_null())
    .group_by("tech_category")
    .agg(count=pl.len())
    .sort("count", descending=True)
)
print(f"\nTechnology distribution:")
print(tech_counts)

# =============================================================================
# Shared color maps (computed once, consistent across both charts)
# =============================================================================

# Tech colors: assign from batlow across ALL categories (including Other), sorted by count.
# Both charts use this same map so colors are identical for every category.
_all_tech = tech_counts.to_pandas()
_batlow_all = load_batlow_colors(len(_all_tech))
TECH_COLOR_MAP = dict(zip(_all_tech['tech_category'], _batlow_all))

# Outcome colors: aligned with Fig1_B_2 seq-detective palette.
#   BB / TB / BT / TT mirror the B-B / T-B / B-T / T-T colors exactly.
#   B  / T  are single-end equivalents: close to BB/TT but visually distinct.
OUTCOME_COLORS = {
    'BB': '#99882c',  # both biological   — olive        (= B-B in Fig1_B_2)
    'TB': '#426f52',  # M1 tech, M2 bio  — forest green  (= T-B in Fig1_B_2)
    'BT': '#0b2c5c',  # M1 bio, M2 tech  — dark navy     (= B-T in Fig1_B_2)
    'TT': '#f29d6c',  # both technical   — coral/peach   (= T-T in Fig1_B_2)
    'B':  '#c4b030',  # SE biological    — lighter olive (close to BB, distinguishable)
    'T':  '#d4603c',  # SE technical     — darker coral  (close to TT, distinguishable)
}

# =============================================================================
# Create multi-level pie chart
# =============================================================================

def create_multilevel_pie(df, output_path, draw_outcomes=True, exclude_categories=None, label_r_overrides=None):
    """Create multi-level polar bar chart (pie chart).

    Parameters
    ----------
    exclude_categories : list[str] | None
        Technology categories to omit entirely from this chart (e.g. ['Other']).
        Their slice area is removed; proportions are recomputed over the remainder.
    """

    # Filter to samples with outcomes, then optionally drop excluded categories
    df_clean = df.filter(pl.col("outcome").is_not_null())
    if exclude_categories:
        df_clean = df_clean.filter(~pl.col("tech_category").is_in(exclude_categories))
    total_samples = len(df_clean)
    print(f"\nPlotting {total_samples} samples with outcomes")

    # Get technology counts
    tech_data = (
        df_clean
        .group_by("tech_category")
        .agg(count=pl.len())
        .sort("count", descending=True)
    ).to_pandas()

    # ── Figure: exact square so the inscribed polar circle matches the spec ──
    # Legend is placed manually in Illustrator; do not include it here.
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.sans-serif'] = ['Arial', 'Helvetica Neue', 'DejaVu Sans']
    fig_d = FIG_MULTI if draw_outcomes else FIG_SIMPLE
    fig, ax = plt.subplots(figsize=fig_d, subplot_kw=dict(projection='polar'))
    _p = 0
    fig.subplots_adjust(left=_p, right=1 - _p, top=1 - _p, bottom=_p)

    # ylim controls radial scale; outer ring spans OUTER_BOTTOM → YLIM_MULTI
    ylim_max = YLIM_MULTI if draw_outcomes else YLIM_SIMPLE

    # Calculate angular positions
    gap_fraction = 0.02
    n_categories = len(tech_data)
    total_gap = gap_fraction * 2 * np.pi
    gap_width = total_gap / max(n_categories, 1)
    available_space = 2 * np.pi - total_gap

    theta_start = np.pi + gap_width
    label_info = []

    # Draw each technology category
    for idx, row in tech_data.iterrows():
        tech_cat = row['tech_category']
        tech_count = row['count']
        tech_proportion = tech_count / total_samples
        tech_width = tech_proportion * available_space

        if tech_width < 0.01:
            continue

        tech_color = TECH_COLOR_MAP.get(tech_cat, '#808080')
        theta_center = theta_start + tech_width / 2

        # Inner ring — technology category (r 0 → 1.0)
        ax.bar(
            theta_center, 1.0, width=tech_width, bottom=0.0,
            color=tech_color, alpha=0.7, edgecolor='white', linewidth=1
        )

        if tech_width > 0.1 and tech_cat != 'CORALL':
            display_name = '10x' if tech_cat == '10x Genomics' else tech_cat
            label_info.append({
                'theta': theta_center,
                'text': f"{display_name}\n({tech_proportion:.1%})",
                'color': tech_color,
            })

        # Outer ring — processing outcomes
        if draw_outcomes:
            outcome_data = (
                df_clean
                .filter(pl.col("tech_category") == tech_cat)
                .group_by("outcome")
                .agg(count=pl.len())
                .sort("count", descending=True)
            ).to_pandas()

            outcome_start = theta_start
            for _, outcome_row in outcome_data.iterrows():
                outcome = outcome_row['outcome']
                outcome_count = outcome_row['count']
                outcome_proportion = outcome_count / tech_count
                outcome_width = tech_width * outcome_proportion

                if outcome_width > 0.005:
                    outcome_theta = outcome_start + outcome_width / 2
                    outcome_color = OUTCOME_COLORS.get(outcome, '#808080')

                    # Outer ring: OUTER_BOTTOM → OUTER_BOTTOM + OUTER_HEIGHT
                    # Physical height = OUTER_HEIGHT / ylim_max × (FIG_MULTI/2) ≈ 0.5 cm
                    ax.bar(
                        outcome_theta, OUTER_HEIGHT, width=outcome_width,
                        bottom=OUTER_BOTTOM,
                        color=outcome_color, alpha=0.9, edgecolor='white', linewidth=0.5
                    )

                    # Label large outcome segments — centred inside the outer ring
                    if outcome_width > 0.15:
                        ax.text(
                            outcome_theta, OUTER_BOTTOM + OUTER_HEIGHT / 2, outcome,
                            ha='center', va='center', fontsize=5,
                            color='white', fontweight='bold'
                        )

                outcome_start += outcome_width

        theta_start += tech_width + gap_width

    # Draw inner-ring labels inside the wedges (no leader lines needed;
    # fine-tuning and the legend are done manually in Illustrator)
    text_dist = 0.35
    last_theta = -1

    for label in label_info:
        theta = label['theta']

        # Alternate radial position when two labels are angularly close
        if abs(theta - last_theta) < 0.628:
            text_dist = -text_dist

        last_theta = theta
        text_r = 0.5 + text_dist   # 0.85 or 0.15 — both inside inner ring

        # Per-label override: caller can force specific labels to a fixed r position
        if label_r_overrides:
            display_name = label['text'].split('\n')[0]
            if display_name in label_r_overrides:
                text_r = label_r_overrides[display_name]

        ax.annotate(
            label['text'],
            xy=(theta, 0.5),
            xytext=(theta, text_r),
            xycoords='data',
            textcoords='data',
            arrowprops=dict(arrowstyle='-', color='black', lw=0.8, alpha=0.8),
            ha='center', va='center', fontweight='bold', fontsize=5,
            bbox=dict(boxstyle="round,pad=0.2", facecolor='white', alpha=0.9,
                      edgecolor='gray', linewidth=0.5)
        )

    # ── Axis cleanup ──
    ax.set_theta_zero_location('N')
    ax.set_theta_direction(-1)
    ax.set_ylim(0, ylim_max)
    ax.set_rticks([])
    ax.set_thetagrids([])
    ax.grid(False)

    # Save at exact figure size — no bbox_inches expansion, no tight_layout
    plt.savefig(output_path, dpi=300)
    print(f"\nSaved plot to {output_path}")
    plt.close()

# Panel B.1a — simple tech pie: all categories (including "Other"); no outer ring.
# DeTCT label pushed out to avoid overlapping 10x.
create_multilevel_pie(
    merged,
    OUTPUT_DIR / "Fig1_B_1_tech_pie.svg",
    draw_outcomes=False,
    exclude_categories=None,
    label_r_overrides={'DeTCT': 0.85},
)

# Panel B.1b — multilevel tech × outcome: "Other" dropped so non-bulk techs
# get detail; outer ring is Seq-Detective per-accession outcome.
# 10x label pushed out to avoid overlapping CEL-Seq.
create_multilevel_pie(
    merged,
    OUTPUT_DIR / "Fig1_B_1_seqdetective_pie.svg",
    draw_outcomes=True,
    exclude_categories=['Other'],
    label_r_overrides={'10x': 0.85},
)

print("\nDone!")
