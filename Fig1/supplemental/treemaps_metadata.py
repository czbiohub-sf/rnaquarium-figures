#!/usr/bin/env python3
"""
Supplemental: Metadata summary treemaps

Treemaps of developmental-stage and tissue-category distributions.
Previously Panel C of Figure 1; moved to supplemental when the 2026-04 layout
replaced Panel C with the pipeline-outcome tables.

Run from the Fig1/ root.
"""

import textwrap
import polars as pl
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.transforms as mtransforms
import squarify
from pathlib import Path

# =============================================================================
# Configuration
# =============================================================================

METADATA_FILE = Path("/hpc/projects/balla_group/sra_experiments/SRA_metadata/dec2025_75k_submitteradded/all_zf_dates_devstage_tissue_tech_curated.tsv")
BATLOW_PALETTE = Path("palette/batlow/DiscretePalettes/batlow100.txt")
OUTPUT_DIR = Path("figures/supplemental")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

CM = 1 / 2.54
PANEL_W = 6.25 * CM    # panel C2: 6.25 × 5 cm
PANEL_H = 5.0 * CM

mpl.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['font.sans-serif'] = ['Arial', 'Helvetica Neue', 'DejaVu Sans']
mpl.rcParams['svg.fonttype'] = 'none'   # SVG text as text, not paths

# Gray overrides for uninformative categories (consistent with Panel D UMAPs)
GRAY_CATS = {
    "Undetermined": "#9E9E9E",
    "All anatomical structures": "#BDBDBD",
}

# =============================================================================
# Load batlow palette
# =============================================================================

def load_batlow_colors(n_colors):
    if BATLOW_PALETTE.exists():
        colors = []
        for line in BATLOW_PALETTE.read_text().strip().split('\n'):
            if not line.strip() or line.startswith('#') or line.startswith('S'):
                continue
            parts = line.split()
            if '#' in line:
                colors.append([p for p in parts if p.startswith('#')][0])
            elif len(parts) >= 3:
                try:
                    colors.append('#{:02x}{:02x}{:02x}'.format(*[int(parts[i]) for i in range(3)]))
                except (ValueError, IndexError):
                    continue
        if colors:
            # Golden-ratio stride: adjacent entries are maximally far apart in palette,
            # so nearby legend items are visually distinct even with many categories.
            # gcd(stride, len) == 1 guarantees no collisions for n_colors < len(colors).
            stride = int(len(colors) * 0.618)
            return [colors[(i * stride) % len(colors)] for i in range(n_colors)]

    return [
        '#03185a', '#1E4B99', '#4C84B1', '#7ABCC7', '#AAF4DD',
        '#be9036', '#c4741e', '#c43e1e', '#9e1e1e', '#5a0a0a',
        '#3a0a0a', '#1a0a0a', '#0a0a0a', '#0B2F8A', '#3568A6',
        '#63A0BC', '#92D8D2', '#c45a1e', '#7a1e1e', '#050505', '#000000'
    ]

def _is_dark(hex_color):
    r, g, b = int(hex_color[1:3], 16), int(hex_color[3:5], 16), int(hex_color[5:7], 16)
    return 0.299 * r + 0.587 * g + 0.114 * b < 128

def _draw_centered_lines(ax, fig, cx, cy, lines, fontsize, linespacing, text_color):
    # Draw each line as a separate text call to avoid multi-line SVG <g> anchor issues.
    # Vertical offsets are in points via ScaledTranslation so they're font-metric-independent.
    n = len(lines)
    for i, line in enumerate(lines):
        y_offset_in = ((n - 1) / 2 - i) * fontsize * linespacing / 72
        trans = ax.transData + mtransforms.ScaledTranslation(0, y_offset_in, fig.dpi_scale_trans)
        ax.text(cx, cy, line, ha='center', va='center',
                fontsize=fontsize, color=text_color, fontweight='bold',
                transform=trans)

# =============================================================================
# Load and process data
# =============================================================================

metadata = pl.read_csv(
    METADATA_FILE,
    separator="\t",
    columns=["run.accession", "devstage_curation_coarse", "tissue_curation_coarse"]
)
print(f"Loaded {len(metadata)} metadata records")

devstage_counts = (
    metadata
    .filter(pl.col("devstage_curation_coarse").is_not_null())
    .group_by("devstage_curation_coarse")
    .agg(count=pl.len())
    .sort("count", descending=True)
)

tissue_counts = (
    metadata
    .filter(pl.col("tissue_curation_coarse").is_not_null())
    .group_by("tissue_curation_coarse")
    .agg(count=pl.len())
    .sort("count", descending=True)
)

# =============================================================================
# Treemap
# =============================================================================

def _wrap_label(text, box_width_units, norm_w, fig_width_inches, fontsize, safety=0.90):
    # coordinate space is 0-norm_w; estimate chars that fit in box_width
    pts_per_unit = fig_width_inches * 72 / norm_w
    chars_per_unit = pts_per_unit / (0.60 * fontsize)
    max_chars = max(6, int(box_width_units * chars_per_unit * safety))
    return textwrap.fill(text, width=max_chars)


def create_treemap(data_df, category_col, output_path, title, figsize=(PANEL_W, PANEL_H), legend=False):
    data = data_df.to_pandas()
    total = data['count'].sum()
    labels = data[category_col].tolist()   # count-descending order
    sizes = data['count'].tolist()
    # Assign batlow to non-gray categories in count-descending order (consistent with Panel D)
    non_gray = [l for l in labels if l not in GRAY_CATS]
    batlow = load_batlow_colors(len(non_gray))
    color_map = {label: batlow[i] for i, label in enumerate(non_gray)}
    color_map.update(GRAY_CATS)
    colors = [color_map[l] for l in labels]

    fig, ax = plt.subplots(figsize=figsize)

    # Use figure aspect ratio so squarify fills the full panel without whitespace
    norm_w = 100 * figsize[0] / figsize[1]
    norm_h = 100
    normed = squarify.normalize_sizes(sizes, norm_w, norm_h)
    rects = squarify.squarify(normed, 0, 0, norm_w, norm_h)

    for rect, label, size, color in zip(rects, labels, sizes, colors):
        x, y, w, h = rect['x'], rect['y'], rect['dx'], rect['dy']
        ax.add_patch(mpatches.Rectangle(
            (x, y), w, h,
            facecolor=color, edgecolor='white', linewidth=1.5
        ))

        pct = size / total
        text_color = 'white' if _is_dark(color) else '#1a1a1a'

        if pct >= 0.10:
            lines = _wrap_label(label, w, norm_w, figsize[0], 7).split('\n') + [f"{pct:.1%}", f"(n={size:,})"]
            _draw_centered_lines(ax, fig, x + w / 2, y + h / 2, lines, 7, 1.4, text_color)
        elif pct >= 0.07:
            lines = _wrap_label(label, w, norm_w, figsize[0], 6).split('\n') + [f"{pct:.1%}"]
            _draw_centered_lines(ax, fig, x + w / 2, y + h / 2, lines, 6, 1.3, text_color)
        elif pct >= 0.03:
            lines = _wrap_label(label, w, norm_w, figsize[0], 5).split('\n') + [f"{pct:.1%}"]
            _draw_centered_lines(ax, fig, x + w / 2, y + h / 2, lines, 5, 1.3, text_color)

    if legend:
        legend_patches = [
            mpatches.Patch(facecolor=colors[i], label=f"{labels[i]} (n={sizes[i]:,})")
            for i in range(len(labels))
        ]
        ax.legend(
            handles=legend_patches,
            loc='upper left', bbox_to_anchor=(1.01, 1.0),
            fontsize=5, framealpha=0.9, borderpad=0.4,
            handlelength=1.0, handleheight=0.8
        )

    ax.set_xlim(0, norm_w)
    ax.set_ylim(0, norm_h)
    ax.axis('off')
    if title:
        ax.set_title(title, fontsize=7, pad=6)

    fig.tight_layout(pad=0)
    plt.savefig(output_path, dpi=300)
    print(f"Saved {output_path}")
    plt.close()

# =============================================================================
# Create plots
# =============================================================================

create_treemap(
    devstage_counts,
    "devstage_curation_coarse",
    OUTPUT_DIR / "Fig1_treemap_devstage.svg",
    None,
)
create_treemap(
    devstage_counts,
    "devstage_curation_coarse",
    OUTPUT_DIR / "Fig1_treemap_devstage_legend.svg",
    "Developmental Stage Distribution",
    legend=True,
)

create_treemap(
    tissue_counts,
    "tissue_curation_coarse",
    OUTPUT_DIR / "Fig1_treemap_tissue.svg",
    None,
)
create_treemap(
    tissue_counts,
    "tissue_curation_coarse",
    OUTPUT_DIR / "Fig1_treemap_tissue_legend.svg",
    "Tissue Type Distribution",
    legend=True,
)

print("\nDone!")
