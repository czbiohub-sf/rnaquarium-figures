#!/usr/bin/env python3
"""
Figure 1 Panel B.5: Pipeline CPU resource usage

Horizontal bar chart of total CPU-hours per pipeline step (completed + failed).
Steps ordered top-to-bottom by pipeline execution order.
"""

import polars as pl
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import csv
from pathlib import Path

# =============================================================================
# Configuration
# =============================================================================

TRACE_FILE = Path("data/75k_unstable/trace-merged-dangerously.txt")
OUTPUT_DIR = Path("figures")
OUTPUT_DIR.mkdir(exist_ok=True)

CM = 1 / 2.54
FIG_SIZE = (6.178 * CM, 6.225 * CM)
DPI = 300

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial', 'Helvetica Neue', 'DejaVu Sans']
plt.rcParams['svg.fonttype'] = 'none'   # SVG text as text, not paths
FONT_LABEL = 6
FONT_TICK  = 5
FONT_LEGEND = 5

# Pipeline step display order and labels
STEP_ORDER = [
    "download",
    "fastp",
    "star_counts",
    "feature_count",
    "kb_negative",
    "hisat2",
    "star",
    "bowtie2",
    "dedup",
    "gsnap",
]
STEP_LABELS = {
    "download":      "Download\n+ Seq-Detective",
    "fastp":         "fastp",
    "star_counts":   "STAR (counts)",
    "feature_count": "featureCounts",
    "kb_negative":   "Kallisto",
    "hisat2":        "HISAT2",
    "star":          "STAR",
    "bowtie2":       "Bowtie2",
    "dedup":         "Deduplication",
    "gsnap":         "GSNAP",
}

COLOR_COMPLETED = "#426f52"   # forest green (batlow)
COLOR_FAILED    = "#f29d6c"   # coral/peach  (batlow)

# =============================================================================
# Load trace
# =============================================================================

rows = []
with open(TRACE_FILE) as f:
    reader = csv.reader(f, delimiter="\t")
    header = next(reader)
    for row in reader:
        if row and row[0] != "task_id":
            rows.append(row)

df = pl.DataFrame(rows, schema=header, orient="row")
print(f"Trace rows loaded: {df.height:,}")

# =============================================================================
# Compute CPU-hours per process × status
# =============================================================================

active = (
    df.filter(pl.col("status").is_in(["COMPLETED", "FAILED"]))
    .filter((pl.col("realtime") != "-") & (pl.col("cpus") != "-"))
    .with_columns([
        pl.col("cpus").cast(pl.Float64),
        pl.col("realtime").cast(pl.Float64),   # milliseconds
    ])
    .with_columns(cpu_hours=pl.col("cpus") * pl.col("realtime") / 3_600_000)
)

cpu_by_step = (
    active.group_by(["process", "status"])
    .agg(pl.col("cpu_hours").sum())
    .pivot(on="status", index="process", values="cpu_hours")
    .rename({"COMPLETED": "completed_h", "FAILED": "failed_h"})
    .with_columns([
        pl.col("completed_h").fill_null(0.0),
        pl.col("failed_h").fill_null(0.0),
    ])
)

# Filter to steps in STEP_ORDER only
cpu_by_step = cpu_by_step.filter(pl.col("process").is_in(STEP_ORDER))

# Convert to k CPU-hours
cpu_by_step = cpu_by_step.with_columns([
    (pl.col("completed_h") / 1000).alias("completed_kh"),
    (pl.col("failed_h")    / 1000).alias("failed_kh"),
])

print("\nCPU-hours by step:")
print(cpu_by_step.sort("completed_h", descending=True))

# =============================================================================
# Plot
# =============================================================================

# Build ordered arrays
rows_d = {r["process"]: r for r in cpu_by_step.to_dicts()}
labels    = [STEP_LABELS[s] for s in STEP_ORDER if s in rows_d]
completed = [rows_d[s]["completed_kh"] for s in STEP_ORDER if s in rows_d]
failed    = [rows_d[s]["failed_kh"]    for s in STEP_ORDER if s in rows_d]
totals    = [c + f for c, f in zip(completed, failed)]
y         = list(range(len(labels)))


def _bar_label(ax, x_pos, y_pos, value, max_val):
    """Annotate bar end if value is large enough to warrant a label."""
    if max_val > 0 and value / max_val >= 0.03:
        ax.text(
            x_pos + max_val * 0.01, y_pos,
            f"{value:.1f}",
            va='center', ha='left', fontsize=FONT_TICK,
        )


def plot_cpu_chart(output_path, mode):
    """
    mode='total'     — single bar per step (completed + failed)
    mode='breakdown' — stacked completed / failed bars with legend
    """
    max_val = max(totals) if totals else 1.0
    fig, ax = plt.subplots(figsize=FIG_SIZE)

    if mode == 'total':
        ax.barh(y, totals, color=COLOR_COMPLETED, height=0.6, zorder=2)
        for i, val in enumerate(totals):
            _bar_label(ax, val, i, val, max_val)
    else:
        ax.barh(y, completed, color=COLOR_COMPLETED, height=0.6, zorder=2)
        ax.barh(y, failed, left=completed, color=COLOR_FAILED, height=0.6, zorder=2)
        for i, val in enumerate(totals):
            _bar_label(ax, val, i, val, max_val)
        legend = ax.legend(
            handles=[
                mpatches.Patch(color=COLOR_COMPLETED, label="Completed"),
                mpatches.Patch(color=COLOR_FAILED,    label="Failed"),
            ],
            fontsize=FONT_LEGEND,
            frameon=True, framealpha=0.7, edgecolor='black', fancybox=False,
            loc='lower right',
        )
        legend.get_frame().set_linewidth(0.5)

    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=FONT_TICK)
    ax.invert_yaxis()
    ax.set_xlabel("CPU-hours (×10³)", fontsize=FONT_LABEL)
    ax.tick_params(axis='x', labelsize=FONT_TICK)
    ax.grid(True, axis='x', alpha=0.15, linewidth=0.5)
    ax.set_axisbelow(True)
    for spine in ax.spines.values():
        spine.set_linewidth(0.5)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    fig.tight_layout(pad=0)
    plt.savefig(output_path, dpi=DPI)
    print(f"Saved to {output_path}")
    plt.close()


plot_cpu_chart(OUTPUT_DIR / "Fig1_B_5_resource_usage_total.svg",     mode='total')
plot_cpu_chart(OUTPUT_DIR / "Fig1_B_5_resource_usage_breakdown.svg", mode='breakdown')
