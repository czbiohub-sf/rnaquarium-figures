#!/usr/bin/env python3
"""
Plot GO enrichment summary from a CSV produced by the pipeline.

Reads {out_prefix}_summary.csv (must contain at least columns:
  - set
  - resolution
  - n_sig_terms_unique
Optionally:
  - n_clusters

Produces the same figure as the pipeline:
  - Left y-axis: unique significant GO terms vs. resolution (solid, filled circles)
  - Right y-axis: number of clusters (dashed, thinner, semi-transparent, open circles)
  - Y-axis labels are annotated with ● (left) and ○ (right)

Usage:
  python plot_enrichment_summary_from_csv.py --csv enrichment_summary.csv --out-png enrichment_summary_plot.png --out-prefix enrichment
"""

from __future__ import annotations

import argparse
import sys
from typing import Optional, List

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def _coerce_summary_df(df: pd.DataFrame) -> pd.DataFrame:
    required_cols = {"set", "resolution", "n_sig_terms_unique"}
    missing = required_cols.difference(df.columns)
    if missing:
        raise ValueError(f"CSV is missing required columns: {sorted(missing)}")

    # Ensure numeric types
    df = df.copy()
    df["resolution"] = pd.to_numeric(df["resolution"], errors="coerce")
    df["n_sig_terms_unique"] = pd.to_numeric(df["n_sig_terms_unique"], errors="coerce")
    if "n_clusters" in df.columns:
        df["n_clusters"] = pd.to_numeric(df["n_clusters"], errors="coerce")

    # Drop rows with NaN in critical columns
    df = df[np.isfinite(df["resolution"]) & np.isfinite(df["n_sig_terms_unique"])].copy()

    # Ensure 'set' is present and string-typed
    if "set" not in df.columns:
        df["set"] = "series"
    else:
        df["set"] = df["set"].astype(str)

    return df

def plot_summary_from_df(summary_df: pd.DataFrame, out_png: Optional[str] = None, out_prefix: Optional[str] = None):
    """
    Create the summary plot with dual y-axes.
    """
    df = _coerce_summary_df(summary_df)
    if df.empty:
        print("WARNING: empty summary_df; nothing to plot.", file=sys.stderr)
        return

    has_clusters = "n_clusters" in df.columns and df["n_clusters"].notna().any()
    labels: List[str] = df["set"].unique().tolist()

    fig, ax = plt.subplots(figsize=(8.0, 5.2))
    ax2 = ax.twinx() if has_clusters else None

    for label in labels:
        sub = df[df["set"] == label].sort_values("resolution")

        # Left axis: unique significant terms (solid line, filled-circle markers)
        line, = ax.plot(
            sub["resolution"].values,
            sub["n_sig_terms_unique"].values,
            marker="o",            # filled circle (closed)
            label=label,
        )

        # Right axis: number of clusters (dashed, thinner, semi-transparent, open-circle markers)
        if has_clusters:
            ax2.plot(
                sub["resolution"].values,
                sub["n_clusters"].values,
                linestyle="--",
                linewidth=1.0,            # thinner dashed line
                alpha=0.5,                # semi-transparent
                marker="o",
                markerfacecolor="none",   # open circle
                markeredgecolor=line.get_color(),
                markeredgewidth=1.0,
                color=line.get_color(),   # match left series color
            )

    ax.set_xlabel("Leiden resolution")

    # UPDATED axis labels:
    # Left: solid line through filled circle
    ax.set_ylabel("Unique significant GO terms (q ≤ threshold)\n—●—")
    # Right: dashed line through open circle (uses spaced en dashes to suggest dashes)
    if has_clusters and ax2 is not None:
        ax2.set_ylabel("Number of clusters\n– – ○ – – ")

    title_main = "GO enrichment vs. clustering resolution (g:Profiler)"
    title_suffix = f"\n{out_prefix}" if out_prefix else ""
    annotation = "\n(dashed = # clusters, right y-axis)" if has_clusters else ""
    ax.set_title(title_main + title_suffix + annotation)

    ax.grid(True, linestyle="--", alpha=0.4)
    ax.legend(title="Series", frameon=False, ncol=1)
    fig.tight_layout()

    if out_png:
        fig.savefig(out_png, dpi=150)
    else:
        plt.show()
    plt.close(fig)



def main():
    p = argparse.ArgumentParser(
        description="Plot enrichment summary from CSV (dual y-axes with clusters).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--csv", required=True, help="Path to {out_prefix}_summary.csv")
    p.add_argument("--out-png", default=None, help="Path to write PNG (if omitted, shows window)")
    p.add_argument("--out-prefix", default=None, help="Optional text appended on title's second line")
    args = p.parse_args()

    df = pd.read_csv(args.csv)
    plot_summary_from_df(df, out_png=args.out_png, out_prefix=args.out_prefix)
    print(f"[INFO] Wrote: {args.out_png}" if args.out_png else "[INFO] Displayed figure window")


if __name__ == "__main__":
    main()
