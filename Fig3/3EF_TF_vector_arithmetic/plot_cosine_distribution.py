#!/usr/bin/env python3
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def parse_cols(arg_vals):
    if not arg_vals:
        return []
    cols = []
    for tok in arg_vals:
        cols.extend([t for t in str(tok).replace(",", " ").split() if t])
    # preserve order, dedupe
    out, seen = [], set()
    for c in cols:
        if c not in seen:
            out.append(c); seen.add(c)
    return out

def parse_marks(mark_args):
    """
    Parse --mark entries of the form:
      --mark 0.42:SOX2  --mark 0.10:'random baseline'
    Returns list of tuples: [(x_float, "label"), ...]
    """
    marks = []
    if not mark_args:
        return marks
    for m in mark_args:
        s = str(m)
        if ":" in s:
            x_str, label = s.split(":", 1)
            x = float(x_str)
            label = label.strip().strip("'").strip('"')
        else:
            # if only a number is provided, use it as the label too
            x = float(s)
            label = s
        marks.append((x, label))
    return marks

def main():
    ap = argparse.ArgumentParser(
        description="Plot distribution(s) from cosine_to_targets.csv with vertical marks and labels."
    )
    ap.add_argument("csv", help="Path to *_cosine_to_targets.csv produced by your earlier script.")
    ap.add_argument("--columns", "-c", nargs="*", default=None,
                    help="Target gene column(s) to plot (space/comma-separated). "
                         "If omitted, uses the first numeric column.")
    ap.add_argument("--bins", type=int, default=30, help="Number of histogram bins (default: 30).")
    ap.add_argument("--density", action="store_true",
                    help="Normalize histogram(s) to probability density.")
    ap.add_argument("--alpha", type=float, default=0.6, help="Bar transparency (default: 0.6).")
    ap.add_argument("--figsize", type=float, nargs=2, default=(8, 5),
                    help="Figure size in inches, e.g. --figsize 8 5")
    ap.add_argument("--title", type=str, default=None, help="Custom plot title.")
    ap.add_argument("--xlabel", type=str, default="Cosine similarity", help="X-axis label.")
    ap.add_argument("--ylabel", type=str, default=None, help="Y-axis label (auto if omitted).")
    ap.add_argument("--xlim", type=float, nargs=2, default=None, help="X limits, e.g. --xlim -0.2 1.0")
    ap.add_argument("--ylim", type=float, nargs=2, default=None, help="Y limits, e.g. --ylim 0 20")
    ap.add_argument("--mark", action="append", default=[],
                    help="Add a vertical line with label at x: 'x:Label'. "
                         "Repeat for multiple marks. Example: --mark 0.42:SOX2 --mark 0.0:baseline")
    ap.add_argument("--mark-y", type=float, default=None,
                    help="Y position for all mark labels (in data units). "
                         "If omitted, labels are placed near the top of the plot.")
    ap.add_argument("--mark-rotation", type=float, default=90.0,
                    help="Text rotation for mark labels (default: 90).")
    ap.add_argument("--out", type=str, default=None,
                    help="Output image path (.png, .pdf, etc). Defaults to alongside the CSV.")
    ap.add_argument("--show", action="store_true", help="Display the plot interactively.")

    args = ap.parse_args()

    csv_path = Path(args.csv)
    if not csv_path.exists():
        raise SystemExit(f"File not found: {csv_path}")

    # Load CSV (handle saved index column if present)
    df = pd.read_csv(csv_path)
    if df.columns[0].startswith("Unnamed"):
        df = df.set_index(df.columns[0])

    # Determine numeric columns
    num_cols = df.select_dtypes(include="number").columns.tolist()
    if len(num_cols) == 0:
        raise SystemExit("No numeric columns found in CSV; nothing to plot.")

    # Columns to plot
    requested = parse_cols(args.columns) if args.columns is not None else []
    if not requested:
        plot_cols = [num_cols[0]]
        print(f"No --columns provided; using first numeric column: {plot_cols[0]}")
    else:
        missing = [c for c in requested if c not in df.columns]
        if missing:
            raise SystemExit(f"Requested column(s) not found in CSV: {missing}")
        plot_cols = requested

    # Build data series list (drop NaNs)
    series = [(c, df[c].dropna().to_numpy()) for c in plot_cols]

    # Plot
    plt.figure(figsize=tuple(args.figsize))
    for name, arr in series:
        if arr.size == 0:
            print(f"Warning: column '{name}' has no data after dropping NaNs; skipping.")
            continue
        plt.hist(arr, bins=args.bins, density=args.density, alpha=args.alpha, label=name, edgecolor="black")

    # Labels & limits
    plt.xlabel(args.xlabel)
    if args.ylabel:
        plt.ylabel(args.ylabel)
    else:
        plt.ylabel("Density" if args.density else "Count")

    title = args.title
    if title is None:
        if len(series) == 1:
            title = f"Distribution of cosine similarity: {series[0][0]}"
        else:
            title = "Distribution of cosine similarity"
    plt.title(title)

    if args.xlim:
        plt.xlim(args.xlim)
    if args.ylim:
        plt.ylim(args.ylim)

    # Vertical lines with labels
    marks = parse_marks(args.mark)
    if marks:
        # For label y position
        ymin, ymax = plt.ylim()
        label_y = args.mark_y if args.mark_y is not None else (ymax * 0.95)
        for x, lab in marks:
            plt.axvline(x, linestyle="--", linewidth=1.5)
            # Keep labels within axes if possible
            if args.xlim:
                xmin, xmax = args.xlim
            else:
                xmin, xmax = plt.xlim()
            x_clamped = min(max(x, xmin), xmax)
            plt.text(
                x_clamped, label_y, lab,
                rotation=args.mark_rotation, ha="center", va="top",
                fontsize=10, bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="none", alpha=0.7)
            )

    if len(series) > 1:
        plt.legend(frameon=False)

    # Output
    if args.out is None:
        suffix = "_".join(plot_cols[:3]) + ("_and_more" if len(plot_cols) > 3 else "")
        out_path = csv_path.with_name(csv_path.stem + f"__dist_{suffix}.png")
    else:
        out_path = Path(args.out)
    plt.tight_layout()
    plt.savefig(out_path, dpi=200)
    print(f"Saved plot to: {out_path}")

    if args.show:
        plt.show()

if __name__ == "__main__":
    main()
