#!/usr/bin/env python3
"""
TMM-normalize a counts matrix, apply Duo's explicit filtering rules, and (by default)
Z-score standardize post-TMM for downstream analyses.

Key changes vs. previous version:
- Replaced CPM-based gene filtering with the requested stepwise filters:
  1) Drop rows with all zeros
  2) Drop columns with all zeros
  3) Drop rows where ALL counts < min_any_count (default 10)
  4) Drop columns where ALL counts < min_any_count (default 10)
  5) Drop rows/columns with > zero_fraction_threshold (default 0.90) zeros
  6) Report NaN/Inf counts
- Output scaling defaults to Z-score standardization (per gene across samples) AFTER TMM.
  You can switch to CPM/log2-CPM or raw TMM counts with --output-scale.

Assumes input is (auto-detected) samples x genes after orientation step.
"""

import argparse
import os
import time
from typing import Literal, Optional

import numpy as np
import pandas as pd

# ---------------------------- I/O helpers ---------------------------------- #

def read_parquet(path: str, index_col: Optional[str] = None) -> pd.DataFrame:
    """Read a Parquet file and set an index, with verbose prints."""
    print(f"Reading parquet: {path}", flush=True)
    df = pd.read_parquet(path)
    print(f"Raw table shape: {df.shape}", flush=True)

    if index_col is not None and index_col in df.columns:
        df = df.set_index(index_col)
        print(f"Set index to column '{index_col}'.", flush=True)
    elif df.index.name is None:
        df.set_index(df.columns[0], inplace=True)
        print(f"No explicit index; used first column as index: '{df.index.name}'.", flush=True)

    print(f"After indexing, shape: {df.shape}", flush=True)
    try:
        print("Head (first 5 rows):", flush=True)
        print(df.head(), flush=True)
    except Exception:
        pass
    return df


def save_output(df: pd.DataFrame, output_path: str) -> None:
    outdir = os.path.dirname(output_path) or "."
    os.makedirs(outdir, exist_ok=True)
    if output_path.lower().endswith(".parquet"):
        df.to_parquet(output_path, index=True)
        print(f"Saved parquet: {output_path}", flush=True)
    else:
        compression = "gzip" if output_path.lower().endswith(".gz") else None
        df.to_csv(output_path, index=True, compression=compression)
        print(f"Saved CSV{' (gzip)' if compression else ''}: {output_path}", flush=True)


# ---------------------------- data wrangling -------------------------------- #

def coerce_orientation(
    df: pd.DataFrame,
    orientation: Literal["auto", "genes_by_cells", "cells_by_genes", "as_is"] = "auto",
) -> pd.DataFrame:
    """Ensure rows=samples/cells and columns=genes for downstream steps."""
    before = df.shape
    if orientation == "as_is":
        print(f"Orientation kept as-is (shape {before}).", flush=True)
        return df
    if orientation == "genes_by_cells":
        df = df.T
        print(f"Orientation 'genes_by_cells' -> transposed to {df.shape}", flush=True)
        return df
    if orientation == "cells_by_genes":
        print(f"Orientation 'cells_by_genes' (shape {before}).", flush=True)
        return df
    # auto heuristic: transpose if rows < cols (typical genes x samples input)
    if df.shape[0] < df.shape[1]:
        df = df.T
        print(f"Orientation 'auto': transposed from {before} to {df.shape}", flush=True)
    else:
        print(f"Orientation 'auto': left as {before}", flush=True)
    return df


def drop_trailing_non_gene_columns(df: pd.DataFrame, n: int) -> pd.DataFrame:
    if n <= 0:
        print("No trailing columns dropped.", flush=True)
        return df
    n_drop = min(n, df.shape[1])
    out = df.iloc[:, : df.shape[1] - n_drop]
    print(f"Dropped last {n_drop} columns: shape {df.shape} -> {out.shape}", flush=True)
    return out


def to_numeric_matrix(df: pd.DataFrame) -> pd.DataFrame:
    print("Coercing to numeric (non-convertible -> NaN -> 0).", flush=True)
    before_cols = df.shape[1]
    df_num = df.apply(pd.to_numeric, errors="coerce")
    dropped_cols = [c for c in df_num.columns if not df_num[c].notna().any()]
    if dropped_cols:
        print(f"Dropping {len(dropped_cols)} fully-NA columns after coercion.", flush=True)
        df_num = df_num.drop(columns=dropped_cols)
    na_ct = int(df_num.isna().sum().sum())
    if na_ct:
        print(f"Filled {na_ct} NA values with 0.", flush=True)
    df_num = df_num.fillna(0)
    print(f"Numeric columns: {before_cols} -> {df_num.shape[1]}", flush=True)
    return df_num


# ------------------------------ gene mapping -------------------------------- #

def map_genes_to_ensembl(
    df_cells_by_genes: pd.DataFrame,
    map_csv: Optional[str] = None,
    gene_col: str = "gene_id",
    ensembl_col: str = "Ensembl_gene_id",
) -> pd.DataFrame:
    """Map column names to Ensembl IDs and sum duplicate mappings."""
    if map_csv is None:
        print("No gene ID mapping provided; skipping.", flush=True)
        return df_cells_by_genes

    print(f"Mapping genes using: {map_csv}", flush=True)
    m = pd.read_csv(map_csv)
    if gene_col not in m.columns or ensembl_col not in m.columns:
        raise ValueError(
            f"Mapping file must contain columns '{gene_col}' and '{ensembl_col}'."
        )

    before_cols = df_cells_by_genes.shape[1]
    mapper = m.set_index(gene_col)[ensembl_col].to_dict()

    df = df_cells_by_genes.copy()
    df.columns = df.columns.map(lambda x: mapper.get(x, x))

    after_map_cols = df.shape[1]
    df = df.loc[:, df.columns.notna()]
    after_drop_na_cols = df.shape[1]

    # Sum duplicate Ensembl IDs
    dup_ct = after_drop_na_cols - len(pd.unique(df.columns))
    df = df.T.groupby(level=0).sum().T
    after_collapse_cols = df.shape[1]

    print(
        "Mapping summary: "
        f"start={before_cols}, after_map={after_map_cols}, "
        f"after_drop_unmapped={after_drop_na_cols}, collapsed_duplicates={dup_ct}, "
        f"final={after_collapse_cols}",
        flush=True,
    )
    return df


# ------------------------------- normalization ------------------------------ #

def print_cpm_row_sums_summary(df: pd.DataFrame, label: str = "") -> None:
    """Summarize per-row sums; for CPM these should be ~1e6."""
    rs = df.sum(axis=1).astype(float)
    q = rs.quantile([0.0, 0.25, 0.5, 0.75, 1.0])
    dev = (rs - 1e6).abs()
    pct_1 = (dev <= 1e4).mean() * 100.0   # within ±1%
    pct_2 = (dev <= 2e4).mean() * 100.0   # within ±2%
    print(
        f"{label}Row sums (CPM expected ≈1e6): "
        f"min={q.iloc[0]:.1f}, 25%={q.iloc[1]:.1f}, median={q.iloc[2]:.1f}, "
        f"75%={q.iloc[3]:.1f}, max={q.iloc[4]:.1f}",
        flush=True,
    )
    print(
        f"{label}Mean |Δ| from 1e6: {dev.mean():.1f}; "
        f"% within ±1%: {pct_1:.1f}%, ±2%: {pct_2:.1f}%",
        flush=True,
    )

def compute_cpm(df_counts: pd.DataFrame) -> pd.DataFrame:
    """Counts per million using row-wise library sizes."""
    libsizes = df_counts.sum(axis=1).astype(float)
    libsizes = libsizes.replace(0, np.nan)
    q = libsizes.quantile([0.0, 0.25, 0.5, 0.75, 1.0])
    print(
        "Library sizes (min,25%,50%,75%,max): "
        f"{q.iloc[0]:.1f}, {q.iloc[1]:.1f}, {q.iloc[2]:.1f}, {q.iloc[3]:.1f}, {q.iloc[4]:.1f}",
        flush=True,
    )
    cpm = df_counts.div(libsizes, axis=0) * 1e6
    return cpm.fillna(0.0)


def compute_logcpm(df_counts: pd.DataFrame, prior_count: float = 0.25) -> pd.DataFrame:
    """edgeR-style log2-CPM: log2(CPM + prior_count)."""
    print(f"Computing log2-CPM with prior_count={prior_count}.", flush=True)
    cpm = compute_cpm(df_counts)
    return np.log2(cpm + prior_count)


def tmm_normalize(df_cells_by_genes: pd.DataFrame) -> pd.DataFrame:
    """TMM normalization using rnanorm; returns a pandas DataFrame.
       rnanorm.TMM().fit_transform(...) returns CPM already, computed with the adjusted (TMM) library sizes. In other words, this function currently outputs TMM-CPM.
    """
    try:
        from rnanorm import TMM
    except ImportError as e:
        raise RuntimeError(
            "rnanorm is required for TMM normalization. Install with: pip install rnanorm"
        ) from e

    print("Normalizing with TMM...", flush=True)
    start = time.time()
    norm = TMM().set_output(transform="pandas").fit_transform(df_cells_by_genes)
    print(f"TMM normalization took {time.time() - start:.2f} s", flush=True)

    # ⇩ Sanity check: row sums should be ~1e6 for CPM
    print_cpm_row_sums_summary(norm, label="After TMM: ")
    return norm


def zscore_standardize(df: pd.DataFrame, by: Literal["gene", "sample"] = "gene") -> pd.DataFrame:
    """Z-score standardize.
    by="gene": standardize each gene (column) across samples (rows). [default]
    by="sample": standardize each sample (row) across genes (columns).
    Zero-variance vectors become all zeros.
    """
    axis = 0 if by == "gene" else 1
    print(f"Z-score standardizing by {by} (axis={axis}).", flush=True)
    if by == "gene":
        mean = df.mean(axis=0)
        std = df.std(axis=0, ddof=0).replace(0, np.nan)
        z = (df - mean) / std
    else:
        mean = df.mean(axis=1)
        std = df.std(axis=1, ddof=0).replace(0, np.nan)
        z = (df.sub(mean, axis=0)).div(std, axis=0)
    z = z.fillna(0.0)
    return z


# ------------------------------ filtering ---------------------------------- #

def filter_counts_explicit(
    df: pd.DataFrame,
    min_any_count: int = 10,
    zero_fraction_threshold: float = 0.90,
) -> pd.DataFrame:
    """Apply the requested filtering steps with verbose prints.

    Steps (on a samples x genes matrix):
    1) Remove rows with all zeros
    2) Remove columns with all zeros
    3) Remove rows where ALL counts < min_any_count
    4) Remove columns where ALL counts < min_any_count
    5) Remove rows and columns with > zero_fraction_threshold zeros
    6) Report NaN/Inf
    """
    print("Filtering count table...", flush=True)

    # 1) rows with all zeros
    print("removing rows with all zeros...", flush=True)
    before = df.shape
    df = df.loc[(df != 0).any(axis=1), :]
    print(f"{before} -> {df.shape}", flush=True)

    # 2) cols with all zeros
    print("removing columns with all zeros...", flush=True)
    before = df.shape
    df = df.loc[:, (df != 0).any(axis=0)]
    print(f"{before} -> {df.shape}", flush=True)

    # 3) rows where all counts < min_any_count
    print("removing rows where all counts < {min_any_count}...".format(min_any_count=min_any_count), flush=True)
    before = df.shape
    df = df.loc[(df >= min_any_count).any(axis=1), :]
    print(f"{before} -> {df.shape}", flush=True)

    # 4) cols where all counts < min_any_count
    print("removing columns where all counts < {min_any_count}...".format(min_any_count=min_any_count), flush=True)
    before = df.shape
    df = df.loc[:, (df >= min_any_count).any(axis=0)]
    print(f"{before} -> {df.shape}", flush=True)

    # 5) rows/cols with > threshold zeros
    print("removing rows and columns in which over {pct}% of the values are zeros...".format(
        pct=int(100 * zero_fraction_threshold)
    ), flush=True)

    row_zero_counts = (df == 0).sum(axis=1)
    col_zero_counts = (df == 0).sum(axis=0)

    row_limit = zero_fraction_threshold * df.shape[1]
    col_limit = zero_fraction_threshold * df.shape[0]

    print(
        f"row count where zeros > {zero_fraction_threshold*100:.0f}% (>{row_limit:.1f} of {df.shape[1]}): "
        f"{(row_zero_counts > row_limit).sum()}",
        flush=True,
    )
    print(
        f"column count where zeros > {zero_fraction_threshold*100:.0f}% (>{col_limit:.1f} of {df.shape[0]}): "
        f"{(col_zero_counts > col_limit).sum()}",
        flush=True,
    )

    rows_to_keep = row_zero_counts <= row_limit
    cols_to_keep = col_zero_counts <= col_limit

    before = df.shape
    df = df.loc[rows_to_keep, cols_to_keep]
    print(f"{before} -> {df.shape}", flush=True)

    # 6) sanity: NaN/Inf
    nan_ct = int(df.isna().sum().sum())
    inf_ct = int(df.isin([np.inf, -np.inf]).sum().sum())
    print(f"Number of nan in the count table: {nan_ct}")
    print(f"Number of inf in the count table: {inf_ct}")

    if df.shape[0] == 0 or df.shape[1] == 0:
        raise ValueError("All data filtered out; adjust thresholds.")

    return df


# ------------------------------- main pipeline ------------------------------ #

def process_counts_table(
    input_parquet_path: str,
    output_path: str,
    orientation: Literal["auto", "genes_by_cells", "cells_by_genes", "as_is"] = "auto",
    index_col: Optional[str] = None,
    drop_last_n_cols: int = 0,
    gene_id_map_path: Optional[str] = None,
    # filtering parameters
    min_any_count: int = 10,
    zero_fraction_threshold: float = 0.90,
    # output scaling and normalization
    output_scale: Literal["zscore", "tmm_counts", "cpm", "logcpm"] = "zscore",
    zscore_by: Literal["gene", "sample"] = "gene",
    prior_count: float = 1,
    do_tmm: bool = True,
) -> None:
    print("=== PIPELINE START ===", flush=True)
    df = read_parquet(input_parquet_path, index_col=index_col)

    print("-- Orientation & typing --", flush=True)
    df = coerce_orientation(df, orientation=orientation)
    df = drop_trailing_non_gene_columns(df, n=drop_last_n_cols)
    #df = to_numeric_matrix(df) # skip numeric conversion

    print("-- Gene ID mapping --", flush=True)
    df = map_genes_to_ensembl(df, gene_id_map_path)

    print("-- Requested explicit filtering --", flush=True)
    df = filter_counts_explicit(
        df,
        min_any_count=min_any_count,
        zero_fraction_threshold=zero_fraction_threshold,
    )
    print(f"Post-filter shape (samples x genes): {df.shape}", flush=True)

    print("-- Pre-normalization NA/Inf check --", flush=True)
    print(f"NaN count: {int(df.isna().sum().sum())} | Inf count: {int(df.isin([np.inf, -np.inf]).sum().sum())}", flush=True)

    if do_tmm:
        print("-- TMM normalization --", flush=True)
        df = tmm_normalize(df)
        print("-- Post-TMM NA/Inf check --", flush=True)
        print(f"NaN count: {int(df.isna().sum().sum())} | Inf count: {int(df.isin([np.inf, -np.inf]).sum().sum())}", flush=True)
    else:
        print("Skipping TMM normalization (per --no-tmm).", flush=True)

    print("-- Output scaling --", flush=True)
    if output_scale == "tmm_counts":
        out = df
        print("Output scale: TMM counts", flush=True)
    elif output_scale == "cpm":
        out = compute_cpm(df)
        print("Output scale: CPM", flush=True)
    elif output_scale == "logcpm":
        #out = compute_logcpm(df, prior_count=prior_count) # deprecated because df is already TMM-CPM

        #df is alreadyTMM-CPM because rnanorm.TMM returns CPM using adjusted lib sizes
        out = np.log2(df + prior_count)   # prior_count in CPM units (e.g. 0.25–1)
        print("Output scale: log2(TMM-CPM + prior_count)", flush=True)
    elif output_scale == "zscore":
        out = zscore_standardize(df, by=zscore_by)
        print(f"Output scale: Z-score by {zscore_by}", flush=True)
    else:
        raise ValueError("Invalid output_scale")

    print("-- Output integrity check --", flush=True)
    print(
        f"NaN count in output: {int(out.isna().sum().sum())} | Inf count: {int(out.isin([np.inf, -np.inf]).sum().sum())}",
        flush=True,
    )
    print(f"Final output shape: {out.shape}", flush=True)

    print(f"Saving to: {output_path}", flush=True)
    save_output(out, output_path)
    print("=== PIPELINE END ===", flush=True)


# --------------------------------- CLI -------------------------------------- #

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=(
            "TMM-normalize counts, apply explicit filters, and emit Z-score (default) or alternative scales."
        )
    )
    p.add_argument("input_parquet", help="Path to input .parquet")
    p.add_argument("output_path", help="Path to output (.parquet or .csv[.gz])")
    p.add_argument(
        "--orientation",
        choices=["auto", "genes_by_cells", "cells_by_genes", "as_is"],
        default="auto",
        help="How to interpret input orientation (default: auto)",
    )
    p.add_argument("--index-col", default=None, help="Optional column to set as index after loading")
    p.add_argument(
        "--drop-last-n-cols",
        type=int,
        default=0,
        help="Drop the last N columns if you know they are non-gene annotations (default: 0)",
    )
    p.add_argument(
        "--gene-id-to-ensembl-id-file",
        type=str,
        default=None,
        help="Optional CSV mapping with columns 'gene_id' and 'Ensembl_gene_id'",
    )
    # filtering params
    p.add_argument(
        "--min-any-count",
        type=int,
        default=10,
        help="Threshold for steps 3-4: keep rows/cols where ANY value >= this (default: 10)",
    )
    p.add_argument(
        "--zero-fraction-threshold",
        type=float,
        default=0.90,
        help="For step 5: drop rows/cols with > this fraction zeros (default: 0.90)",
    )
    # scaling / normalization
    p.add_argument(
        "--output-scale",
        choices=["zscore", "tmm_counts", "cpm", "logcpm"],
        default="tmm_counts",
        help="What to write: Z-score, TMM counts(default), CPM, or log2-CPM",
    )
    p.add_argument(
        "--zscore-by",
        choices=["gene", "sample"],
        default="gene",
        help="Z-score across samples for each gene (default) or across genes for each sample",
    )
    p.add_argument("--prior-count", type=float, default=1, help="Prior count added in log2-CPM")
    p.add_argument("--no-tmm", action="store_true", help="Skip TMM normalization")
    return p.parse_args()


if __name__ == "__main__":
    args = parse_args()
    process_counts_table(
        input_parquet_path=args.input_parquet,
        output_path=args.output_path,
        orientation=args.orientation,
        index_col=args.index_col,
        drop_last_n_cols=args.drop_last_n_cols,
        gene_id_map_path=args.gene_id_to_ensembl_id_file,
        min_any_count=args.min_any_count,
        zero_fraction_threshold=args.zero_fraction_threshold,
        output_scale=args.output_scale,
        zscore_by=args.zscore_by,
        prior_count=args.prior_count,
        do_tmm=not args.no_tmm,
    )
