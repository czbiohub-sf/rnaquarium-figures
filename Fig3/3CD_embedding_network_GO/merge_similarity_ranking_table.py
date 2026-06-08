#!/usr/bin/env python3
"""
Merge four CSVs on 'gene', prefix columns by filename, and sort by the
embedding cosine similarity column (descending). Writes the merged CSV.

Usage:
  python merge_gene_tables.py file1.csv file2.csv file3.csv file4.csv -o merged.csv
"""

from __future__ import annotations
import argparse
from functools import reduce
from pathlib import Path
import pandas as pd
import numpy as np
from pandas.api.types import is_numeric_dtype, is_string_dtype


def try_convert_numeric(s: pd.Series) -> pd.Series:
    """
    Try converting to numeric without deprecated errors='ignore'.
    If conversion raises, return the original series unchanged.
    """
    if is_numeric_dtype(s):
        return s
    try:
        s2 = s
        if is_string_dtype(s2):
            s2 = s2.str.strip()
        return pd.to_numeric(s2)  # may raise on mixed values; we'll catch below
    except Exception:
        return s


def to_nullable_int(s: pd.Series) -> pd.Series:
    """
    Convert a Series to pandas nullable integer Int64, preserving missing values.
    - Coerce non-numeric to NaN.
    - If any non-NaN values have fractional parts, round to nearest integer.
    """
    x = pd.to_numeric(s, errors="coerce")
    # Identify fractional values (e.g., 3.2). Ignore NaNs.
    frac_mask = (x.notna()) & ((x % 1) != 0)
    if frac_mask.any():
        x = x.round()
    return x.astype("Int64")  # nullable integer dtype


def is_rank_col(colname: str) -> bool:
    # After prefixing, rank columns look like "<stem>_rank".
    # Use a case-insensitive check on the last token.
    return colname != "gene" and colname.split("_")[-1].lower() == "rank"


def read_with_prefix(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    if "gene" not in df.columns:
        raise ValueError(f"'gene' column not found in {path}")

    # Normalize column names
    df.columns = [c.strip() for c in df.columns]
    df["gene"] = df["gene"].astype(str)

    prefix = Path(path).stem
    prefix = prefix.removesuffix('__correlation_ranked')
    prefix = prefix.removesuffix('_similarity_ranked')
    print(prefix)
    
    rename_map = {}
    for c in df.columns:
        if c == "gene":
            continue
        cl = c.strip().lower()
        if cl.endswith("rank"):
            rename_map[c] = f"{prefix}_rank"
        elif cl.endswith("correlation"):
            rename_map[c] = f"{prefix}_correlation"
        elif cl.endswith("cosine_similarity"):
            rename_map[c] = f"{prefix}_cosine_similarity"
        else:
            print(f"edge case WARNING: prefix: {prefix}, c: {c}")
            #rename_map[c] = f"{prefix}_{c.strip()}"

    df = df.rename(columns=rename_map)

    # Generic numeric conversion attempt for non-'gene' columns
    for c in df.columns:
        if c != "gene":
            df[c] = try_convert_numeric(df[c])

    # Force rank columns to nullable integer
    for c in df.columns:
        if is_rank_col(c):
            df[c] = to_nullable_int(df[c])

    return df


def main():
    parser = argparse.ArgumentParser(description="Merge CSVs on 'gene' with filename-prefixed columns.")
    parser.add_argument("inputs", nargs=4, help="Four input CSV files")
    parser.add_argument("-o", "--output", default="merged.csv", help="Output CSV path (default: merged.csv)")
    args = parser.parse_args()

    input_paths = [Path(p) for p in args.inputs]
    for p in input_paths:
        if not p.exists():
            raise FileNotFoundError(f"Input file not found: {p}")

    dfs = [read_with_prefix(p) for p in input_paths]
    merged = reduce(lambda left, right: pd.merge(left, right, on="gene", how="outer"), dfs)

    # Find the (prefixed) embedding cosine similarity column
    embed_cols = [c for c in merged.columns if c.endswith("_embedding_cosine_similarity")]
    if not embed_cols:
        raise RuntimeError("No column ending with 'embedding_cosine_similarity' found after renaming.")
    sort_col = embed_cols[0]

    # Ensure the sort column is numeric; invalids -> NaN so they sort last
    merged[sort_col] = pd.to_numeric(merged[sort_col], errors="coerce")

    merged = merged.sort_values(by=sort_col, ascending=False, na_position="last")

    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    merged.to_csv(args.output, index=False)
    print(f"Wrote: {Path(args.output).resolve()}")


if __name__ == "__main__":
    main()
