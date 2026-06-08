#!/usr/bin/env python3
import argparse
import numpy as np
from pathlib import Path

def load_vector_csv(path: str) -> np.ndarray:
    """
    Load a single embedding vector from a CSV file.
    - Accepts files that have either:
      * exactly one row of values, or
      * a header row (e.g., 0,1,2,...) followed by one row of values.
    - If multiple rows are present, the LAST non-empty row is used.
    """
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f"File not found: {path}")

    # genfromtxt is tolerant of mixed int/float and missing values.
    arr = np.genfromtxt(path, delimiter=",")
    if arr.ndim == 0:
        raise ValueError(f"Could not parse any numeric data from {path}")

    # If two (or more) rows are present, use the last row as the vector.
    vec = arr[-1] if arr.ndim > 1 else arr
    vec = np.asarray(vec, dtype=np.float64)

    if np.isnan(vec).any():
        raise ValueError(f"Found NaNs in vector parsed from {path}")

    return vec

def cosine_similarity(a: np.ndarray, b: np.ndarray) -> float:
    a = np.asarray(a, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)
    if a.shape != b.shape:
        raise ValueError(f"Vector length mismatch: {a.shape} vs {b.shape}")
    na = np.linalg.norm(a)
    nb = np.linalg.norm(b)
    if na == 0.0 or nb == 0.0:
        raise ValueError("Cosine similarity undefined for zero-norm vector(s).")
    return float(np.dot(a, b) / (na * nb))

def main():
    parser = argparse.ArgumentParser(
        description="Compute cosine similarity between two embedding CSVs."
    )
    parser.add_argument("csv1", help="Path to first CSV")
    parser.add_argument("csv2", help="Path to second CSV")
    args = parser.parse_args()

    v1 = load_vector_csv(args.csv1)
    v2 = load_vector_csv(args.csv2)

    cos = cosine_similarity(v1, v2)
    print(f"Cosine similarity: {cos:.8f}")

if __name__ == "__main__":
    main()
