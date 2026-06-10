#!/usr/bin/env python3
"""
Extract seq-detective metrics from per-run JSON files for contour plotting.
Outputs a parquet file with mapping_rate, sparsity, etc. for each mate.

Optimized to avoid NFS metadata operations - uses run ID list from judgement file.
"""

import argparse
import json
import polars as pl
from pathlib import Path
import sys

JUDGEMENT_PATH = Path("data/75k_unstable/seq-detective-judgement-summary-all.txt")
OUTPUT_PATH = Path("data/75k_unstable/seqdetective_metrics.parquet")

_parser = argparse.ArgumentParser(description=__doc__)
_parser.add_argument(
    "--seqdetective-dir", type=Path,
    default=Path("data/75k_unstable/host_mapping/download"),
    help="Pipeline output download/ tree holding per-run seq-detective JSON files.",
)
_parser.add_argument(
    "--rosetta-path", type=Path,
    default=Path("data/metadata/zf_rosetta.tsv"),
    help="Rosetta TSV mapping run IDs to accessions.",
)
_args, _ = _parser.parse_known_args()

SEQDETECTIVE_DIR = _args.seqdetective_dir
ROSETTA_PATH = _args.rosetta_path


def parse_json_file(json_path: str) -> dict | None:
    """Parse a single seq-detective JSON file and extract metrics."""
    try:
        with open(json_path) as f:
            data = json.load(f)
    except (json.JSONDecodeError, IOError, FileNotFoundError):
        return None

    accession = data.get("accession", "").replace("_subsample", "")
    readfiles = data.get("readfiles", [])

    if not readfiles:
        return None

    result = {"id": accession}

    # Extract metrics for each mate
    for i, rf in enumerate(readfiles[:2]):  # max 2 mates
        mate_suffix = f"_m{i+1}" if len(readfiles) > 1 else ""
        mapping = rf.get("mapping", {})
        qc = mapping.get("qc", {})

        result[f"mapping_rate{mate_suffix}"] = qc.get("mapping_rate")
        result[f"nofeature_rate{mate_suffix}"] = qc.get("nofeature_rate")
        result[f"sparsity{mate_suffix}"] = qc.get("sparsity")
        result[f"pos_strand_rate{mate_suffix}"] = qc.get("pos_strand_rate")
        result[f"readlen{mate_suffix}"] = rf.get("readlen")

    result["n_mates"] = len(readfiles)
    return result


def main():
    # Get run IDs from judgement file (avoids expensive directory scan)
    print(f"Loading run IDs from {JUDGEMENT_PATH}...", flush=True)
    sd = pl.read_csv(
        JUDGEMENT_PATH,
        separator="\t",
        has_header=False,
        new_columns=["id", "file1", "file2", "j1", "j2", "reason"]
    )
    # Strip _subsample suffix to get actual run IDs
    run_ids = sd.with_columns(
        pl.col("id").str.strip_suffix("_subsample")
    )["id"].unique().to_list()

    print(f"Found {len(run_ids)} unique run IDs", flush=True)

    # Process each run by constructing path directly (no metadata calls)
    all_results = []
    missing = 0

    for i, run_id in enumerate(run_ids):
        json_path = f"{SEQDETECTIVE_DIR}/{run_id}/seq-detective-stats.json"
        result = parse_json_file(json_path)
        if result:
            all_results.append(result)
        else:
            missing += 1

        if (i + 1) % 10000 == 0:
            print(f"  Processed {i + 1}/{len(run_ids)}, extracted {len(all_results)}, missing {missing}...", flush=True)

    print(f"Total: extracted {len(all_results)} runs, {missing} missing JSON files", flush=True)

    # Convert to polars DataFrame
    df = pl.DataFrame(all_results)

    # Join with rosetta for bioproject mapping
    rosetta = pl.read_csv(ROSETTA_PATH, separator="\t")
    df = df.join(
        rosetta.select(["run.accession", "bioproject"]),
        left_on="id",
        right_on="run.accession",
        how="left"
    )
    print(f"Joined with rosetta: {df.filter(pl.col('bioproject').is_not_null()).height} runs have bioproject", flush=True)

    # Save to parquet
    OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    df.write_parquet(OUTPUT_PATH)
    print(f"Saved to {OUTPUT_PATH}", flush=True)

    # Print summary
    print("\nSummary:", flush=True)
    print(f"  Total runs: {df.height}", flush=True)
    print(f"  Single-end (1 mate): {df.filter(pl.col('n_mates') == 1).height}", flush=True)
    print(f"  Paired-end (2 mates): {df.filter(pl.col('n_mates') == 2).height}", flush=True)

    # Print column null counts
    print("\nColumn null counts:", flush=True)
    for col in df.columns:
        null_count = df[col].null_count()
        if null_count > 0:
            print(f"  {col}: {null_count} nulls", flush=True)


if __name__ == "__main__":
    main()
