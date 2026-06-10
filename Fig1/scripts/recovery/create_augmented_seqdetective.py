#!/usr/bin/env python3
"""
Create augmented seq-detective judgement summary that includes recovered entries.
Reads from both the original seq-detective summary and newly recovered download directories.
"""

import argparse
import polars as pl
from pathlib import Path

# Paths
ORIGINAL_SEQDETECTIVE = Path("data/75k_unstable/seq-detective-judgement-summary-all.txt")
RECOVERABLE_CSV = Path("data/75k_unstable/recoverable_files.csv")
OUTPUT_PATH = Path("data/75k_unstable/seq-detective-judgement-summary-augmented.txt")

_parser = argparse.ArgumentParser(description=__doc__)
_parser.add_argument(
    "--download-base", type=Path,
    default=Path("data/75k_unstable/host_mapping/download"),
    help="Pipeline output download/ tree holding recovered seq-detective judgements.",
)
_args, _ = _parser.parse_known_args()

DOWNLOAD_BASE = _args.download_base


def load_original_seqdetective():
    """Load original seq-detective summary."""
    print("Loading original seq-detective data...")

    df = pl.read_csv(
        ORIGINAL_SEQDETECTIVE,
        separator="\t",
        has_header=False,
        new_columns=["id", "file1", "file2", "judgement1", "judgement2", "reason"]
    )

    # Strip _subsample suffix to get run ID
    df = df.with_columns(
        pl.col("id").str.strip_suffix("_subsample").alias("run_id")
    )

    print(f"  Original entries: {df.height}")
    return df


def load_recovered_seqdetective():
    """Load seq-detective judgements from recovered download directories."""
    print("\nLoading recovered seq-detective judgements...")

    # Load recoverable files
    recoverable = pl.read_csv(RECOVERABLE_CSV)

    # Filter to runs with download work directories
    with_download = recoverable.filter(
        (pl.col("work_dir_download") != "") &
        (pl.col("work_dir_download").is_not_null())
    )

    print(f"  Checking {with_download.height} runs with download directories...")

    recovered_entries = []
    found = 0
    missing = 0

    for idx, row in enumerate(with_download.iter_rows(named=True)):
        if (idx + 1) % 500 == 0:
            print(f"    Progress: {idx + 1}/{with_download.height} runs ({found} found, {missing} missing)")

        run_id = row["run_id"]

        # Check if seq-detective judgement file exists in download directory
        judgement_file = DOWNLOAD_BASE / run_id / "seq-detective-judgement.txt"

        if judgement_file.exists():
            try:
                with open(judgement_file) as f:
                    line = f.read().strip()

                # Parse tab-separated line
                parts = line.split("\t")
                if len(parts) >= 5:
                    recovered_entries.append({
                        "id": parts[0],
                        "file1": parts[1],
                        "file2": parts[2] if len(parts) > 2 else "",
                        "judgement1": parts[3] if len(parts) > 3 else "",
                        "judgement2": parts[4] if len(parts) > 4 else "",
                        "reason": parts[5] if len(parts) > 5 else "",
                        "run_id": run_id
                    })
                    found += 1
            except Exception as e:
                print(f"    WARNING: Failed to read {judgement_file}: {e}")
                missing += 1
        else:
            missing += 1

    print(f"\n  Recovered entries: {found}")
    print(f"  Missing entries: {missing}")

    if not recovered_entries:
        return pl.DataFrame(schema={
            "id": pl.Utf8,
            "file1": pl.Utf8,
            "file2": pl.Utf8,
            "judgement1": pl.Utf8,
            "judgement2": pl.Utf8,
            "reason": pl.Utf8,
            "run_id": pl.Utf8
        })

    return pl.DataFrame(recovered_entries)


def merge_seqdetective_data(original_df, recovered_df):
    """Merge original and recovered seq-detective data."""
    print("\nMerging seq-detective data...")

    # Get existing run IDs
    existing_ids = set(original_df["run_id"].to_list())
    print(f"  Original unique runs: {len(existing_ids)}")

    # Filter recovered to only new runs
    new_entries = recovered_df.filter(
        ~pl.col("run_id").is_in(list(existing_ids))
    )
    print(f"  New entries from recovery: {new_entries.height}")

    # Combine
    augmented = pl.concat([original_df, new_entries])

    # Sort by run_id
    augmented = augmented.sort("run_id")

    print(f"  Total augmented entries: {augmented.height}")

    return augmented


def save_augmented_seqdetective(df, output_path):
    """Save augmented seq-detective summary."""
    print(f"\nSaving augmented seq-detective data to {output_path}...")

    # Select columns in original format (without run_id)
    output_df = df.select(["id", "file1", "file2", "judgement1", "judgement2", "reason"])

    # Write as TSV without header
    output_df.write_csv(
        output_path,
        separator="\t",
        include_header=False
    )

    print(f"  Saved {output_df.height} entries")


def main():
    print("Creating augmented seq-detective judgement summary")
    print("="*70)

    # Load data
    original_df = load_original_seqdetective()
    recovered_df = load_recovered_seqdetective()

    # Merge
    augmented_df = merge_seqdetective_data(original_df, recovered_df)

    # Save
    save_augmented_seqdetective(augmented_df, OUTPUT_PATH)

    print("\n" + "="*70)
    print("SUMMARY")
    print("="*70)
    print(f"Original entries: {original_df.height}")
    print(f"Recovered entries (new): {augmented_df.height - original_df.height}")
    print(f"Total augmented entries: {augmented_df.height}")
    print(f"Output file: {OUTPUT_PATH}")


if __name__ == "__main__":
    main()
