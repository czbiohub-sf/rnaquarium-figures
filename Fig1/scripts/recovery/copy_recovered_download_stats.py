#!/usr/bin/env python3
"""
Copy recovered download stats files to versioned output directory.
Checks for existing files and only copies missing ones to avoid clobbering.
"""

import argparse
import polars as pl
from pathlib import Path
import shutil
import os
from datetime import datetime

# Paths
RECOVERABLE_CSV = Path("data/75k_unstable/recoverable_files.csv")
LOG_FILE = Path("data/75k_unstable/download_copy.log")

_parser = argparse.ArgumentParser(description=__doc__)
_parser.add_argument(
    "--dest-base", type=Path,
    default=Path("data/75k_unstable/host_mapping/download"),
    help="Pipeline output download/ tree to copy recovered per-run stats into.",
)
_args, _ = _parser.parse_known_args()

DEST_BASE = _args.dest_base

# Files to copy from download work directories
DOWNLOAD_FILES = [
    "info.txt",
    "stats.txt",
    "seq-detective-judgement.txt",
    "seq-detective-stats.json"
]


def log_message(msg, log_file):
    """Write message to both stdout and log file."""
    timestamp = datetime.now().isoformat()
    log_line = f"[{timestamp}] {msg}"
    print(log_line)
    with open(log_file, "a") as f:
        f.write(log_line + "\n")


def copy_download_files():
    """Copy download stats files from recovered work directories."""

    # Initialize log
    with open(LOG_FILE, "w") as f:
        f.write(f"Download File Recovery Log - Started at {datetime.now().isoformat()}\n")
        f.write("="*70 + "\n\n")

    log_message("Loading recoverable files data...", LOG_FILE)
    df = pl.read_csv(RECOVERABLE_CSV)

    # Filter to runs with download work directories
    with_download = df.filter(
        (pl.col("work_dir_download") != "") &
        (pl.col("work_dir_download").is_not_null())
    )

    log_message(f"Found {with_download.height} runs with recoverable download work directories", LOG_FILE)

    # Track statistics
    stats = {
        "total_runs": with_download.height,
        "runs_processed": 0,
        "runs_skipped_existing": 0,
        "runs_copied": 0,
        "files_copied": 0,
        "files_skipped": 0,
        "errors": 0
    }

    # Process each run
    for idx, row in enumerate(with_download.iter_rows(named=True)):
        run_id = row["run_id"]
        work_dir = Path(row["work_dir_download"])

        if (idx + 1) % 100 == 0:
            log_message(
                f"Progress: {idx + 1}/{with_download.height} runs | "
                f"Copied: {stats['runs_copied']} runs, {stats['files_copied']} files | "
                f"Errors: {stats['errors']}",
                LOG_FILE
            )

        stats["runs_processed"] += 1

        # Check if source work directory exists
        if not work_dir.exists():
            log_message(f"WARNING: Work directory not found: {work_dir}", LOG_FILE)
            stats["errors"] += 1
            continue

        # Destination directory
        dest_dir = DEST_BASE / run_id

        # Skip if destination directory already exists (assume it's complete)
        if dest_dir.exists():
            stats["runs_skipped_existing"] += 1
            stats["files_skipped"] += len(DOWNLOAD_FILES)
            continue

        # Create destination directory
        try:
            dest_dir.mkdir(parents=True, exist_ok=True)
        except Exception as e:
            log_message(f"ERROR: Failed to create {dest_dir}: {e}", LOG_FILE)
            stats["errors"] += 1
            continue

        # Copy all files for this run
        run_copied = False
        for filename in DOWNLOAD_FILES:
            src_file = work_dir / filename
            dest_file = dest_dir / filename

            # Check if source exists
            if not src_file.exists():
                log_message(f"WARNING: Source file not found: {src_file}", LOG_FILE)
                stats["errors"] += 1
                continue

            # Copy the file
            try:
                shutil.copy2(src_file, dest_file)
                stats["files_copied"] += 1
                run_copied = True
            except Exception as e:
                log_message(f"ERROR: Failed to copy {src_file} to {dest_file}: {e}", LOG_FILE)
                stats["errors"] += 1

        if run_copied:
            stats["runs_copied"] += 1

    # Final summary
    log_message("\n" + "="*70, LOG_FILE)
    log_message("COPY SUMMARY", LOG_FILE)
    log_message("="*70, LOG_FILE)
    log_message(f"Total runs with download work dirs: {stats['total_runs']}", LOG_FILE)
    log_message(f"Runs processed: {stats['runs_processed']}", LOG_FILE)
    log_message(f"Runs with all files already existing: {stats['runs_skipped_existing']}", LOG_FILE)
    log_message(f"Runs with files copied: {stats['runs_copied']}", LOG_FILE)
    log_message(f"Total files copied: {stats['files_copied']}", LOG_FILE)
    log_message(f"Total files skipped (already exist): {stats['files_skipped']}", LOG_FILE)
    log_message(f"Errors encountered: {stats['errors']}", LOG_FILE)
    log_message(f"\nCompleted at {datetime.now().isoformat()}", LOG_FILE)

    return stats


if __name__ == "__main__":
    print("Starting download file recovery...")
    print("This may take a while due to NFS operations.")
    print(f"Progress will be logged to: {LOG_FILE}")
    print()

    stats = copy_download_files()

    print("\n" + "="*70)
    print("Recovery complete!")
    print(f"Files copied: {stats['files_copied']}")
    print(f"Runs updated: {stats['runs_copied']}")
    print(f"Runs already complete: {stats['runs_skipped_existing']}")
    print(f"Check {LOG_FILE} for detailed log")
