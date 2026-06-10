#!/usr/bin/env python3
"""
Find which dropout runs have recoverable files in nextflow work directories.
NFS-optimized: builds index once, then looks up hashes.
"""

import argparse
import polars as pl
from pathlib import Path
import os
import json
from collections import defaultdict

# Data paths
STATS_PATH = Path("data/75k_unstable/stats-merged.csv")
SEQDETECTIVE_PATH = Path("data/75k_unstable/seq-detective-judgement-summary-all.txt")
TRACE_PATH = Path("data/75k_unstable/trace-merged-dangerously.txt")

_parser = argparse.ArgumentParser(description=__doc__)
_parser.add_argument(
    "--nf-tmp-base", type=Path,
    default=Path("data/75k_unstable/nf_work"),
    help="Root of the Nextflow work-directory tree to scan for recoverable files.",
)
_args, _ = _parser.parse_known_args()

NF_TMP_BASE = _args.nf_tmp_base


def load_data():
    """Load datasets."""
    print("Loading datasets...")

    # Load stats to identify complete runs
    stats = pl.read_csv(STATS_PATH)
    stats_ids = set(stats["id"].to_list())
    print(f"  Stats (complete runs): {len(stats_ids)}")

    # Load seq-detective to identify runs that started
    sd = pl.read_csv(
        SEQDETECTIVE_PATH,
        separator="\t",
        has_header=False,
        new_columns=["id", "file1", "file2", "judgement1", "judgement2", "reason"]
    ).with_columns(
        pl.col("id").str.strip_suffix("_subsample")
    )
    print(f"  Seq-Detective: {sd.height}")

    # Load trace with all task information
    trace = pl.read_csv(TRACE_PATH, separator="\t", null_values="-")
    print(f"  Trace: {trace.height} task entries")

    return stats_ids, sd, trace


def build_work_dir_index():
    """
    Build index of all work directories once (NFS-friendly).
    Returns dict mapping truncated hash -> full path.

    Directory structure:
      /nf_tmp/XX/YYYYYY... where hash in trace is "XX/YYYYYY"
    """
    print(f"\n{'='*70}")
    print("BUILDING WORK DIRECTORY INDEX")
    print(f"{'='*70}")
    print(f"Scanning {NF_TMP_BASE}...")
    print("(This is a one-time scan of the NFS directory tree)")

    index = {}

    # List all first-level subdirectories (2-char prefixes)
    try:
        prefixes = [d for d in os.listdir(NF_TMP_BASE)
                   if os.path.isdir(NF_TMP_BASE / d) and len(d) == 2]
        print(f"Found {len(prefixes)} prefix directories")
    except Exception as e:
        print(f"Error listing {NF_TMP_BASE}: {e}")
        return index

    # For each prefix, list all work directories
    for i, prefix in enumerate(sorted(prefixes)):
        if (i + 1) % 16 == 0:
            print(f"  Scanned {i + 1}/{len(prefixes)} prefixes... ({len(index)} work dirs found)")

        prefix_dir = NF_TMP_BASE / prefix
        try:
            # List all directories under this prefix
            work_dirs = [d for d in os.listdir(prefix_dir)
                        if os.path.isdir(prefix_dir / d)]

            for work_dir in work_dirs:
                # Truncated hash is "prefix/first_6_chars"
                # Full hash in directory name might be longer
                truncated = f"{prefix}/{work_dir[:6]}"
                full_path = str(prefix_dir / work_dir)

                # Store mapping (there might be collisions, but unlikely)
                if truncated in index:
                    # If collision, keep the first one found
                    pass
                else:
                    index[truncated] = full_path

        except Exception as e:
            print(f"Error listing {prefix_dir}: {e}")
            continue

    print(f"\nCompleted: indexed {len(index)} work directories")
    return index


def analyze_recoverable_files(stats_ids, sd, trace, work_dir_index):
    """Analyze which files can be recovered from work directories."""

    # Get dropout runs (in seq-detective but not in final stats)
    dropout_runs = sd.filter(~pl.col("id").is_in(list(stats_ids)))
    dropout_ids = set(dropout_runs["id"].to_list())

    print(f"\n{'='*70}")
    print(f"ANALYZING {len(dropout_ids)} DROPOUT RUNS")
    print(f"{'='*70}")

    # Get trace entries for dropout runs
    dropout_trace = trace.filter(
        pl.col("tag").is_in(list(dropout_ids))
    )

    print(f"Trace entries for dropout runs: {dropout_trace.height}")

    # Analyze by run
    results = []
    dirs_found = 0
    dirs_searched = 0

    # Track overall statistics
    process_file_counts = defaultdict(int)
    process_found_counts = defaultdict(int)

    for idx, run_id in enumerate(sorted(dropout_ids)):
        if (idx + 1) % 500 == 0:
            print(f"  Processed {idx + 1}/{len(dropout_ids)} runs... ({dirs_found} work dirs found)")

        run_trace = dropout_trace.filter(pl.col("tag") == run_id)

        if run_trace.height == 0:
            sd_row = dropout_runs.filter(pl.col("id") == run_id)
            results.append({
                "run_id": run_id,
                "completed_processes": [],
                "failed_processes": [],
                "work_dirs_found": {},
                "work_dirs_missing": {},
                "has_recoverable_files": False,
                "seq_detective_judgement": sd_row["judgement1"][0] if sd_row.height > 0 else "",
                "seq_detective_reason": sd_row["reason"][0] if sd_row.height > 0 else "",
            })
            continue

        # Identify completed and failed processes
        completed = run_trace.filter(
            pl.col("status").is_in(["COMPLETED", "CACHED"])
        )
        failed = run_trace.filter(
            pl.col("status").is_in(["FAILED", "ABORTED"])
        )

        completed_processes = completed["process"].unique().to_list() if completed.height > 0 else []
        failed_processes = failed["process"].unique().to_list() if failed.height > 0 else []

        # Find work directories using pre-built index
        work_dirs_found = {}
        work_dirs_missing = {}

        for row in run_trace.iter_rows(named=True):
            if row["hash"] and row["status"] in ["COMPLETED", "CACHED", "FAILED"]:
                process = row["process"]
                hash_val = row["hash"]

                dirs_searched += 1

                # Look up in index (fast dictionary lookup, no filesystem access)
                if hash_val in work_dir_index:
                    work_dirs_found[process] = work_dir_index[hash_val]
                    dirs_found += 1
                    process_found_counts[process] += 1
                else:
                    work_dirs_missing[process] = hash_val

                process_file_counts[process] += 1

        # Get seq-detective info
        sd_row = dropout_runs.filter(pl.col("id") == run_id)
        judgement = f"{sd_row['judgement1'][0]}{sd_row['judgement2'][0] or ''}"
        reason = sd_row["reason"][0]

        results.append({
            "run_id": run_id,
            "completed_processes": completed_processes,
            "failed_processes": failed_processes,
            "work_dirs_found": work_dirs_found,
            "work_dirs_missing": work_dirs_missing,
            "has_recoverable_files": bool(work_dirs_found),
            "seq_detective_judgement": judgement,
            "seq_detective_reason": reason,
        })

    print(f"\n  Completed: {len(dropout_ids)}/{len(dropout_ids)} runs")
    print(f"  Total work directories searched: {dirs_searched}")
    print(f"  Total work directories found: {dirs_found}")

    return results, process_file_counts, process_found_counts


def summarize_results(results, process_file_counts, process_found_counts):
    """Print summary statistics."""

    print(f"\n{'='*70}")
    print("SUMMARY")
    print(f"{'='*70}")

    # Overall recovery statistics
    with_files = sum(1 for r in results if r["has_recoverable_files"])
    without_files = len(results) - with_files

    print(f"\nRuns with recoverable files: {with_files}/{len(results)} ({100*with_files/len(results):.1f}%)")
    print(f"Runs without recoverable files: {without_files}/{len(results)} ({100*without_files/len(results):.1f}%)")

    # Recovery by process type
    print(f"\n{'='*70}")
    print("FILE RECOVERY BY PROCESS TYPE")
    print(f"{'='*70}")
    print(f"{'Process':<25s} {'Found':<10s} {'Total':<10s} {'Recovery %':<10s}")
    print("-" * 70)

    for process in sorted(process_file_counts.keys(), key=lambda p: process_found_counts.get(p, 0), reverse=True):
        found = process_found_counts.get(process, 0)
        total = process_file_counts[process]
        pct = 100 * found / total if total > 0 else 0
        print(f"{process:<25s} {found:<10d} {total:<10d} {pct:<10.1f}")

    # Categorize by number of recovered work dirs
    recovery_levels = defaultdict(int)
    for r in results:
        n_found = len(r["work_dirs_found"])
        recovery_levels[n_found] += 1

    print(f"\n{'='*70}")
    print("RECOVERY LEVELS")
    print(f"{'='*70}")
    print(f"{'Work Dirs Found':<20s} {'Number of Runs':<15s}")
    print("-" * 70)
    for n_dirs in sorted(recovery_levels.keys(), reverse=True):
        print(f"{n_dirs:<20d} {recovery_levels[n_dirs]:<15d}")


def save_results(results):
    """Save detailed results to files."""

    output_dir = Path("data/75k_unstable")
    output_dir.mkdir(parents=True, exist_ok=True)

    # Save as JSON for detailed inspection
    json_path = output_dir / "recoverable_files_detailed.json"
    with open(json_path, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nSaved detailed results to {json_path}")

    # Create a flat table for easier analysis
    flat_records = []
    for r in results:
        base_record = {
            "run_id": r["run_id"],
            "seq_detective_judgement": r["seq_detective_judgement"],
            "seq_detective_reason": r["seq_detective_reason"],
            "num_completed_processes": len(r["completed_processes"]),
            "num_work_dirs_found": len(r["work_dirs_found"]),
            "num_work_dirs_missing": len(r["work_dirs_missing"]),
            "has_recoverable_files": r["has_recoverable_files"],
            "completed_processes": ",".join(r["completed_processes"]),
            "processes_with_files": ",".join(r["work_dirs_found"].keys()),
        }

        # Add work directory paths for key stages
        for stage in ["download", "fastp", "kb_negative", "star_counts", "host_cram", "stats_csv"]:
            base_record[f"work_dir_{stage}"] = r["work_dirs_found"].get(stage, "")

        flat_records.append(base_record)

    # Save as CSV
    df = pl.DataFrame(flat_records)
    csv_path = output_dir / "recoverable_files.csv"
    df.write_csv(csv_path)
    print(f"Saved flat table to {csv_path}")

    # Save as Parquet
    parquet_path = output_dir / "recoverable_files.parquet"
    df.write_parquet(parquet_path)
    print(f"Saved flat table to {parquet_path}")

    return df


def main():
    # Load data
    stats_ids, sd, trace = load_data()

    # Build work directory index (single NFS traversal)
    work_dir_index = build_work_dir_index()

    # Analyze using the index (no additional filesystem access)
    results, process_file_counts, process_found_counts = analyze_recoverable_files(
        stats_ids, sd, trace, work_dir_index
    )

    # Summarize and save
    summarize_results(results, process_file_counts, process_found_counts)
    df = save_results(results)

    print(f"\n{'='*70}")
    print("NEXT STEPS")
    print(f"{'='*70}")
    print("\n1. Review recoverable_files.csv to prioritize recovery efforts")
    print("2. Focus on runs with has_recoverable_files=true")
    print("3. Use work_dir_* columns to locate specific output files")
    print("4. Check which processes have highest recovery rates")


if __name__ == "__main__":
    main()
