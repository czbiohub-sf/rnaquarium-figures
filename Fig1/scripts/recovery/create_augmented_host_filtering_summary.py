#!/usr/bin/env python3
"""
Create augmented host-filtering summary that includes dropout run statistics.
"""

import argparse
import polars as pl
from pathlib import Path

# Paths
_parser = argparse.ArgumentParser(description=__doc__)
_parser.add_argument(
    "--original-summary", type=Path,
    default=Path("data/75k_unstable/host-filtering.summary.txt"),
    help="Original pipeline host-filtering summary to augment with dropout stats.",
)
_args, _ = _parser.parse_known_args()

ORIGINAL_SUMMARY = _args.original_summary
STATS_WITH_DROPOUTS = Path("data/75k_unstable/stats-with-dropouts-enhanced.csv")
SEQDETECTIVE_AUGMENTED = Path("data/75k_unstable/seq-detective-judgement-summary-augmented.txt")
OUTPUT_PATH = Path("data/75k_unstable/host-filtering.summary.after-recovery.txt")


def load_original_summary():
    """Load original host-filtering summary."""
    print("Loading original host-filtering summary...")

    summary = {}
    with open(ORIGINAL_SUMMARY) as f:
        for line in f:
            line = line.strip()
            if not line or "\t" not in line:
                continue

            key, value = line.split("\t", 1)
            try:
                summary[key] = int(value)
            except ValueError:
                summary[key] = value

    print(f"  Original summary loaded with {len(summary)} entries")
    return summary


def calculate_augmented_summary(stats_df, original_summary):
    """Calculate augmented summary statistics."""
    print("\nCalculating augmented summary...")

    # Count runs with data at each stage
    stats = {
        "input_accessions": original_summary.get("input_accessions", 0),
        "local_fastq_runs": original_summary.get("local_fastq_runs", 0),
    }

    # Runs with various stages available
    stats["runs_available_ok"] = stats_df.filter(
        pl.col("starting_reads").is_not_null()
    ).height

    # Only count runs with complete unmapped output (final_reads column exists and not null)
    stats["runs_in_unmapped"] = stats_df.filter(
        pl.col("final_reads").is_not_null()
    ).height

    # For dropouts, count how many reached each stage
    stats["runs_with_fastp"] = stats_df.filter(
        pl.col("fastp_reads_after").is_not_null()
    ).height

    stats["runs_with_kallisto"] = stats_df.filter(
        pl.col("kallisto_aligned").is_not_null()
    ).height

    stats["runs_with_hisat2"] = stats_df.filter(
        pl.col("hisat2_aligned").is_not_null()
    ).height

    stats["runs_with_star"] = stats_df.filter(
        pl.col("star_aligned_unique").is_not_null()
    ).height

    stats["runs_with_bowtie2"] = stats_df.filter(
        pl.col("bowtie2_aligned").is_not_null()
    ).height

    stats["runs_with_dedup"] = stats_df.filter(
        pl.col("dedup_reads_after").is_not_null()
    ).height

    # Sum read counts across all runs (including dropouts)
    # Use null-safe sum
    stats["num_starting_reads"] = int(stats_df["starting_reads"].fill_null(0).sum())
    stats["fastp_reads_after"] = int(stats_df["fastp_reads_after"].fill_null(0).sum())
    stats["kallisto_reads_after"] = int(stats_df["kallisto_unaligned"].fill_null(0).sum())  # Unaligned = what goes to next stage
    stats["hisat2_reads_after"] = int(stats_df["hisat2_unaligned"].fill_null(0).sum())
    stats["star_reads_after"] = int(stats_df["star_unaligned"].fill_null(0).sum())
    stats["bowtie_reads_after"] = int(stats_df["bowtie2_unaligned"].fill_null(0).sum())
    stats["dedup_reads_after"] = int(stats_df["dedup_reads_after"].fill_null(0).sum())
    stats["num_unmapped_reads"] = int(stats_df["final_reads"].fill_null(0).sum())

    # Calculate dropout counts
    stats["dropouts_after_download"] = stats_df.filter(
        pl.col("starting_reads").is_not_null() &
        pl.col("fastp_reads_after").is_null()
    ).height

    stats["dropouts_after_fastp"] = stats_df.filter(
        pl.col("fastp_reads_after").is_not_null() &
        pl.col("kallisto_aligned").is_null()
    ).height

    stats["dropouts_after_kallisto"] = stats_df.filter(
        pl.col("kallisto_aligned").is_not_null() &
        pl.col("hisat2_aligned").is_null()
    ).height

    stats["dropouts_after_hisat2"] = stats_df.filter(
        pl.col("hisat2_aligned").is_not_null() &
        pl.col("star_aligned_unique").is_null()
    ).height

    stats["dropouts_after_star"] = stats_df.filter(
        pl.col("star_aligned_unique").is_not_null() &
        pl.col("bowtie2_aligned").is_null()
    ).height

    stats["dropouts_after_bowtie2"] = stats_df.filter(
        pl.col("bowtie2_aligned").is_not_null() &
        pl.col("dedup_reads_after").is_null()
    ).height

    stats["dropouts_after_dedup"] = stats_df.filter(
        pl.col("dedup_reads_after").is_not_null() &
        pl.col("final_reads").is_null()
    ).height

    print(f"  Calculated {len(stats)} summary statistics")

    return stats


def save_augmented_summary(stats, output_path):
    """Save augmented summary."""
    print(f"\nSaving augmented summary to {output_path}...")

    with open(output_path, "w") as f:
        # Write in same format as original
        for key, value in stats.items():
            f.write(f"{key}\t{value}\n")

    print(f"  Saved {len(stats)} entries")


def print_comparison(original, augmented):
    """Print comparison of original vs augmented."""
    print("\n" + "="*70)
    print("COMPARISON: ORIGINAL vs AUGMENTED")
    print("="*70)
    print(f"{'Metric':<30s} {'Original':<15s} {'Augmented':<15s} {'Change':<15s}")
    print("-"*70)

    # Compare key metrics
    key_metrics = [
        "runs_available_ok",
        "runs_in_unmapped",
        "num_starting_reads",
        "fastp_reads_after",
        "kallisto_reads_after",
        "hisat2_reads_after",
        "star_reads_after",
        "bowtie_reads_after",
        "dedup_reads_after",
        "num_unmapped_reads",
    ]

    for metric in key_metrics:
        orig_val = original.get(metric, 0)
        aug_val = augmented.get(metric, 0)
        change = aug_val - orig_val

        # Format with commas for readability
        if isinstance(orig_val, int):
            orig_str = f"{orig_val:,}"
            aug_str = f"{aug_val:,}"
            change_str = f"+{change:,}" if change >= 0 else f"{change:,}"
        else:
            orig_str = str(orig_val)
            aug_str = str(aug_val)
            change_str = f"{change}"

        print(f"{metric:<30s} {orig_str:<15s} {aug_str:<15s} {change_str:<15s}")

    # Show dropout statistics
    print("\n" + "="*70)
    print("DROPOUT STATISTICS (NEW IN AUGMENTED)")
    print("="*70)
    dropout_metrics = [
        "dropouts_after_download",
        "dropouts_after_fastp",
        "dropouts_after_kallisto",
        "dropouts_after_hisat2",
        "dropouts_after_star",
        "dropouts_after_bowtie2",
        "dropouts_after_dedup",
    ]

    for metric in dropout_metrics:
        value = augmented.get(metric, 0)
        print(f"{metric:<40s} {value:,}")


def main():
    print("Creating augmented host-filtering summary")
    print("="*70)

    # Load original summary
    original_summary = load_original_summary()

    # Load augmented stats
    print("\nLoading stats-with-dropouts.csv...")
    stats_df = pl.read_csv(STATS_WITH_DROPOUTS)
    print(f"  Loaded {stats_df.height} runs")

    # Calculate augmented summary
    augmented_summary = calculate_augmented_summary(stats_df, original_summary)

    # Save
    save_augmented_summary(augmented_summary, OUTPUT_PATH)

    # Print comparison
    print_comparison(original_summary, augmented_summary)

    print("\n" + "="*70)
    print("Output file: " + str(OUTPUT_PATH))
    print("="*70)


if __name__ == "__main__":
    main()
