#!/usr/bin/env python3
"""
Create augmented stats CSV with full statistics from all recovered dropout stages.
Parses stats from all available work directories (download, fastp, kallisto, aligners, etc.)
"""

import polars as pl
from pathlib import Path
import json
import re

# Paths
STATS_MERGED = Path("data/75k_unstable/stats-merged.csv")
RECOVERABLE_JSON = Path("data/75k_unstable/recoverable_files_detailed.json")
OUTPUT_PATH = Path("data/75k_unstable/stats-with-dropouts-enhanced.csv")


def parse_sam_flags(stats_file):
    """Parse SAM flag statistics file (hisat2, bowtie2, gsnap)."""
    if not stats_file.exists():
        return None

    total = multi = aligned = unaligned = mixed = 0

    try:
        with open(stats_file) as f:
            for line in f:
                parts = line.strip().split()
                if len(parts) != 2:
                    continue

                count = int(parts[0])
                sambits = int(parts[1])

                # Ignore supplementary, QC fail, second-in-pair
                if (sambits & 0x80) or (sambits & 0x200) or (sambits & 0x800):
                    continue

                # Track secondary mappings
                if sambits & 0x100:
                    multi += count
                    continue

                total += count

                # Unaligned
                if ((sambits & 0x5) == 0x4) or ((sambits & 0xD) > 0x4):
                    unaligned += count

                # Aligned (concordant)
                if not ((sambits & 0x4) or (sambits & 0x8)):
                    aligned += count

                # Mixed
                if ((sambits & 0xD) == 0x5) or ((sambits & 0xD) == 0x9):
                    mixed += count

        unique = aligned - multi
        return {
            "total": total,
            "aligned": aligned,
            "multi": multi,
            "unique": unique,
            "unaligned": unaligned,
            "mixed": mixed
        }
    except Exception as e:
        print(f"  Warning: Failed to parse {stats_file}: {e}")
        return None


def parse_star_log(log_file):
    """Parse STAR Log.final.out file."""
    if not log_file.exists():
        return None

    try:
        with open(log_file) as f:
            content = f.read()

        # Extract stats using regex
        stats = {}

        patterns = {
            "reads_before": r"Number of input reads \|\s+(\d+)",
            "avg_len": r"Average input read length \|\s+(\d+)",
            "unique": r"Uniquely mapped reads number \|\s+(\d+)",
            "multialign": r"Number of reads mapped to multiple loci \|\s+(\d+)",
            "too_short": r"Number of reads unmapped: too short \|\s+(\d+)"
        }

        for key, pattern in patterns.items():
            match = re.search(pattern, content)
            if match:
                stats[key] = int(match.group(1))

        if "reads_before" in stats and "unique" in stats and "multialign" in stats:
            stats["unaligned"] = stats["reads_before"] - stats["unique"] - stats["multialign"]

        return stats
    except Exception as e:
        print(f"  Warning: Failed to parse {log_file}: {e}")
        return None


def parse_kallisto_json(json_file):
    """Parse kallisto/kb run_info.json."""
    if not json_file.exists():
        return None

    try:
        with open(json_file) as f:
            data = json.load(f)

        stats = {}
        if "n_targets" in data:
            stats["targets"] = data["n_targets"]
        if "n_processed" in data:
            stats["reads_before"] = data["n_processed"]
        if "n_pseudoaligned" in data:
            stats["aligned"] = data["n_pseudoaligned"]
        if "n_unique" in data:
            stats["aligned_unique"] = data["n_unique"]

        if "reads_before" in stats and "aligned" in stats:
            stats["unaligned"] = stats["reads_before"] - stats["aligned"]

        return stats
    except Exception as e:
        print(f"  Warning: Failed to parse {json_file}: {e}")
        return None


def parse_fastp_stats(stats_file):
    """Parse fastp stats.txt."""
    if not stats_file.exists():
        return None

    try:
        with open(stats_file) as f:
            content = f.read()

        stats = {}

        # Extract reads before/after filtering
        patterns = {
            "reads_before": r"Read1 before filtering:\s+total reads:\s+(\d+)",
            "reads_after": r"Read1 after filtering:\s+total reads:\s+(\d+)",
            "reads_too_short": r"reads failed due to too short:\s+(\d+)",
            "reads_trimmed": r"reads with adapter trimmed:\s+(\d+)"
        }

        for key, pattern in patterns.items():
            match = re.search(pattern, content)
            if match:
                stats[key] = int(match.group(1))

        return stats
    except Exception as e:
        print(f"  Warning: Failed to parse {stats_file}: {e}")
        return None


def parse_download_stats(download_dir):
    """Parse download directory stats (starting reads, fastp, etc.)."""
    result = {
        "starting_reads": None,
        "r1_median_len": None,
        "r2_median_len": None,
        "single_end": None,
        "fastp_reads_before": None,
        "fastp_reads_after": None,
        "fastp_reads_too_short": None,
        "fastp_reads_trimmed": None,
    }

    # Parse stats.txt for starting reads
    stats_file = download_dir / "stats.txt"
    if stats_file.exists():
        try:
            with open(stats_file) as f:
                content = f.read()
                match = re.search(r"spots read\s*:\s*([\d,]+)", content, re.IGNORECASE)
                if match:
                    result["starting_reads"] = int(match.group(1).replace(",", ""))
        except Exception as e:
            print(f"  Warning: Failed to parse {stats_file}: {e}")

    # Parse seq-detective-stats.json
    seqdet_file = download_dir / "seq-detective-stats.json"
    if seqdet_file.exists():
        try:
            with open(seqdet_file) as f:
                data = json.load(f)

            if "readfiles" in data and len(data["readfiles"]) > 0:
                readfile = data["readfiles"][0]

                # Read length
                if "readlen" in readfile:
                    result["r1_median_len"] = readfile["readlen"]

                # Single/paired end
                result["single_end"] = len(data["readfiles"]) == 1
                if len(data["readfiles"]) > 1:
                    result["r2_median_len"] = data["readfiles"][1].get("readlen")

        except Exception as e:
            print(f"  Warning: Failed to parse {seqdet_file}: {e}")

    return result


def parse_dedup_output(dedup_dir):
    """Parse dedup stats from .command.out file."""
    cmd_out = dedup_dir / ".command.out"
    if not cmd_out.exists():
        return None

    try:
        with open(cmd_out) as f:
            content = f.read()

        stats = {}

        # Extract stats: "total reads: 1000062" and "unique reads: 507809"
        match_total = re.search(r"total reads:\s*(\d+)", content)
        match_unique = re.search(r"unique reads:\s*(\d+)", content)

        if match_total:
            stats["reads_before"] = int(match_total.group(1))
        if match_unique:
            stats["reads_after"] = int(match_unique.group(1))

        return stats
    except Exception as e:
        print(f"  Warning: Failed to parse {cmd_out}: {e}")
        return None


def parse_dropout_run(run_id, work_dirs):
    """Parse all available stats for a dropout run."""
    stats = {
        "id": run_id,
        # Will be filled in from work directories
    }

    # Parse download directory (starting reads, fastp)
    if "download" in work_dirs:
        download_dir = Path(work_dirs["download"])
        download_stats = parse_download_stats(download_dir)
        stats.update(download_stats)

    # Parse fastp work directory
    if "fastp" in work_dirs:
        fastp_dir = Path(work_dirs["fastp"])
        fastp_stats = parse_fastp_stats(fastp_dir / "stats.txt")
        if fastp_stats:
            stats["fastp_reads_before"] = fastp_stats.get("reads_before")
            stats["fastp_reads_after"] = fastp_stats.get("reads_after")
            stats["fastp_reads_too_short"] = fastp_stats.get("reads_too_short")
            stats["fastp_reads_trimmed"] = fastp_stats.get("reads_trimmed")

    # Parse kallisto (kb_negative)
    if "kb_negative" in work_dirs:
        kb_dir = Path(work_dirs["kb_negative"])
        kb_stats = parse_kallisto_json(kb_dir / "run_info.json")
        if kb_stats:
            stats["kallisto_reads_before"] = kb_stats.get("reads_before")
            stats["kallisto_aligned"] = kb_stats.get("aligned")
            stats["kallisto_aligned_unique"] = kb_stats.get("aligned_unique")
            stats["kallisto_unaligned"] = kb_stats.get("unaligned")
            stats["kallisto_targets"] = kb_stats.get("targets")

    # Parse HISAT2
    if "hisat2" in work_dirs:
        hisat2_dir = Path(work_dirs["hisat2"])
        hisat2_stats = parse_sam_flags(hisat2_dir / "hisat2.stats.txt")
        if hisat2_stats:
            stats["hisat2_reads_before"] = hisat2_stats["total"]
            stats["hisat2_aligned"] = hisat2_stats["aligned"]
            stats["hisat2_multialign"] = hisat2_stats["multi"]
            stats["hisat2_aligned_unique"] = hisat2_stats["unique"]
            stats["hisat2_unaligned"] = hisat2_stats["unaligned"]
            stats["hisat2_mixed"] = hisat2_stats["mixed"]

    # Parse STAR
    if "star" in work_dirs:
        star_dir = Path(work_dirs["star"])
        star_stats = parse_star_log(star_dir / "Log.final.out")
        if star_stats:
            stats["star_reads_before"] = star_stats.get("reads_before")
            stats["star_avg_len"] = star_stats.get("avg_len")
            stats["star_aligned_unique"] = star_stats.get("unique")
            stats["star_multialign"] = star_stats.get("multialign")
            stats["star_unaligned"] = star_stats.get("unaligned")
            stats["star_too_short"] = star_stats.get("too_short")

    # Parse Bowtie2
    if "bowtie2" in work_dirs:
        bowtie2_dir = Path(work_dirs["bowtie2"])
        bowtie2_stats = parse_sam_flags(bowtie2_dir / "bowtie2.stats.txt")
        if bowtie2_stats:
            stats["bowtie2_reads_before"] = bowtie2_stats["total"]
            stats["bowtie2_aligned"] = bowtie2_stats["aligned"]
            stats["bowtie2_multialign"] = bowtie2_stats["multi"]
            stats["bowtie2_aligned_unique"] = bowtie2_stats["unique"]
            stats["bowtie2_unaligned"] = bowtie2_stats["unaligned"]
            stats["bowtie2_mixed"] = bowtie2_stats["mixed"]

    # Parse Dedup
    if "dedup" in work_dirs:
        dedup_dir = Path(work_dirs["dedup"])
        dedup_stats = parse_dedup_output(dedup_dir)
        if dedup_stats and "reads_before" in dedup_stats:
            stats["dedup_reads_before"] = dedup_stats["reads_before"]
            stats["dedup_reads_after"] = dedup_stats.get("reads_after")

    # Parse GSNAP
    if "gsnap" in work_dirs:
        gsnap_dir = Path(work_dirs["gsnap"])
        gsnap_stats = parse_sam_flags(gsnap_dir / "gsnap.stats.txt")
        if gsnap_stats:
            stats["gsnap_reads_before"] = gsnap_stats["total"]
            stats["gsnap_aligned"] = gsnap_stats["aligned"]
            stats["gsnap_multialign"] = gsnap_stats["multi"]
            stats["gsnap_aligned_unique"] = gsnap_stats["unique"]
            stats["gsnap_unaligned"] = gsnap_stats["unaligned"]
            stats["gsnap_mixed"] = gsnap_stats["mixed"]
            stats["final_reads"] = gsnap_stats["unaligned"]

    return stats


def main():
    print("Creating enhanced stats-with-dropouts CSV")
    print("="*70)

    # Load complete runs
    print("Loading complete runs...")
    complete_df = pl.read_csv(STATS_MERGED)
    print(f"  Complete runs: {complete_df.height}")

    # Load dropout work directories
    print("\nLoading dropout run work directories...")
    with open(RECOVERABLE_JSON) as f:
        dropout_data = json.load(f)
    print(f"  Dropout runs: {len(dropout_data)}")

    # Parse dropout stats
    print("\nParsing dropout run statistics...")
    dropout_stats = []
    success = 0
    failed = 0

    for idx, entry in enumerate(dropout_data):
        if (idx + 1) % 500 == 0:
            print(f"  Progress: {idx + 1}/{len(dropout_data)} ({success} success, {failed} failed)")

        run_id = entry["run_id"]
        work_dirs = entry["work_dirs_found"]

        if not work_dirs:
            failed += 1
            continue

        stats = parse_dropout_run(run_id, work_dirs)

        # Debug specific runs
        if run_id in ["ERR1427405", "SRR11293496"]:
            print(f"  DEBUG {run_id}: starting_reads={stats.get('starting_reads')}, dedup_before={stats.get('dedup_reads_before')}, dedup_after={stats.get('dedup_reads_after')}")

        if stats.get("starting_reads") is not None:
            dropout_stats.append(stats)
            success += 1
        else:
            failed += 1

    print(f"\n  Successfully parsed: {success} dropout runs")
    print(f"  Failed: {failed} dropout runs")

    # Create DataFrame with full schema inference
    # infer_schema_length=None tells Polars to scan all dicts to determine schema
    dropout_df = pl.DataFrame(dropout_stats, infer_schema_length=None) if dropout_stats else None

    # Debug: check if dedup columns exist
    if dropout_df is not None:
        print(f"\nDEBUG: dropout_df columns with 'dedup': {[c for c in dropout_df.columns if 'dedup' in c]}")
        err_row = dropout_df.filter(pl.col('id') == 'ERR1427405')
        if err_row.height > 0:
            print(f"DEBUG: ERR1427405 in dropout_df: dedup_before={err_row['dedup_reads_before'][0] if 'dedup_reads_before' in dropout_df.columns else 'NO COL'}")

    # Merge with complete runs
    print("\nMerging with complete runs...")

    if dropout_df is not None and not dropout_df.is_empty():
        # Add missing columns to both DataFrames
        all_cols = list(set(complete_df.columns) | set(dropout_df.columns))

        for col in all_cols:
            if col not in complete_df.columns:
                complete_df = complete_df.with_columns(pl.lit(None).alias(col))
            if col not in dropout_df.columns:
                dropout_df = dropout_df.with_columns(pl.lit(None).alias(col))

        # Reorder both to match
        complete_df = complete_df.select(all_cols)
        dropout_df = dropout_df.select(all_cols)

        # Combine
        augmented_df = pl.concat([complete_df, dropout_df])
    else:
        augmented_df = complete_df

    augmented_df = augmented_df.sort("id")

    print(f"  Complete runs: {complete_df.height}")
    print(f"  Dropout runs with stats: {len(dropout_stats)}")
    print(f"  Total augmented runs: {augmented_df.height}")

    # Save
    print(f"\nSaving to {OUTPUT_PATH}...")
    augmented_df.write_csv(OUTPUT_PATH)
    print(f"  Saved {augmented_df.height} runs")

    # Show sample
    print("\nSample dropout runs with enhanced stats:")
    if dropout_df is not None and not dropout_df.is_empty():
        sample = dropout_df.filter(
            pl.col("kallisto_aligned").is_not_null() |
            pl.col("hisat2_aligned").is_not_null()
        ).select([
            "id", "starting_reads", "fastp_reads_after",
            "kallisto_aligned", "hisat2_aligned", "star_aligned_unique"
        ]).head(5)
        print(sample)

    print("\n" + "="*70)
    print("COMPLETE!")
    print(f"Output: {OUTPUT_PATH}")


if __name__ == "__main__":
    main()
