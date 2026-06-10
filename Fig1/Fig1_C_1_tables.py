#!/usr/bin/env python3
"""
Figure 1 Panel C.1: Pipeline processing outcome tables

Three tables summarizing RNAquarium run and read counts at each pipeline stage,
plus the metatranscriptome contig breakdown. Written to HTML via great_tables;
screenshot/embed for layout.

Breaks down by Seq-Detective technology group (SE, PE, T-filt).
"""

import polars as pl
from pathlib import Path
from great_tables import GT

_GT_OPTIONS = dict(
    table_font_names=["Arial", "Helvetica Neue", "sans-serif"],
    table_font_size="6pt",
    heading_title_font_size="7pt",
    column_labels_font_size="6pt",
    data_row_padding="2pt",
)

# =============================================================================
# Configuration
# =============================================================================

STATS_FILE = Path("data/75k_unstable/stats-with-dropouts-enhanced.csv")
STATS_MERGED = Path("data/75k_unstable/stats-merged.csv")
SEQDETECTIVE_FILE = Path("data/75k_unstable/seq-detective-judgement-summary-augmented.txt")
TRACE_FILE = Path("data/75k_unstable/trace-merged-dangerously.txt")
HOST_SUMMARY = Path("data/75k_unstable/host-filtering.summary.after-recovery.txt")
OUTPUT_DIR = Path("figures")
OUTPUT_DIR.mkdir(exist_ok=True)

# =============================================================================
# Load data
# =============================================================================

# Load seq-detective judgements
sd_judgement = pl.read_csv(
    SEQDETECTIVE_FILE,
    separator="\t",
    has_header=False,
    new_columns=['id', 'file1', 'file2', 'judgement1', 'judgement2', 'reason']
)

# Map seq-detective judgements to technology groups
# B=usable SE, T=filtered SE, BB=usable PE, TT=low-quality PE, BT/TB=mixed quality
judgement_mapping = {
    "B": "SE",
    "T": "SE",
    "BB": "PE",
    "TT": "PE",
    "TB": "T-filt",
    "BT": "T-filt"
}

judgement_map = sd_judgement.with_columns(
    id=pl.col("id").str.strip_suffix("_subsample"),
    tech_group=pl.concat_str(
        pl.col("judgement1"),
        pl.col("judgement2"),
        ignore_nulls=True
    ).replace(judgement_mapping)
).select("id", "tech_group")

print(f"Loaded {len(judgement_map)} seq-detective judgements")
print(judgement_map.group_by("tech_group").len().sort("tech_group"))

# Load stats file with all pipeline metrics
filter_stats = pl.read_csv(STATS_FILE)
print(f"Loaded {len(filter_stats)} pipeline statistics records")

# Load complete runs and trace for dropout analysis
stats_merged = pl.read_csv(STATS_MERGED)
complete_ids = set(stats_merged["id"].to_list())
n_complete = len(complete_ids)

# Load and deduplicate trace (keep best status per tag/process)
if TRACE_FILE.exists():
    trace = pl.read_csv(TRACE_FILE, separator="\t", null_values="-")
    trace_clean = (
        trace
        .filter(pl.col("task_id") != "task_id")
        .with_columns(pl.col("task_id").cast(pl.Int64))
        .with_columns(
            pl.col("status").replace({
                "COMPLETED": 0, "CACHED": 0, "FAILED": 1, "ABORTED": 2
            }).cast(pl.Int32).alias("status_rank")
        )
    )
    deduped_trace = (
        trace_clean
        .sort(["status_rank", "task_id"], descending=[False, True])
        .unique(subset=["tag", "process"], keep="first")
        .filter(pl.col("tag").is_in(set(judgement_map["id"].to_list())))  # exclude strays
    )
else:
    deduped_trace = None

# =============================================================================
# Calculate run counts by stage and technology group
# =============================================================================

# Join stats with technology groups
stats_with_tech = filter_stats.join(judgement_map, on="id", how="inner")

# Define pipeline stages based on column presence
# Each stage maps to the OUTPUT of that stage (reads that PASSED the filter)
# For aligners: *_unaligned = reads that did NOT align = passed through
stage_definitions = {
    "starting": "starting_reads",
    "seq_detective": "fastp_reads_before",  # after download, before fastp
    "fastp": "fastp_reads_after",           # after fastp quality filtering
    "kb_negative": "kallisto_unaligned",    # after kallisto host filtering (unaligned = non-host)
    "hisat2": "hisat2_unaligned",           # after hisat2 (unaligned = non-host)
    "star": "star_unaligned",               # after star (unaligned = non-host)
    "bowtie2": "bowtie2_unaligned",         # after bowtie2 (unaligned = non-host)
    "dedup": "dedup_reads_after",           # after deduplication
    "final": "final_reads"                  # final unmapped reads for metatranscriptome
}

run_counts_data = []
for stage, col in stage_definitions.items():
    stage_counts = (
        stats_with_tech
        .filter(pl.col(col).is_not_null())
        .group_by("tech_group")
        .agg(count=pl.len())
        .with_columns(step=pl.lit(stage))
    )
    run_counts_data.append(stage_counts)

run_counts = (
    pl.concat(run_counts_data)
    .pivot(
        index="step",
        on="tech_group",
        values="count"
    )
    .with_columns(
        Total=pl.sum_horizontal(pl.col("SE"), pl.col("PE"), pl.col("T-filt"))
    )
    # Order stages in pipeline order
    .with_columns(
        step_order=pl.col("step").replace({
            "starting": 1, "seq_detective": 2, "fastp": 3, "kb_negative": 4,
            "hisat2": 5, "star": 6, "bowtie2": 7, "dedup": 8, "final": 9
        })
    )
    .sort("step_order")
    .drop("step_order")
)

print("\nRun counts by stage:")
print(run_counts)

# =============================================================================
# Calculate read counts by stage and technology group
# =============================================================================

read_counts_data = []
for stage, col in stage_definitions.items():
    stage_reads = (
        stats_with_tech
        .filter(pl.col(col).is_not_null())
        .group_by("tech_group")
        .agg(reads=pl.col(col).sum())
        .with_columns(step=pl.lit(stage))
    )
    read_counts_data.append(stage_reads)

read_counts = (
    pl.concat(read_counts_data)
    .pivot(
        index="step",
        on="tech_group",
        values="reads"
    )
    .with_columns(
        Total=pl.sum_horizontal(pl.col("SE"), pl.col("PE"), pl.col("T-filt"))
    )
    # Order stages in pipeline order
    .with_columns(
        step_order=pl.col("step").replace({
            "starting": 1, "seq_detective": 2, "fastp": 3, "kb_negative": 4,
            "hisat2": 5, "star": 6, "bowtie2": 7, "dedup": 8, "final": 9
        })
    )
    .sort("step_order")
    .drop("step_order")
)

print("\nRead counts by stage:")
print(read_counts)

# =============================================================================
# Create combined table with runs and reads
# =============================================================================

RUN_NAMEMAP = {"SE": "se_runs", "PE": "pe_runs", "T-filt": "tfilt_runs", "Total": "total_runs"}
READ_NAMEMAP = {"SE": "se_reads", "PE": "pe_reads", "T-filt": "tfilt_reads", "Total": "total_reads"}

wide_df = (
    run_counts.rename(RUN_NAMEMAP)
    .join(read_counts.rename(READ_NAMEMAP), on="step")
)

print("\nCombined table:")
print(wide_df)

# =============================================================================
# Generate HTML tables with great_tables
# =============================================================================

# Full table with all technology groups
full_table = (
    GT(wide_df)
    .tab_header("Zebrafish RNAquarium Pipeline Filtering Statistics")
    .fmt_integer(columns=pl.selectors.ends_with("runs"))
    .fmt_number(
        columns=pl.selectors.ends_with("reads"),
        compact=True,
        pattern="{x}",
        scale_by=1,
        n_sigfig=3
    )
    .data_color(columns=pl.selectors.ends_with("reads"))
    .tab_options(**_GT_OPTIONS)
)

html_content = full_table.as_raw_html()
(OUTPUT_DIR / "Fig1_C_1_run_read_filtering.html").write_text(html_content)
print(f"\nSaved full table to {OUTPUT_DIR / 'Fig1_C_1_run_read_filtering.html'}")

# Simplified totals table
totals_table = (
    GT(wide_df.select("step", "total_runs", "total_reads"))
    .tab_header("Zebrafish Run & Read Filtering")
    .cols_label(
        step="Pipeline Step",
        total_runs="Runs",
        total_reads="Reads"
    )
    .fmt_integer(columns="total_runs")
    .fmt_number(
        columns="total_reads",
        compact=True,
        pattern="{x}",
        scale_by=1,
        n_sigfig=3
    )
    .data_color(columns="total_reads", palette="Blues", domain=[1e9, 2e12])
    .tab_options(**_GT_OPTIONS)
)

html_content = totals_table.as_raw_html()
(OUTPUT_DIR / "Fig1_C_1_run_read_totals.html").write_text(html_content)
print(f"Saved totals table to {OUTPUT_DIR / 'Fig1_C_1_run_read_totals.html'}")

# =============================================================================
# Optional: Create contig breakdown table (Part II metrics)
# =============================================================================

# These are hardcoded summary statistics from Part II (metatranscriptome assembly)
contigs_df = pl.DataFrame({
    "step": [
        "rnaSPAdes",
        "NT Danio + human filter",
        "NT or NR hits",
        "&nbsp;&nbsp;&nbsp;&nbsp;NT BLAST hits",
        "&nbsp;&nbsp;&nbsp;&nbsp;NR Diamond hits",
        "BBDuk flagged removed",
        "&nbsp;&nbsp;&nbsp;&nbsp;Dark Matter (non-NT/NR)"
    ],
    "contigs": [
        89.8e6,
        48.0e6,
        37.9e6,
        30.8e6,
        28.7e6,
        34.7e6,
        10.1e6
    ]
})

contigs_table = (
    GT(contigs_df)
    .tab_header("Zebrafish Transcript Breakdown")
    .cols_label(
        step="Part II Step",
        contigs="Transcripts"
    )
    .fmt_number(
        columns="contigs",
        compact=True,
        pattern="{x}",
        scale_by=1,
        n_sigfig=3
    )
    .data_color(columns="contigs", palette="Oranges", domain=[1e6, 1e8])
    .tab_options(**_GT_OPTIONS)
)

html_content = contigs_table.as_raw_html()
(OUTPUT_DIR / "Fig1_C_1_contigs_breakdown.html").write_text(html_content)
print(f"Saved contigs table to {OUTPUT_DIR / 'Fig1_C_1_contigs_breakdown.html'}")

# =============================================================================
# Trace-adjusted run counts table (complete + dropout per stage, by tech group)
# =============================================================================

if deduped_trace is not None:
    tech_group_map = dict(zip(judgement_map["id"].to_list(), judgement_map["tech_group"].to_list()))

    complete_tech = {}
    for run_id in complete_ids:
        tg = tech_group_map.get(run_id)
        if tg:
            complete_tech[tg] = complete_tech.get(tg, 0) + 1

    # Each row = output of that step = unique tags seen by the *next* process
    # (regardless of whether the next process completed).
    # read_col = column representing reads output from that step.
    STAGE_NEXT_PROCESS = [
        ("Input",         "download",    "starting_reads"),
        ("Seq-Detective", "fastp",       "fastp_reads_before"),
        ("fastp",         "kb_negative", "fastp_reads_after"),
        ("Kallisto",      "hisat2",      "kallisto_unaligned"),
        ("HISAT2",        "star",        "hisat2_unaligned"),
        ("STAR",          "bowtie2",     "star_unaligned"),
        ("Bowtie2",       "dedup",       "bowtie2_unaligned"),
        ("dedup",         "gsnap",       "dedup_reads_after"),
        ("GSNAP (final)", "stats_csv",   "final_reads"),
    ]

    trace_rows = []
    for display, next_process, read_col in STAGE_NEXT_PROCESS:
        total_runs = deduped_trace.filter(
            pl.col("process") == next_process
        )["tag"].n_unique()

        total_reads = (
            stats_with_tech
            .filter(pl.col(read_col).is_not_null())[read_col]
            .sum()
        ) if read_col in stats_with_tech.columns else None

        trace_rows.append({
            "step": display,
            "total_runs": total_runs,
            "total_reads": total_reads,
        })

    trace_df = pl.DataFrame(trace_rows)

    trace_table = (
        GT(trace_df)
        .tab_header("Zebrafish Run & Read Filtering")
        .cols_label(step="Pipeline Step", total_runs="Runs", total_reads="Reads")
        .fmt_integer(columns="total_runs")
        .fmt_number(
            columns="total_reads",
            compact=True,
            pattern="{x}",
            scale_by=1,
            n_sigfig=3
        )
        .data_color(columns="total_reads", palette="Blues", domain=[1e9, 2e12])
        .tab_options(**_GT_OPTIONS)
    )

    html_content = trace_table.as_raw_html()
    (OUTPUT_DIR / "Fig1_C_1_trace_adjusted.html").write_text(html_content)
    print(f"\nSaved trace-adjusted table to {OUTPUT_DIR / 'Fig1_C_1_trace_adjusted.html'}")
    print("\nTrace-adjusted run counts:")
    print(trace_df)

print("\nDone!")
