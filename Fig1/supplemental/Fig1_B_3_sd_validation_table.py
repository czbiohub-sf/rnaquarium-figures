#!/usr/bin/env python3
"""
Figure 1 Panel B.3: Seq-Detective validation table

Cross-tabulates submitter file format (semantic_name) against SRA-reported
layout (B = single-end, BB = paired-end) and Seq-Detective reclassification,
grouped by metadata-inferred technology.  Exported as one HTML table per
technology group.
"""

import polars as pl
from pathlib import Path
from great_tables import GT, loc, style
import great_tables.loc as gtloc

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

ANNOT_FILE  = Path("scripts/metadata/annotated_metadata_combined.csv")
SD_FILE     = Path("data/75k_unstable/seq-detective-judgement-summary-augmented.txt")
OUTPUT_DIR  = Path("figures")
OUTPUT_DIR.mkdir(exist_ok=True)

# Technologies to emit tables for (in display order)
TECHNOLOGIES = [
    "10x",
    "smartseq",
    "celseq",
    "microwellseq",
    "indrops",
    "scirnaseq",
    "fluidigm",
    "dropseq",
]

# Column groups: (SRA layout, SD judgement, display label)
# SRA layout inferred from whether seq-detective was given 1 or 2 files.
COL_ORDER = [
    ("B",  "B",  "B → B"),
    ("B",  "T",  "B → T"),
    ("BB", "BB", "BB → BB"),
    ("BB", "TB", "BB → TB"),
    ("BB", "BT", "BB → BT"),
    ("BB", "TT", "BB → TT"),
]
COL_KEYS   = [f"{sra}_{sd}" for sra, sd, _ in COL_ORDER]
COL_LABELS = {f"{sra}_{sd}": lbl for sra, sd, lbl in COL_ORDER}

# =============================================================================
# Load & join
# =============================================================================

annot = pl.read_csv(ANNOT_FILE)

sd = (
    pl.read_csv(
        SD_FILE,
        separator="\t",
        has_header=False,
        new_columns=["id", "file1", "file2", "judgement1", "judgement2", "reason"],
    )
    .with_columns(
        id=pl.col("id").str.strip_suffix("_subsample"),
        sd_judgement=pl.concat_str(
            pl.col("judgement1"), pl.col("judgement2"), ignore_nulls=True
        ),
        sra_layout=pl.when(
            pl.col("file2").is_null() | (pl.col("file2") == "")
        ).then(pl.lit("B")).otherwise(pl.lit("BB")),
    )
    .select("id", "sd_judgement", "sra_layout")
)

joined = annot.join(sd, left_on="accession", right_on="id", how="inner")

# =============================================================================
# Simplify semantic_name → display format
# =============================================================================

def simplify_name(s: str | None) -> str:
    if s is None or s in ("NA", ""):
        return "unknown/NA"
    tokens = s.split()
    n = len(tokens)
    types = set(tokens)
    if types == {"fastq"}:
        return "fastq" if n == 1 else f"fastq ×{n}"
    if types == {"bam"}:
        return "bam" if n == 1 else f"bam ×{n}"
    if types == {"cram"}:
        return "cram" if n == 1 else f"cram ×{n}"
    if "10X" in s:
        return "10X bam"
    if types == {"bam", "cram"}:
        return "bam+cram"
    if types == {"bam", "fastq"}:
        return "bam+fastq"
    if "SOLiD_native" in types:
        return "SOLiD"
    if "Illumina" in s:
        return "Illumina native"
    if types == {"nanopore"}:
        return "nanopore"
    if "pacbio" in s.lower():
        return "pacbio"
    return s[:25]

joined = joined.with_columns(
    fmt=pl.col("semantic_name").map_elements(simplify_name, return_dtype=pl.String)
)

# =============================================================================
# Build pivot table for one technology
# =============================================================================

def build_pivot(tech: str) -> pl.DataFrame | None:
    sub = joined if tech == "all" else joined.filter(pl.col("technology") == tech)
    if sub.height == 0:
        return None

    pivot = (
        sub
        .with_columns(
            col_key=pl.concat_str(pl.col("sra_layout"), pl.lit("_"), pl.col("sd_judgement"))
        )
        .group_by("fmt", "col_key")
        .len()
        .pivot(index="fmt", on="col_key", values="len")
        .fill_null(0)
    )

    # Ensure all expected columns are present
    for key in COL_KEYS:
        if key not in pivot.columns:
            pivot = pivot.with_columns(pl.lit(0).cast(pl.Int32).alias(key))

    pivot = pivot.with_columns(
        total=pl.sum_horizontal([pl.col(k) for k in COL_KEYS])
    ).sort("total", descending=True)

    # Convert counts → percentage strings (blank for zero)
    for key in COL_KEYS:
        pivot = pivot.with_columns(
            pl.when(pl.col(key) == 0)
            .then(pl.lit(""))
            .otherwise(
                pl.concat_str(
                    (pl.col(key) * 100 / pl.col("total")).cast(pl.Int32).cast(pl.String),
                    pl.lit("%"),
                )
            )
            .alias(key)
        )

    return pivot.select(["fmt"] + COL_KEYS + ["total"])

# =============================================================================
# Render one GT table
# =============================================================================

def make_table(tech: str, pivot: pl.DataFrame) -> GT:
    n_total = joined.height if tech == "all" else joined.filter(pl.col("technology") == tech).height

    col_labels = {"fmt": "Submitter format", "total": "n"} | COL_LABELS

    gt = (
        GT(pivot)
        .tab_header(f"{tech}  (n = {n_total:,})")
        .tab_spanner(label="SRA = B (single-end)", columns=["B_B", "B_T"])
        .tab_spanner(label="SRA = BB (paired-end)", columns=["BB_BB", "BB_TB", "BB_BT", "BB_TT"])
        .cols_label(**col_labels)
        # Right-align percentage columns
        .cols_align(align="right", columns=COL_KEYS)
        .cols_align(align="right", columns="total")
        # Shade non-empty cells to make pattern visible
        .tab_style(
            style=style.fill(color="#dbeafe"),   # light blue
            locations=loc.body(columns="B_B",    rows=pl.col("B_B")    != ""),
        )
        .tab_style(
            style=style.fill(color="#fee2e2"),   # light red (T-filtered)
            locations=loc.body(columns="B_T",    rows=pl.col("B_T")    != ""),
        )
        .tab_style(
            style=style.fill(color="#dcfce7"),   # light green (proper PE)
            locations=loc.body(columns="BB_BB",  rows=pl.col("BB_BB")  != ""),
        )
        .tab_style(
            style=style.fill(color="#fef9c3"),   # light yellow (TB/BT mixed)
            locations=loc.body(columns="BB_TB",  rows=pl.col("BB_TB")  != ""),
        )
        .tab_style(
            style=style.fill(color="#fef9c3"),
            locations=loc.body(columns="BB_BT",  rows=pl.col("BB_BT")  != ""),
        )
        .tab_style(
            style=style.fill(color="#fee2e2"),   # light red (both T-filtered)
            locations=loc.body(columns="BB_TT",  rows=pl.col("BB_TT")  != ""),
        )
        .tab_options(**_GT_OPTIONS)
    )
    return gt

# =============================================================================
# Main
# =============================================================================

for tech in ["all"] + TECHNOLOGIES:
    pivot = build_pivot(tech)
    if pivot is None:
        print(f"  {tech}: no data, skipping")
        continue

    gt = make_table(tech, pivot)
    out_path = OUTPUT_DIR / f"Fig1_B_3_sd_validation_{tech}.html"
    out_path.write_text(gt.as_raw_html())
    print(f"Saved {out_path}  ({pivot.height} format rows)")

print("\nDone!")
