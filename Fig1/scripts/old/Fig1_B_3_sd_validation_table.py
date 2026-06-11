#!/usr/bin/env python3
"""
Figure 1 Panel B.3: Seq-Detective Validation Table

Validates Seq-Detective (SD) performance by checking whether it assigns the
correct B/T outcome given the known technology and file submission format.
Ground truth is based on domain knowledge about library chemistry, not SD output.

Outputs:
  figures/Fig1_B_3_sd_validation_table.html
  figures/Fig1_B_3_10x_file_type_investigation.html
"""

import polars as pl
import pandas as pd
from pathlib import Path
from great_tables import GT

# =============================================================================
# Configuration
# =============================================================================

ANNOTATED_META = Path("scripts/metadata/annotated_metadata_combined.csv")
CURATED_FILE = Path("data/zf-core-v2-74K_problematic_with_bulk.kmers.allcols.csv")
SD_FILE = Path("data/75k_unstable/seq-detective-judgement-summary-augmented.txt")
OUTPUT_DIR = Path("figures")

_GT_OPTIONS = dict(
    table_font_names=["Arial", "Helvetica Neue", "sans-serif"],
    table_font_size="6pt",
    heading_title_font_size="7pt",
    column_labels_font_size="6pt",
    data_row_padding="2pt",
)

# Expected SD outcomes per (tech_canonical, semantic_name_type)
# Based on library chemistry domain knowledge
EXPECTED_OUTCOMES = {
    ("10x",           "fastq x2"):  "TB",   # R1=barcode+UMI (Technical), R2=cDNA (Biological)
    ("10x",           "fastq x3"):  "B",    # SRA serves R1+R2 as single bio stream; I1 discarded
    ("10x",           "BAM"):       "B",    # SRA marks BAM as biological
    ("10x",           "fastq x4+"): "TB",   # typical 10x multi-lane submission
    ("10x",           "fastq x1"):  "B",    # single stream, biological
    ("smartseq",      "fastq x2"):  "BB",   # both mates are full-length cDNA
    ("smartseq",      "fastq x1"):  "B",
    ("bulk_or_other", "fastq x2"):  "BB",
    ("bulk_or_other", "fastq x1"):  "B",
    ("bulk_depleted", "fastq x2"):  "BB",
    ("bulk_depleted", "fastq x1"):  "B",
}

TECH_DISPLAY = {
    "10x":           "10x Chromium",
    "smartseq":      "Smart-seq",
    "bulk_or_other": "Bulk RNA-seq",
    "bulk_depleted": "rRNA-depleted",
}

# =============================================================================
# Helper functions
# =============================================================================

def map_tech_canonical(tech_value: str | None) -> str | None:
    """Map curated technology string to canonical label."""
    if tech_value is None or (isinstance(tech_value, float)):
        return None
    v = str(tech_value).lower()
    if "10x" in v or "chromium" in v:
        return "10x"
    if "smart-seq" in v or "smart_seq" in v:
        return "smartseq"
    # Note: "SMART-seq" and "Smart-seq" both match above via lowercase
    if "detct" in v or "corall" in v or "rrna" in v or "ribodepletion" in v:
        return "bulk_depleted"
    return "bulk_or_other"


def classify_semantic_name(name: str | None) -> str:
    """Simplify SRA file format string to a coarse category."""
    if name is None or str(name).strip() == "" or str(name).strip().upper() == "NA":
        return "unknown"
    name_lower = str(name).lower()
    if "bam" in name_lower:
        return "BAM"
    count_fastq = name_lower.count("fastq")
    if count_fastq == 1:
        return "fastq x1"
    if count_fastq == 2:
        return "fastq x2"
    if count_fastq == 3:
        return "fastq x3"
    if count_fastq >= 4:
        return "fastq x4+"
    return "other"


# =============================================================================
# 1. Load data
# =============================================================================

print("Loading annotated metadata...")
meta = pl.read_csv(ANNOTATED_META, infer_schema_length=10000).select(
    ["accession", "semantic_name", "conflict_flags"]
)
print(f"  Loaded {len(meta):,} rows from annotated metadata")

print("Loading curated annotations...")
curated_pd = pd.read_csv(CURATED_FILE, low_memory=False, usecols=["accession", "technology"])
curated_pd["tech_canonical"] = curated_pd["technology"].apply(map_tech_canonical)
curated = pl.from_pandas(curated_pd[["accession", "tech_canonical"]])
print(f"  Loaded {len(curated):,} rows from curated file")
print(f"  tech_canonical distribution:")
print(curated["tech_canonical"].value_counts(sort=True))

print("Loading Seq-Detective judgements...")
sd = pl.read_csv(
    SD_FILE,
    separator="\t",
    has_header=False,
    new_columns=["id", "file1", "file2", "j1", "j2", "reason"]
).with_columns(
    id=pl.col("id").str.strip_suffix("_subsample"),
    sd_outcome=pl.concat_str(pl.col("j1"), pl.col("j2"), ignore_nulls=True)
).select(["id", "sd_outcome"])
print(f"  Loaded {len(sd):,} SD judgements")
print(f"  SD outcome distribution:")
print(sd["sd_outcome"].value_counts(sort=True))

# =============================================================================
# 2. Join all three datasets
# =============================================================================

print("\nJoining datasets...")
# Start from curated (has our ground truth tech labels)
joined = (
    curated
    .join(meta.rename({"accession": "accession"}), on="accession", how="left")
    .join(sd.rename({"id": "accession"}), on="accession", how="left")
)
print(f"  Joined table: {len(joined):,} rows")

# =============================================================================
# 3. Compute semantic_name_type
# =============================================================================

print("Computing semantic_name_type...")
# Use polars map_elements for the Python classification function
joined = joined.with_columns(
    semantic_name_type=pl.col("semantic_name").map_elements(
        classify_semantic_name, return_dtype=pl.String
    )
)
print("  semantic_name_type distribution:")
print(joined["semantic_name_type"].value_counts(sort=True))

# =============================================================================
# 4. Compute accuracy statistics
# =============================================================================

print("\nComputing accuracy statistics...")

# Only keep rows with curated tech, SD outcome, and non-null semantic type
analysis = joined.filter(
    pl.col("tech_canonical").is_not_null()
    & pl.col("sd_outcome").is_not_null()
    & (pl.col("semantic_name_type") != "unknown")
    & (pl.col("semantic_name_type") != "other")
)
print(f"  Analysis rows (with tech + SD outcome + known format): {len(analysis):,}")

# Add expected outcome column
# Build mapping as lists for polars replace
keys_tech = [k[0] for k in EXPECTED_OUTCOMES]
keys_fmt = [k[1] for k in EXPECTED_OUTCOMES]
vals_exp = list(EXPECTED_OUTCOMES.values())

expected_map = {f"{t}|{f}": e for (t, f), e in EXPECTED_OUTCOMES.items()}
analysis = analysis.with_columns(
    expected=pl.concat_str(
        pl.col("tech_canonical"), pl.lit("|"), pl.col("semantic_name_type")
    ).replace(expected_map)
).filter(
    # Keep only rows where we have a defined expected outcome
    # (i.e. the concatenated key was in EXPECTED_OUTCOMES)
    pl.concat_str(
        pl.col("tech_canonical"), pl.lit("|"), pl.col("semantic_name_type")
    ).is_in(list(expected_map.keys()))
)

# All possible SD outcome columns
ALL_OUTCOMES = ["B", "T", "BB", "BT", "TB", "TT"]

# Group by tech + format
groups = (
    analysis
    .group_by(["tech_canonical", "semantic_name_type", "expected"])
    .agg([
        pl.len().alias("N"),
        *[
            pl.col("sd_outcome").filter(pl.col("sd_outcome") == outcome).len().alias(outcome)
            for outcome in ALL_OUTCOMES
        ],
        pl.col("sd_outcome").filter(
            pl.col("sd_outcome") == pl.col("expected")
        ).len().alias("correct_count")
    ])
    .filter(pl.col("N") >= 10)
    .with_columns(
        accuracy=(pl.col("correct_count") / pl.col("N") * 100).round(1),
        Technology=pl.col("tech_canonical").replace(TECH_DISPLAY),
    )
    .rename({"semantic_name_type": "Format", "expected": "Expected"})
    .sort(["tech_canonical", "N"], descending=[False, True])
)

print("\nAccuracy table:")
print(groups.select(["Technology", "Format", "Expected", "N", "accuracy"] + ALL_OUTCOMES))

# Print per-tech summary
print("\nAccuracy summary by technology:")
for tech in ["10x", "smartseq", "bulk_or_other", "bulk_depleted"]:
    tech_rows = groups.filter(pl.col("tech_canonical") == tech)
    if len(tech_rows) > 0:
        total_n = tech_rows["N"].sum()
        weighted_acc = (tech_rows["N"] * tech_rows["accuracy"]).sum() / total_n
        print(f"  {TECH_DISPLAY.get(tech, tech)}: {weighted_acc:.1f}% weighted accuracy (N={total_n:,})")

# =============================================================================
# 5. Build great_tables HTML — Main Validation Table
# =============================================================================

print("\nBuilding main validation table...")

table_df = groups.select(
    ["Technology", "Format", "Expected", "N", "accuracy"] + ALL_OUTCOMES
).rename({"accuracy": "Accuracy (%)"})

gt_main = (
    GT(table_df)
    .tab_header("Seq-Detective Validation: Accuracy by Technology and File Format")
    .tab_spanner(
        label="SD Outcome Counts",
        columns=ALL_OUTCOMES
    )
    .fmt_integer(columns=["N"] + ALL_OUTCOMES)
    .fmt_number(columns=["Accuracy (%)"], decimals=1)
    .data_color(
        columns=["Accuracy (%)"],
        palette=["#d62728", "#2ca02c"],
        domain=[0, 100]
    )
    .tab_source_note(
        source_note=(
            "Expected outcomes based on library chemistry: "
            "10x barcode reads (R1) are technical, cDNA reads (R2) biological; "
            "Smart-seq produces full-length cDNA in both mates."
        )
    )
    .tab_options(**_GT_OPTIONS)
)

out_main = OUTPUT_DIR / "Fig1_B_3_sd_validation_table.html"
out_main.write_text(gt_main.as_raw_html())
print(f"  Saved main table to {out_main}")

# =============================================================================
# 6. Secondary table: 10x file type investigation
# =============================================================================

print("Building 10x investigation table...")

tenx_df = (
    analysis
    .filter(pl.col("tech_canonical") == "10x")
    .group_by(["semantic_name_type", "sd_outcome"])
    .agg(pl.len().alias("count"))
    .pivot(index="semantic_name_type", on="sd_outcome", values="count")
    .sort("semantic_name_type")
)

# Ensure all outcome columns exist (fill missing with 0)
for col in ALL_OUTCOMES:
    if col not in tenx_df.columns:
        tenx_df = tenx_df.with_columns(pl.lit(0).cast(pl.UInt32).alias(col))

# Add total column
tenx_df = tenx_df.with_columns(
    N=pl.sum_horizontal([pl.col(c) for c in ALL_OUTCOMES])
).rename({"semantic_name_type": "File Format"})

# Fill nulls with 0
for col in ALL_OUTCOMES + ["N"]:
    tenx_df = tenx_df.with_columns(pl.col(col).fill_null(0))

tenx_df = tenx_df.sort("N", descending=True)

print("  10x file type investigation:")
print(tenx_df)

gt_tenx = (
    GT(tenx_df)
    .tab_header("10x Chromium: SD Outcome by SRA File Format")
    .tab_spanner(
        label="SD Outcome Counts",
        columns=ALL_OUTCOMES
    )
    .fmt_integer(columns=["N"] + ALL_OUTCOMES)
    .data_color(
        columns=["TB"],
        palette=["#f7f7f7", "#2ca02c"],
    )
    .tab_source_note(
        source_note=(
            "All 10x runs in curated dataset. "
            "TB = mate1 Technical (barcode/UMI), mate2 Biological (cDNA) — the expected outcome for 10x fastq ×2. "
            "B = single biological stream (expected for BAM, fastq ×1, fastq ×3)."
        )
    )
    .tab_options(**_GT_OPTIONS)
)

out_tenx = OUTPUT_DIR / "Fig1_B_3_10x_file_type_investigation.html"
out_tenx.write_text(gt_tenx.as_raw_html())
print(f"  Saved 10x investigation table to {out_tenx}")

print("\nDone!")
