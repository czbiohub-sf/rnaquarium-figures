#!/usr/bin/env python3
"""
Figure 1 Panel B.3: Seq-Detective validation table

Agreement between Seq-Detective read-type judgements and human-curated
ground truth (b_rtype) across 156 manually annotated accessions spanning
diverse sequencing technologies and file formats.

Simple table: one row per technology, columns = n and agreement %.
"""

import polars as pl
from pathlib import Path
from great_tables import GT, loc, style

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

BENCHMARK_FILE = Path("data/sd_benchmark_set.tsv")
SD_FILE        = Path("data/75k_unstable/seq-detective-judgement-summary-augmented.txt")
OUTPUT_DIR     = Path("figures")
OUTPUT_DIR.mkdir(exist_ok=True)

# Display order for technology rows (unlisted techs appended alphabetically)
TECH_ORDER = [
    "bulk", "10x", "smartseq", "celseq", "indrops", "microwellseq",
    "scirnaseq", "dropseq", "DeTCT", "TCR",
]

# =============================================================================
# Load & join
# =============================================================================

benchmark = pl.read_csv(BENCHMARK_FILE, separator="\t")

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
    )
    .select("id", "sd_judgement")
)

joined = (
    benchmark
    .join(sd, left_on="accession", right_on="id", how="left")
    .with_columns(
        b_rtype_norm=pl.col("b_rtype").str.replace_all(" ", "").fill_null("missing")
    )
    .filter(
        pl.col("technology_actual").is_not_null() & pl.col("b_rtype").is_not_null()
    )
    .with_columns(
        match=(pl.col("sd_judgement") == pl.col("b_rtype_norm"))
    )
)

print(joined["accession","technology_actual","b_rtype","sd_judgement"].write_csv())

# =============================================================================
# Aggregate by technology
# =============================================================================

# TCR-seq uses read layouts that we don't model; exclude from validation.
joined = joined.filter(pl.col("technology_actual") != "TCR")

agg = (
    joined
    .group_by("technology_actual")
    .agg(n=pl.len(), agree=pl.col("match").sum())
)

# Collapse technologies with n<10 into "other"
agg = agg.with_columns(
    technology_actual=pl.when(pl.col("n") < 10)
    .then(pl.lit("other"))
    .otherwise(pl.col("technology_actual"))
).group_by("technology_actual").agg(
    n=pl.col("n").sum(),
    agree=pl.col("agree").sum(),
).with_columns(
    pct=(pl.col("agree") * 100.0 / pl.col("n")).round(1)
)

# Sort by TECH_ORDER, other and Total at end
known = {t: str(i).zfill(3) for i, t in enumerate(TECH_ORDER)}
agg = agg.with_columns(
    _order=pl.col("technology_actual").replace_strict(
        known | {"other": "998"}, default="999"
    )
).sort(["_order", "technology_actual"]).drop("_order")

# Totals row
totals = agg.clear().extend(pl.DataFrame({
    "technology_actual": ["Total"],
    "n":     pl.Series([joined.height], dtype=agg["n"].dtype),
    "agree": pl.Series([int(joined["match"].sum())], dtype=agg["agree"].dtype),
    "pct":   [round(100.0 * joined["match"].sum() / joined.height, 1)],
}))
table_df = pl.concat([agg, totals])

print(table_df)

# =============================================================================
# great_tables render
# =============================================================================

gt = (
    GT(table_df.drop("agree"))
    .tab_header("Seq-Detective read-type validation")
    .cols_label(
        technology_actual="Technology",
        n="n",
        pct="Agreement %",
    )
    .fmt_number(columns="pct", decimals=1)
    .cols_align(align="right", columns=["n", "pct"])
    .data_color(
        columns="pct",
        rows=pl.col("technology_actual") != "Total",
        palette=["#ffffff", "#28A745"],  # white → dark green
        domain=[75.0, 100.0],
    )
    # Bold the totals row
    .tab_style(
        style=style.text(weight="bold"),
        locations=loc.body(rows=pl.col("technology_actual") == "Total"),
    )
    .tab_options(**_GT_OPTIONS)
)

out_path = OUTPUT_DIR / "Fig1_B_3_sd_validation.html"
out_path.write_text(gt.as_raw_html())
print(f"\nSaved {out_path}")
