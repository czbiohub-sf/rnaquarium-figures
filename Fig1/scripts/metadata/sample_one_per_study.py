#!/usr/bin/env python3
"""
Sample one accession per unique study and combine annotations with full metadata.

Reads:
  - annotated_metadata_combined.csv: parsed/annotated metadata
  - full_metadata_combined.csv: raw full text metadata

Outputs:
  - sample_one_per_study_combined.csv: one run per study with both annotations and full metadata
"""

import pandas as pd
import argparse
from pathlib import Path

def main():
    parser = argparse.ArgumentParser(
        description="Sample one run per study and combine annotations with full metadata"
    )
    parser.add_argument(
        "--annotated",
        default="annotated_metadata_combined.csv",
        help="Path to annotated metadata CSV"
    )
    parser.add_argument(
        "--full",
        default="full_metadata_combined.csv",
        help="Path to full metadata CSV"
    )
    parser.add_argument(
        "--output",
        default="sample_one_per_study_combined.csv",
        help="Output file path"
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=None,
        help="Random seed for reproducibility"
    )
    args = parser.parse_args()

    # Load both files
    print(f"Loading annotated metadata from {args.annotated}...")
    annotated = pd.read_csv(args.annotated)

    print(f"Loading full metadata from {args.full}...")
    full = pd.read_csv(args.full)

    # Sample one accession per study (randomly)
    print(f"Sampling one run per study...")
    if args.seed is not None:
        sampled = annotated.groupby("study").sample(n=1, random_state=args.seed)
    else:
        sampled = annotated.groupby("study").sample(n=1)

    print(f"Sampled {len(sampled)} runs from {sampled['study'].nunique()} unique studies")

    # Merge with full metadata
    print(f"Merging with full metadata...")
    merged = sampled.merge(
        full[["accession"] + [col for col in full.columns if col not in sampled.columns]],
        on="accession",
        how="left"
    )

    # Reorder columns: basic metadata, then full text metadata, then annotations
    # Annotation columns (these are the parsed/classified outputs)
    annotation_cols = [
        "bulk", "droncseq", "dropseq", "fluidigm", "indrops", "marsseq", "matqseq",
        "quartzseq", "splitseq", "superseq", "microwellseq", "scirnaseq", "celseq",
        "cytoseq", "seqwell", "strtseq", "sortseq", "icell8", "iclip", "454",
        "scartrace", "nebnext", "nextera", "ribozero", "10x", "generic-scrnaseq",
        "smartseq", "lexogen", "trueseq", "bias_3prime", "bias_5prime",
        "bias_fullength", "stranded", "unstranded", "sel_polya_hint",
        "sel_ribozero_hint", "sel_random_hint", "sel_smallrna_hint", "assay_type",
        "sc_or_bulk_src", "layout", "selection_class_src", "platform_family",
        "instrument_generation", "instrument_model_slug", "generic-scrnaseq-only",
        "technology", "sc_or_bulk", "read_bias", "selection_class", "conflict_flags"
    ]

    # Basic metadata columns (from annotated, not annotations)
    basic_cols = [
        "accession", "experiment", "pool_member", "published_date", "total_spots",
        "semantic_name", "study", "library_strategy", "library_source",
        "library_selection", "library_layout_tag", "platform", "instrument_model"
    ]

    # Full text metadata columns (from full_metadata_combined)
    full_text_cols = [col for col in merged.columns
                      if col not in basic_cols and col not in annotation_cols]

    # Reorder
    reordered_cols = basic_cols + full_text_cols + annotation_cols
    merged = merged[[col for col in reordered_cols if col in merged.columns]]

    # Write output
    print(f"Writing output to {args.output}...")
    merged.to_csv(args.output, index=False)
    print(f"Done! Output has {len(merged)} rows and {len(merged.columns)} columns")

if __name__ == "__main__":
    main()
