#!/usr/bin/env python3
"""
Create a table of 10x technology filtering outcomes by semantic_name.

Reads annotated_metadata_combined.csv and seq-detective-judgement-summary-all.txt
to create a cross-tabulation showing how many samples of each outcome
appear in each data format (semantic_name).
"""

import pandas as pd
import numpy as np
from pathlib import Path

# File paths
script_dir = Path(__file__).parent
metadata_file = script_dir / "annotated_metadata_combined.csv"
seqdet_file = script_dir.parent.parent / "data" / "75k_unstable" / "seq-detective-judgement-summary-all.txt"

print(f"Reading metadata from {metadata_file}")
df_meta = pd.read_csv(metadata_file)

print(f"Reading seq-detective results from {seqdet_file}")
# Parse seq-detective file - always 6 columns, but structure differs:
# Single-end: run_id, filename, empty, outcome, empty, description
# Paired-end: run_id, mate1, mate2, outcome1, outcome2, description
# Detect paired-end by checking if column 3 looks like a filename

rows = []
with open(seqdet_file) as f:
    for line in f:
        fields = line.rstrip('\n').split('\t')
        run_id = fields[0]

        # Check if column 3 (index 2) is a second filename (paired-end)
        is_paired = fields[2] and ('fastq' in fields[2].lower() or 'bam' in fields[2].lower())

        if is_paired:
            # Paired-end: outcomes in fields 4 and 5 (indices 3 and 4)
            outcome = fields[3] + fields[4]
        else:
            # Single-end: outcome in field 4 (index 3)
            outcome = fields[3]

        rows.append({'run_id': run_id, 'outcome': outcome})

df_sd = pd.DataFrame(rows)

# Filter for 10x technology
df_10x = df_meta[df_meta['technology'] == '10x'].copy()
print(f"Found {len(df_10x)} 10x samples")

# Strip _subsample suffix from accession to match seq-detective run_id
df_sd['run_id_base'] = df_sd['run_id'].str.replace('_subsample$', '', regex=True)

# Merge seq-detective outcomes with metadata
df_merged = df_10x.merge(
    df_sd[['run_id_base', 'outcome']],
    left_on='accession',
    right_on='run_id_base',
    how='left'
)

print(f"Matched {df_merged['outcome'].notna().sum()} samples with seq-detective results")

# Simplify semantic_name for better table presentation
def simplify_semantic_name(name):
    if pd.isna(name):
        return 'unknown'
    name_str = str(name).strip()
    if 'bam' in name_str.lower():
        return '10x genomics bam'
    # Count fastq occurrences and create concise label
    fastq_count = name_str.lower().count('fastq')
    if fastq_count > 0:
        if fastq_count == 1:
            return 'fastq'
        elif fastq_count == 2:
            return 'fastq ×2'
        elif fastq_count == 3:
            return 'fastq ×3'
        else:
            return f'fastq ×{fastq_count}'
    return name_str

df_merged['semantic_name_simplified'] = df_merged['semantic_name'].apply(simplify_semantic_name)

# Create pivot table: rows = outcome, columns = semantic_name
print("\nCreating pivot table...")
pivot = pd.crosstab(
    df_merged['outcome'],
    df_merged['semantic_name_simplified'],
    margins=True
)

# Reorder outcomes: single-end first (B, T), then paired-end (BB, BT, TB, TT)
outcome_order = ['B', 'T', 'BB', 'BT', 'TB', 'TT', 'All']
outcome_order = [o for o in outcome_order if o in pivot.index]
pivot = pivot.loc[outcome_order]

print("\nPivot table (outcome vs semantic_name):")
print(pivot)

# Save to CSV
output_file = script_dir / "10x_outcome_by_semantic_name.csv"
pivot.to_csv(output_file)
print(f"\nSaved to {output_file}")

# Alternative format: semantic_name with outcome breakdown
print("\n" + "="*60)
print("Alternative format: breakdown by outcome for each semantic_name")
print("="*60)

# Create a more detailed format showing outcomes side-by-side
detailed_pivot = df_merged.groupby(['semantic_name_simplified', 'outcome']).size().unstack(fill_value=0)

# Reorder outcome columns: single-end first (B, T), then paired-end (BB, BT, TB, TT)
outcome_order = ['B', 'T', 'BB', 'BT', 'TB', 'TT']
outcome_order = [o for o in outcome_order if o in detailed_pivot.columns]
detailed_pivot = detailed_pivot[outcome_order]

print(detailed_pivot)

output_file2 = script_dir / "10x_outcome_detailed.csv"
detailed_pivot.to_csv(output_file2)
print(f"\nSaved detailed version to {output_file2}")

# Print summary statistics
print("\n" + "="*60)
print("Summary statistics")
print("="*60)
print(f"Total 10x samples with metadata: {len(df_10x)}")
print(f"Samples with seq-detective outcomes: {df_merged['outcome'].notna().sum()}")
print(f"\nOutcome distribution:")
print(df_merged['outcome'].value_counts(dropna=False))
print(f"\nSemantic name distribution (simplified):")
print(df_merged['semantic_name_simplified'].value_counts())
