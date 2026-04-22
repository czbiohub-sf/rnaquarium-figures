#!/bin/bash

DATA_FILE="data/75k_unstable/all_zf_datescurated_withGEO.tsv"
OUT_DIR="unique_values"

# Create output directory
mkdir -p "$OUT_DIR"

# Read the header line and get column names
IFS=$'\t' read -r -a headers < <(head -1 "$DATA_FILE")

# For each column
for i in "${!headers[@]}"; do
    col_num=$((i + 1))
    col_name="${headers[$i]}"

    # Extract unique values from column and write to file
    cut -f "${col_num}" "$DATA_FILE" | tail -n +2 | sort -u > "$OUT_DIR/${col_name}.txt"

    echo "Processed column ${col_name} (column ${col_num})"
done

echo "Done! Unique values saved to $OUT_DIR/"
