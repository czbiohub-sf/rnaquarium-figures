#!/usr/bin/env python3
import pandas as pd
from collections import defaultdict
import re

# Read the CSV file
csv_path = 'sample_one_per_study_combined.csv'
df = pd.read_csv(csv_path)

# Limit to first 1000 rows
df = df.head(1000)

discrepancies = defaultdict(list)
conflict_check = defaultdict(lambda: {"flagged": 0, "unflagged": 0})

# Single-cell technology keywords
sc_keywords = {
    'quartz-seq': 'quartzseq',
    'drop-seq': 'dropseq',
    'indrops': 'indrops',
    'mars-seq': 'marsseq',
    'droplet': 'dropseq',
    'smartseq': 'smartseq',
    'smart-seq': 'smartseq',
    'smart-seq2': 'smartseq',
    '10x': '10x',
    'chromium': '10x',
    'sci-rna-seq': 'scirnaseq',
    'ddseq': 'dropseq',
    'fluidigm': 'fluidigm',
    'c1': 'fluidigm',
    'dronc-seq': 'droncseq',
    'split-seq': 'splitseq',
}

selection_keywords = {
    'poly(a)': 'poly_a',
    'polya': 'poly_a',
    'poly-a': 'poly_a',
    'polyA': 'poly_a',
    'ribozero': 'ribozero',
    'ribo-zero': 'ribozero',
    'rrna depletion': 'ribozero',
    'random': 'random_priming',
    'random priming': 'random_priming',
    "3' bias": '3prime',
    "3-prime": '3prime',
    "3 prime": '3prime',
    "5' bias": '5prime',
    "5-prime": '5prime',
}

print("=" * 100)
print("DISCREPANCY ANALYSIS: First 1000 rows")
print("=" * 100)

for idx, row in df.iterrows():
    row_num = idx + 2  # Account for header and 0-indexing
    accession = row['accession']

    # Get metadata fields
    lib_source = str(row['library_source']).upper() if pd.notna(row['library_source']) else ''
    lib_strategy = str(row['library_strategy']).lower() if pd.notna(row['library_strategy']) else ''
    lib_selection = str(row['library_selection']).lower() if pd.notna(row['library_selection']) else ''
    lib_construction = str(row['library_construction_protocol']).lower() if pd.notna(row['library_construction_protocol']) else ''
    title = str(row['title']).lower() if pd.notna(row['title']) else ''
    alias = str(row['alias']).lower() if pd.notna(row['alias']) else ''
    description = str(row['description']).lower() if pd.notna(row['description']) else ''

    # Get annotation fields
    assay_type = str(row['assay_type']).lower() if pd.notna(row['assay_type']) else ''
    sc_or_bulk_src = str(row['sc_or_bulk_src']).lower() if pd.notna(row['sc_or_bulk_src']) else ''
    sc_or_bulk = str(row['sc_or_bulk']).lower() if pd.notna(row['sc_or_bulk']) else ''
    technology = str(row['technology']).lower() if pd.notna(row['technology']) else ''
    selection_class_src = str(row['selection_class_src']).lower() if pd.notna(row['selection_class_src']) else ''
    selection_class = str(row['selection_class']).lower() if pd.notna(row['selection_class']) else ''
    read_bias = str(row['read_bias']).lower() if pd.notna(row['read_bias']) else ''
    conflict_flags = str(row['conflict_flags']).lower() if pd.notna(row['conflict_flags']) else ''

    # Check for single-cell tech keywords in metadata text
    sc_tech_in_text = []
    for keyword, tech_code in sc_keywords.items():
        if keyword in lib_construction or keyword in title or keyword in alias or keyword in description:
            sc_tech_in_text.append((keyword, tech_code))

    # Check 1: SC technology in text but bulk source
    if sc_tech_in_text and lib_source == 'TRANSCRIPTOMIC' and sc_or_bulk == 'bulk':
        has_conflict_flag = 'sc_tech' in conflict_flags if conflict_flags else False
        conflict_check['sc_tech_bulk_mismatch'][("flagged" if has_conflict_flag else "unflagged")] += 1
        if not has_conflict_flag:
            discrepancies['sc_tech_in_text_but_bulk_classification'].append({
                'accession': accession,
                'row': row_num,
                'sc_techs_found': [t[0] for t in sc_tech_in_text],
                'library_source': lib_source,
                'technology_annotated': technology,
                'sc_or_bulk': sc_or_bulk,
                'lib_construction': lib_construction[:100]
            })

    # Check 2: Metatranscriptomic source but rna_seq assay_type
    if lib_source == 'METATRANSCRIPTOMIC' and assay_type == 'rna_seq':
        discrepancies['metatranscriptomic_as_rna_seq'].append({
            'accession': accession,
            'row': row_num,
            'library_source': lib_source,
            'assay_type': assay_type,
            'title': title[:80]
        })

    # Check 3: library_selection field vs selection_class inferred from text
    selection_in_text = []
    for keyword, sel_code in selection_keywords.items():
        if keyword in lib_construction or keyword in description:
            selection_in_text.append((keyword, sel_code))

    # If text mentions selection but library_selection field is unspecified/unknown
    if selection_in_text and lib_selection in ['unspecified', 'unknown', 'other', '']:
        # Check if selection_class was correctly inferred
        inferred_selections = [s[1] for s in selection_in_text]
        if selection_class not in inferred_selections and selection_class != 'unknown':
            discrepancies['selection_discrepancy'].append({
                'accession': accession,
                'row': row_num,
                'library_selection': lib_selection,
                'selection_in_text': [s[0] for s in selection_in_text],
                'selection_class_inferred': selection_class,
                'lib_construction': lib_construction[:100]
            })

    # Check 4: library_selection field mismatch with annotated selection_class_src
    lib_sel_to_annotated = {
        'polya': 'poly_a',
        'poly-a': 'poly_a',
        'polyA': 'poly_a',
        'ribozero': 'ribozero',
        'random': 'random_priming',
        'oligo-dt': 'poly_a',
        'cdna': 'cdna_unspecified',
    }

    lib_sel_normalized = lib_selection.lower().replace('-', '').replace('(', '').replace(')', '')
    for key, expected_val in lib_sel_to_annotated.items():
        if key.replace('-', '') in lib_sel_normalized:
            if selection_class_src != expected_val and selection_class_src != 'unknown':
                discrepancies['library_selection_annotation_mismatch'].append({
                    'accession': accession,
                    'row': row_num,
                    'library_selection_field': lib_selection,
                    'selection_class_src_annotated': selection_class_src,
                    'expected': expected_val
                })
            break

# Print results
print("\n1. SINGLE-CELL TECH IN TEXT BUT BULK CLASSIFICATION (unflagged conflicts)")
print(f"   Count: {len(discrepancies['sc_tech_in_text_but_bulk_classification'])}")
if discrepancies['sc_tech_in_text_but_bulk_classification']:
    for d in discrepancies['sc_tech_in_text_but_bulk_classification'][:20]:
        print(f"   - {d['accession']} (row {d['row']}): {d['sc_techs_found']}")
        print(f"     Protocol: {d['lib_construction']}")

print("\n2. METATRANSCRIPTOMIC SOURCE CLASSIFIED AS RNA_SEQ")
print(f"   Count: {len(discrepancies['metatranscriptomic_as_rna_seq'])}")
if discrepancies['metatranscriptomic_as_rna_seq']:
    for d in discrepancies['metatranscriptomic_as_rna_seq'][:20]:
        print(f"   - {d['accession']} (row {d['row']})")
        print(f"     Title: {d['title']}")

print("\n3. SELECTION CLASS DISCREPANCIES (unspecified field but inferred)")
print(f"   Count: {len(discrepancies['selection_discrepancy'])}")
if discrepancies['selection_discrepancy']:
    for d in discrepancies['selection_discrepancy'][:20]:
        print(f"   - {d['accession']} (row {d['row']}): field='{d['library_selection']}' vs inferred='{d['selection_class_inferred']}'")
        print(f"     Text mentions: {d['selection_in_text']}")

print("\n4. LIBRARY_SELECTION TO ANNOTATION FIELD MISMATCH")
print(f"   Count: {len(discrepancies['library_selection_annotation_mismatch'])}")
if discrepancies['library_selection_annotation_mismatch']:
    for d in discrepancies['library_selection_annotation_mismatch'][:20]:
        print(f"   - {d['accession']} (row {d['row']}): field='{d['library_selection_field']}' → annotated='{d['selection_class_src_annotated']}' (expected {d['expected']})")

# Summary
print("\n" + "=" * 100)
print("SUMMARY")
print("=" * 100)
print(f"Total rows analyzed: {len(df)}")
print(f"Rows with sc-tech in text but bulk classification (unflagged): {len(discrepancies['sc_tech_in_text_but_bulk_classification'])}")
print(f"Rows with metatranscriptomic as rna_seq: {len(discrepancies['metatranscriptomic_as_rna_seq'])}")
print(f"Rows with selection class discrepancies: {len(discrepancies['selection_discrepancy'])}")
print(f"Rows with library_selection annotation mismatch: {len(discrepancies['library_selection_annotation_mismatch'])}")
