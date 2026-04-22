# Text-Based Sequencing Technology Classification Guide

## What You Have

Two comprehensive systems for classifying sequencing technologies and library layouts from SRA metadata:

### 1. **REGEX_RULES.md** — Reference Encyclopedia
Complete regex patterns organized by technology class. Use this to:
- Understand available classification patterns
- Build custom rules
- Map keywords to ontology levels (assay_type, method, platform, bias, etc.)

### 2. **apply_regex_classification.py** — Working Implementation
Python classifier with:
- `TechClassifier` class with pre-built regex patterns
- Methods for classifying `experiment.library_construction_protocol` text
- Methods for parsing `experiment.spot_descriptor` XML
- Example usage demonstrating real data classification

---

## Quick Start

### Test the Classifier
```bash
python apply_regex_classification.py
```

Output shows real protocol examples classified:
- **SMART-Seq v4**: ✓ Detected as single_cell + bias_fullength
- **MARS-seq**: ✓ Detected as single_cell + 3' bias + poly-A selection
- **Spot descriptors**: ✓ Parsed read structure, technical reads, layout

### Use in Your Workflow

```python
from apply_regex_classification import TechClassifier
import pandas as pd

classifier = TechClassifier()

# Classify a single protocol
protocol_text = "10x Genomics Chromium Single Cell 3' v3.1..."
results = classifier.classify_protocol(protocol_text)
# → {'single_cell': [('sc_10x', pattern, 'high')], ...}

# Parse spot descriptor XML
xml = '<SPOT_DESCRIPTOR>...<READ_CLASS>Technical Read</READ_CLASS>...</SPOT_DESCRIPTOR>'
features = classifier.classify_spot_descriptor(xml)
# → {'spot_length': 150, 'num_reads': 2, 'has_technical_read': True, ...}

# Apply to a full dataframe
df = pd.read_csv('data.tsv', sep='\t')
classifier_results = df['experiment.library_construction_protocol'].apply(
    classifier.classify_protocol
)
```

---

## Key Finding: High-Signal Columns

| Column | Uniqueness | Signal Strength | Use For |
|--------|------------|-----------------|---------|
| **experiment.library_construction_protocol** | 1,663 unique | Very High | Kit detection, bias, selection |
| **experiment.spot_descriptor** | 63 unique | Very High | Read structure, barcode presence |
| **experiment.library_selection** | 17 unique | High | RNA enrichment (poly-A vs random) |
| **experiment.platform** | 11 unique | High | Sequencer platform (Illumina, PacBio, ONT) |
| **experiment.instrument_model** | 46 unique | Medium | Specific machine model |
| **experiment.library_strategy** | 13 unique | High | Assay type (RNA-Seq, ChIP-Seq, ATAC-seq) |
| **experiment.library_source** | 3 unique | High | Bulk vs TRANSCRIPTOMIC SINGLE CELL |

---

## Classification Hierarchy

### Level 1: **Assay Type** (from `experiment.library_strategy`)
```
RNA-Seq (most common)
├─ Bulk RNA-Seq
└─ Single-cell RNA-Seq
ChIP-Seq, ATAC-seq, miRNA-Seq, ...
```

### Level 2: **Method** (from `experiment.library_construction_protocol`)
```
Single-Cell Methods:
├─ Droplet-based: 10x Genomics, Drop-seq
├─ Well-plate barcoding: CEL-Seq, MARS-seq
├─ Full-length capture: SMART-Seq v4, SMARTer
└─ Other: Microfluidics, Split-pool

Bulk Methods:
├─ PolyA-biased: TruSeq, Illumina standard
├─ 3'-biased: Lexogen QuantSeq
├─ Unbiased: NEBNext, BGISeq
└─ Other: Ribozero, rRNA depletion
```

### Level 3: **Technical Bias** (from `experiment.library_construction_protocol` + `experiment.spot_descriptor`)
```
Direction:
├─ 3' bias (oligo-dT capture, QuantSeq)
├─ 5' bias (rare)
└─ Full-length (SMART-Seq, random priming)

Structure:
├─ Paired-end (2 reads from insert)
├─ Single-end (1 read from insert)
└─ Barcode read structure (detected in spot_descriptor)
```

### Level 4: **Selection** (from `experiment.library_selection`)
```
mRNA enrichment:
├─ Poly-A tail: cDNA_oligo_dT, PolyA
└─ Random priming: cDNA_randomPriming, RANDOM PCR

RNA type filtering:
├─ Small RNA: CAGE, miRNA
├─ rRNA depletion: Inverse rRNA, Ribozero
└─ Other: ATAC-seq, ChIP-Seq specific
```

---

## Example Classifications

### **10x Genomics Single-Cell**
```
Protocol text: "...10x Genomics Chromium Single Cell 3' GEM library...Chromium Next GEM..."
    ✓ Assay: Single-cell RNA-Seq
    ✓ Method: Droplet-based (10x Genomics)
    ✓ Bias: 3' (gene expression)
    ✓ Selection: Poly-A (implicit in 10x protocol)
    ✓ Layout: Paired-end (R1 = barcode+UMI, R2 = transcript)

Spot descriptor: <READ_CLASS>Technical Read</READ_CLASS> + <READ_CLASS>Application Read</READ_CLASS>
    ✓ Barcode present: True
    ✓ Read 1: Technical (barcode/UMI)
    ✓ Read 2: Application (transcript)
```

### **SMART-Seq v4 (Bulk)**
```
Protocol text: "...SMART-Seq v4 Ultra Low Input RNA Kit...cDNA synthesis...amplification..."
    ✓ Assay: RNA-Seq (could be bulk or single-cell depending on starting material)
    ✓ Method: Full-length capture
    ✓ Bias: Full-length (full transcript coverage)
    ✓ Selection: Not poly-A specific (captures whole transcriptome)
    ✓ Library: Nextera XT (fragmented after amplification)
```

### **Lexogen QuantSeq (3'-Biased Bulk)**
```
Protocol text: "...Lexogen QuantSeq 3' FWD library kit...100 bp, single end..."
    ✓ Assay: RNA-Seq
    ✓ Method: 3'-biased bulk
    ✓ Bias: 3' Forward (3' end sequencing)
    ✓ Selection: Poly-A (implicit in 3' protocol)
    ✓ Layout: Single-end
```

---

## Regex Pattern Examples

### For Rule-Building

**Detect 10x Chromium:**
```regex
(?i)(10x|10X|chromium).{0,80}(gem|droplet|v[23]|3[′'].*gene)
```
Matches: "10x Genomics Chromium", "Chromium Single Cell 3' v3.1", etc.

**Detect 3' Bias:**
```regex
(?i)(3[′']|three\s*prime|oligo.?dt|quantseq.*3[′'])
```
Matches: "3'", "oligo-dT", "QuantSeq 3' FWD", etc.

**Detect Technical Reads (Spot Descriptor):**
```regex
<READ_CLASS>Technical\s*Read</READ_CLASS>
```
Indicates presence of barcode/UMI read

**Detect Multiple Reads (Paired-End):**
```regex
(?i)(<READ_SPEC>.*){2,}
```
Indicates paired-end or multiplexed structure

---

## Next Steps

1. **Load your full dataset** (currently reading `unique_values/` files)
   ```python
   df = pd.read_csv('data/75k_unstable/all_zf_datescurated_withGEO.tsv', sep='\t')
   ```

2. **Apply classifier** to add classification columns:
   ```python
   classifier = TechClassifier()
   df['protocol_classification'] = df['experiment.library_construction_protocol'].apply(
       classifier.classify_protocol
   )
   df['spot_features'] = df['experiment.spot_descriptor'].apply(
       classifier.classify_spot_descriptor
   )
   ```

3. **Create ontology mapping** from regex matches → standardized terms:
   - `sc_10x` → assay_type="single-cell", method="droplet", platform="10x", ...
   - `bulk_lexogen` → assay_type="bulk", method="3'-biased", selection="poly-a", ...

4. **Validate against known samples** (use a subset where you've manually verified the tech)

5. **Handle edge cases** (overlapping protocols, ambiguous text, manufacturer variations)

---

## Troubleshooting

### Pattern doesn't match expected text
- Protocols have high variance in description (spaces, hyphens, capitalization)
- Use `(?i)` for case-insensitive matching
- Account for spacing: `smart.?seq` matches "smart-seq" and "smartseq"
- Test patterns in Python: `re.search(pattern, text, re.IGNORECASE)`

### Spot descriptor parsing fails
- Some entries are "NA" or malformed XML
- Check for `pd.isna()` before parsing
- Use try/except around XML parsing

### Multiple matches for one sample
- This is expected (a sample can be both "10x" AND "3' bias" AND "paired-end")
- Aggregate by category to avoid double-counting
- Use confidence levels (high/medium/low) to prioritize

---

## Files in This Directory

- **REGEX_RULES.md** — Regex reference encyclopedia (this document)
- **apply_regex_classification.py** — Working classifier implementation
- **unique_values/*** — Unique value files from SRA metadata (source data)
  - `experiment.library_construction_protocol.txt` (1,663 entries)
  - `experiment.spot_descriptor.txt` (63 entries)
  - Other columns for validation

---

## Citation Notes

If you use these classification rules in published work:
- SRA metadata comes from [NCBI SRA](https://www.ncbi.nlm.nih.gov/sra)
- Regex patterns designed for this ZF dataset (may need tuning for other organisms/projects)
- Document any custom rules you add for reproducibility
