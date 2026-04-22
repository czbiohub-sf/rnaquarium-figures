# Regex Rules for Sequencing Technology & Library Layout Classification

Generated from free text fields in `unique_values/`

## 1. Library Kit & Protocol Detection
**Source Column:** `experiment.library_construction_protocol`

### Single-Cell Technologies
```regex
# 10x Genomics (Chromium-based)
(?i)(10x|10X|10-x).{0,50}(chromium|gem|droplet)
(?i)(chromium|gem).{0,30}(single\s*cell|scRNA|3')
(?i)(10x|10X).{0,50}(3'|5'|gene\s*expression)

# CELseq (in-tube barcoding)
(?i)(celseq|cel-seq)
(?i)celseq.{0,50}(barcode|pool|plate)

# MARS-seq (in-well barcoding)
(?i)(mars[-_]?seq|mars\s+seq)

# SMART-Seq (full-length capture, single-cell)
(?i)(smart[-_]?seq|SMARTseq|SMART-seq2|SMARTer)
(?i)smart.{0,30}(ultra|low\s*input)

# Drop-seq / Microfluidics
(?i)(drop[-_]?seq|microfluidic)

# Other single-cell methods
(?i)(geminim|split[-_]?pool|sci[-_]?RNA)
```

### Bulk / Tissue-Level Technologies
```regex
# TruSeq (Illumina bulk RNA-seq)
(?i)(trueseq|true-seq|truseq)
(?i)trueseq.{0,30}(stranded|total|mrna|ribo)

# Lexogen QuantSeq (3' bias, 3' end sequencing)
(?i)(lexogen|quantseq|quant[-_]?seq)
(?i)(quantseq|lexogen).{0,30}(3'|5'|fwd|rev)

# NEBNext
(?i)(nebnext|neb[-_]?next)
(?i)nebnext.{0,30}(ultra|directional|rna)

# Nextera (DNA/cDNA fragmentation)
(?i)(nextera|nextera\s*xt)

# BGISeq / DNBSEQ (BGI platform)
(?i)(bgiseq|bgi[-_]?seq|dnbseq)
```

### Library Selection/Enrichment
```regex
# Poly-A selection (mRNA enrichment)
(?i)(poly\s*a|polya|oligo\s*dt|oligo-dt)
(?i)(m?rna.*enriched|polya.*selection)

# Random priming (whole transcriptome)
(?i)(random.*prim|random.*hexamer)

# Ribosomal RNA depletion
(?i)(ribo.*zero|ribozero|rrna.*depl)

# Small RNA selection
(?i)(small\s*rna|mirna|trna)

# Other enrichments
(?i)(cage|atac[-_]?seq|chip[-_]?seq)
```

## 2. Bias & Protocol Direction
**Source Column:** `experiment.library_construction_protocol`

### 3' vs 5' Bias
```regex
# 3' bias indicators (3' end capture, common for scRNA & poly-A selection)
(?i)(3'|3'.*bias|3'.*end|3'.*gene.*expression|3'.*forward|three.*prime)
(?i)(oligo.*dt|poly.*a|quantseq.*3)
(?i)(smart[-_]?seq.*3'|3'.*rna[-_]?seq)

# 5' bias indicators (5' end capture, rarer)
(?i)(5'|5'.*bias|5'.*end|5'.*forward|five.*prime)

# Full-length bias (whole transcript)
(?i)(full[-_]?length|whole.*length|cdna.*length|transcript.*length)

# Not strand-specific
(?i)(unstranded|non[-_]?stranded)

# Strand-specific
(?i)(stranded|strand.*specific|directional)
```

### Read Configuration
```regex
# Paired-end
(?i)(paired?[-_]?end|pair[-_]?end|pe\b|2x[0-9]+)

# Single-end
(?i)(single[-_]?end|se\b|1x[0-9]+)

# Technical reads (adapters, barcodes, UMIs embedded)
(?i)(umi|barcode|adapter|technical\s*read|index)
```

## 3. Spot Descriptor Patterns
**Source Column:** `experiment.spot_descriptor`

### XML Structure Rules
```regex
# Technical reads present (indicates barcodes or adapters)
<READ_CLASS>Technical\s*Read</READ_CLASS>

# Read type configuration
<READ_TYPE>(Forward|Reverse|Adapter)</READ_TYPE>

# Multiple reads (paired, multiplex barcode structure)
(?i)(<READ_SPEC>.*){2,}

# Single-end reads (only one READ_SPEC)
(?i)^<SPOT_DESCRIPTOR>(?:(?!</SPOT_DESCRIPTOR>).)*<READ_SPEC>.*Forward.*</READ_SPEC>(?:(?!</SPOT_DESCRIPTOR>).)*</SPOT_DESCRIPTOR>$

# Spot length (insert size + adapters)
<SPOT_LENGTH>([0-9]+)</SPOT_LENGTH>

# BASE_COORD patterns indicating barcode position
<BASE_COORD>([0-9]+)</BASE_COORD>
```

### Read Configuration from Descriptor
```regex
# UMI/barcode read (short technical read before main insert)
<READ_SPEC>.*<READ_CLASS>Technical.*<BASE_COORD>[1-4]</BASE_COORD>.*</READ_SPEC>
(?i)<READ_SPEC>.*adapter.*</READ_SPEC>

# Long insert with short second read (possibly barcode)
<SPOT_LENGTH>([5-9][0-9]{2,})</SPOT_LENGTH>.*<READ_SPEC>.*<READ_INDEX>1

# Asymmetric paired reads (one long, one short - typical of barcode + insert)
<BASE_COORD>[1-5]</BASE_COORD>.*<READ_INDEX>1.*<BASE_COORD>([6-9][0-9]+|[1-9][0-9]{2,})
```

## 4. Combined Classification Patterns

### 10x Genomics Single-Cell (Droplet-based)
```regex
# High specificity pattern
(?i)10x.*(?:chromium|gem|droplet).*(?:v[23]|v3|3'|gene.*exp)

# Alternative indicators
(?i)chromium.{0,50}(?:10x|single.*cell|3'.*barcode)
```

### CEL-Seq / MARS-Seq (Well-Plate Barcoding)
```regex
(?i)(?:cel-seq|mars-seq).{0,100}(?:plate|well|barcode|384)
```

### SMART-Seq (Full-Length, Low Input)
```regex
(?i)smart[-_]?seq.*(?:ultra|v[234]|low.*input)
```

### Bulk RNA-Seq with 3' Bias
```regex
(?i)(?:lexogen|quantseq).{0,50}(?:3'|single.*end|forward)
```

### Standard Illumina Bulk (TruSeq, NEBNext)
```regex
(?i)(?:trueseq|nebnext).{0,50}(?:stranded|rna.*seq|polya)
```

## 5. Quick Lookup Map

| Pattern | Column | Indicator |
|---------|--------|-----------|
| `10x.*chromium` | library_construction_protocol | Single-cell droplet |
| `celseq\|mars-seq` | library_construction_protocol | Well-plate barcode |
| `smart-seq.*ultra` | library_construction_protocol | Full-length single-cell |
| `lexogen.*quantseq.*3'` | library_construction_protocol | 3' bias bulk |
| `trueseq.*stranded` | library_construction_protocol | Standard bulk RNA-seq |
| `oligo.?dt` | library_construction_protocol | PolyA selection (3' bias) |
| `random.*prim` | library_construction_protocol | Random priming (whole tx) |
| `<READ_CLASS>Technical` | spot_descriptor | Barcode/UMI present |
| Paired `<BASE_COORD>` | spot_descriptor | Paired-end or multiplex |
| `<SPOT_LENGTH>[5-9][0-9]{2}` | spot_descriptor | Long insert (>500 bp) |

## Notes

- Regexes use `(?i)` for case-insensitive matching
- Free text has significant variation (spaces, hyphens, capitalization)
- Some protocols overlap (e.g., SMART-Seq v4 can be used for bulk or single-cell)
- `spot_descriptor` XML is highly structured; extract via tag matching rather than full regex
- `library_construction_protocol` is unstructured; regex captures fragmented keywords
