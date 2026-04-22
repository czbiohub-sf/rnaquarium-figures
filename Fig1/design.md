# RNAquarium "75k" Figure 1

## Panel A: Pipeline Overview (simplified flow)
 - (from external SVG source)

## Panel B: Key metrics [CURRENT]
 - table of input/output metrics
   - table 1: host filter: run & read filtering
   - table 2: metatranscriptome: contig breakdown
 - pie chart of sequencing technologies / seq-detective filtering outcomes
   - pie 1: sequencing technology categories (simple)
   - pie 2: multi-level pie chart expands non-"Other" technologies from pie 1, with outcome ring
 - scatter/hexbin of mapped reads (mate1 vs mate2) colored by seq-detective per-mate classification
 - [moved to pipeline diagram] computational resource usage per step
 - [abandoned] distributions of gene sparsity (min sparsity vs mate1/mate2 mapping ratio)

## Panel C: Transcriptome clustering (formerly Panel D)
 - umap: development stage
 - umap: tissue type
 - legends drawn in a separate step (see below)

## Supplemental: Metadata treemaps (formerly Panel C)
 - tree map of datasets by developmental stage
 - tree map of datasets by tissue type

## Layout
total panel area dimensions; plots may be smaller
all units in cm unless specified
```
[A 16.0 x 5.22]
[B1 6.85 x 3.53][C2 devstage umap 6.44 x 7.125][C2 devstage legend 3.0 x 7.125]
[B2 6.85 x 5.00][C2 tissue umap 6.44 x 7.125][C2 tissue legend 3.0 x 7.125]
[B2 6.85 x 5.75]
```
 - font: Arial
 - minimum font size: 5pt
 - maximum font size: 7pt
 - B2 nests two pies; smaller has radius ~=0.9cm, larger has total radius ~=1.8cm (including outer ring)
 - generic table row top/bottom padding at most 2pt
 - top (5.2cm A) devoted to pipeline diagram
 - below left column (6.85cm B) devoted to run statistics and seq-detective
 - below right column (9.44cm C) devoted to transcriptome output

## text covering
 - optimization at multiple levels 1) ordering of tools matters (ex. efficient aligners first; k-mer tools before slow blast 2) resource profiling and job dispatch 3) smart parameterization + bug fixes to boost compute efficiency
 - developed Seq-detective to handle data format/tech-platform heterogeneity
 - progress monitoring, resume & failure recovery
 - part 1 processing overview, metrics and outputs
 - part 2 processing overview, metrics and outputs

# general rules
 - polars for data processing
 - tables made with `great_tables`, export to HTML, which will be screenshotted, embedded, or
   converted to SVG later.
 - most plots export as SVG for manual editing; all SVG text is rendered as `<text>` elements (`svg.fonttype = 'none'`), never as paths
 - pie charts drawn as polar bar plots to allow multi-level information
 - color palettes: batlow, palette data in `palette/batlow`
   - batlow (main sequential gradient)
   - batlow/CategoricalPalettes/
   - batlow/DiscretePalettes/
 - distribution plots need to draw all samples, so clustering algorithm must be suitable for
   thousand distributions + 77k samples.
 - we expect ~8 meaningful clusters
 - main figures scripts should be "ipynb notebook-like" i.e. with analogies to notebook block grouping and
   heading comments for each block.
 - avoid LLM-isms.  don't over-comment.  don't excessively print status (important checkpoints, statistics okay)
 - **FILE NAMING**: All output files use panel-specific prefixes (e.g., `Fig1_B_1_*`, `Fig1_B_2_*`, etc.)
 - **CONTOUR PLOTS**: Draw single contour level per cluster (`levels=1`) for clarity
 - **ANNOTATION TEXT**: All annotation text must have `horizontalalignment='center'`
 - **BOUNDING DIMENSIONS**: Panel dimensions are exact — draw all the way to the boundary with no added margins

## Panel B: Table: pipeline filtering statistics
status: COMPLETE
old notebook: `filtering.ipynb`
new script path: `Fig1_B_1_filtering.py`

**Data source**: `data/75k_unstable/stats-with-dropouts-enhanced.csv`
- Use OUTPUT columns for each pipeline stage (reads that PASSED the filter)
- For aligners (kallisto, hisat2, star, bowtie2): use `*_unaligned` columns
  - `*_unaligned` = reads that did NOT align to host = passed through to next stage
  - This matches `host-filtering.summary.after-recovery.txt` calculation method
- Counts ALL runs with data at each stage (including partial/dropout runs)

**Stage column mapping**:
```
starting      → starting_reads
seq_detective → fastp_reads_before  (after download)
fastp         → fastp_reads_after   (after quality filtering)
kb_negative   → kallisto_unaligned  (after host filtering)
hisat2        → hisat2_unaligned    (after hisat2 alignment)
star          → star_unaligned      (after star alignment)
bowtie2       → bowtie2_unaligned   (after bowtie2 alignment)
dedup         → dedup_reads_after   (after deduplication)
final         → final_reads         (final unmapped reads)
```

**Outputs**:
- `Fig1_B_1_filtering_stats_full.html`
- `Fig1_B_1_filtering_stats_totals.html`
- `Fig1_B_1_contigs_breakdown.html`

## Panel B: Plots: SRA sequencing technologies pie charts
status: REFINING
old notebook: `create_technology_plots.py`
new script path: `Fig1_B_2_technologies.py`

 - pie chart of sequencing technologies
   - pie 1 (simple): all technology categories including "Other" — no outcome ring
   - pie 2 (multilevel): named technologies only — **"Other" excluded** (this chart shows detail)
     - (inner pie) technology label; proportions recomputed over non-Other samples
     - (outer ring) seq-detective filtering outcome

**Color consistency rules**:
- Tech category colors: shared `TECH_COLOR_MAP` computed once from full batlow palette
  across all categories sorted by count — identical in both charts.
- Outcome colors: aligned with Fig1_B_3 seq-detective palette (see below).
  Single-end outcomes (B / T) are close to but visually distinct from BB / TT.
  ```
  BB → #99882c  (olive,         = B-B in Fig1_B_3)
  TB → #426f52  (forest green,  = T-B in Fig1_B_3)
  BT → #0b2c5c  (dark navy,     = B-T in Fig1_B_3)
  TT → #f29d6c  (coral/peach,   = T-T in Fig1_B_3)
  B  → #c4b030  (lighter olive, close to BB, distinguishable)
  T  → #d4603c  (darker coral,  close to TT, distinguishable)
  ```

**Outputs**:
- `Fig1_B_2_technology_outcomes_multilevel.svg`
- `Fig1_B_2_technology_categories_simple.svg`

## Panel B: Read/mate distributions (mapping)
status: WORKING (seq-detective view implemented)
old notebook: `create_contour_plots.py`

### Approach 1: Clustering view (exploratory)
script: `Fig1_B_3_distribution_mapping.py`

 - distributions of mapped reads (mate1 vs mate2) (seq-detective metric)
   - hexbin plots, clustering by bioproject distributions.
   - must draw all samples.  log scale density on single plot with outlines and annotation of the
     clusters.
     - separate cluster plots as supplemental figure
   - **single contour level** per cluster for clarity

**Data source**: `data/75k_unstable/seqdetective_metrics.parquet`
- Uses Wasserstein distance for bioproject clustering
- Hierarchical clustering with ward linkage
- Cuts to ~8 clusters as expected

**Outputs**:
- `Fig1_B_3_mapping_mate1_vs_mate2_hexbin.svg` (main figure)
- `Fig1_B_3_mapping_dendrogram.svg` (supplemental)
- `Fig1_B_3_mapping_cluster_panels.svg` (supplemental)

### Approach 2: Seq-Detective per-mate classification view (QUESTIONING)
script: `Fig1_B_3_distribution_mapping_seqdetective.py`

Purpose: Show "this is the shape of the input, which is why we need seq-detective to filter it
before RNA-seq processing"

 - distributions of mapped reads (mate1 vs mate2) colored by seq-detective **per-mate** classification
 - seq-detective classifies each mate independently as Biological (B) or Technical (T)
 - scatter/hexbin plots showing all ~48k PE samples in 4 categories:
   - **B-B** (both biological): #6F7845 (olive) - 32,596 samples
   - **T-B** (M1 technical, M2 biological): #D89E50 (orange) - 13,386 samples
   - **B-T** (M1 biological, M2 technical): #36535F (blue-gray) - 873 samples
   - **T-T** (both technical): #F6A986 (salmon) - 1,502 samples
 - supplemental: 2x2 grid of hexbin density plots (one per category)
 - supplemental: breakdown by filtering reason for technical samples

**Styling (Nature-appropriate)**:
- Colors: Custom palette (authoritative source: `BATLOW_COLORS` dict in script)
  - B-B: #99882c (olive/yellow-green, both biological)
  - B-T: #0b2c5c (dark navy, M1 bio M2 tech, minor issue)
  - T-B: #426f52 (forest green, M1 tech M2 bio, common issue)
  - T-T: #f29d6c (coral/peach, both technical)
- Consistent dimensions: Main (8×8), Grid (12×12), Panel (16×10)
- Consistent hexbin gridsize: 50
- Shared scales: all axes 0-1
- Equal aspect ratios: all subplots use `ax.set_aspect('equal')` for consistent hexbin sizing
- Font sizes: Title 6pt, Labels 6pt, Legend 6pt, Ticks 5pt
- Clean styling: thin spines (0.5pt), subtle grid (α=0.15)
- Legend: Short labels ("B B (n=...)"), large markers (2.5×), semitransparent background (α=0.7)
- Main plot title: "Seq-Detective Mate Filtering Outcomes for Mate-Paired Runs"

**Data sources**:
- `data/75k_unstable/seqdetective_metrics.parquet` (mapping rates)
- `data/75k_unstable/seq-detective-judgement-summary-augmented.txt` (per-mate classifications)
  - Format: ID \t mate1_file \t mate2_file \t grade1 \t grade2 \t reason

**Per-mate classification totals** (48,357 PE samples):
- B-B: 32,596 (67.4%) - both mates biological
- T-B: 13,386 (27.7%) - mate1 technical, mate2 biological
- T-T: 1,502 (3.1%) - both mates technical
- B-T: 873 (1.8%) - mate1 biological, mate2 technical

**Common filtering reasons**:
- Biological: "biological fallback assumption", "mate1-mate2 similar by mapping diff"
- Technical: "mate1 technical by mapping diff", "mates < 9% mapping rate", "sc-like readlen"

**Outputs**:
- `Fig1_B_3_mapping_seqdetective_categories.svg` (main: scatter colored by category)
- `Fig1_B_3_mapping_seqdetective_density_grid.svg` (supplemental: 2x2 hexbin grid)
- `Fig1_B_3_mapping_seqdetective_technical_reasons.svg` (supplemental: technical reasons) 

## Computational resource usage (moved to pipeline diagram, Panel A)
status: MOVED — resource counts are now summarized directly in the pipeline diagram SVG (Panel A),
not presented as a separate Panel B figure.
script `Fig1_B_5_resource_usage.py` kept for reference / supplemental use.


## Panel C: Transcriptome UMAPs
status: REFINING
new script path: `Fig1_D_1_umap.py`

reference umap notebook: `/hpc/projects/data.science/duo.peng/RQ_umaps/step2.umap.expression.ipynb`

**Data source**: `data/75k_unstable/75k_unstable_anndata_zfin_aliases_metadata.log2tmmcpm.h5ad`
  (symlink → `/hpc/projects/balla_group/sra_experiments/versioned_zf_output/75k_unstable/adata_obj/75k_unstable_anndata_zfin_aliases_metadata.log2tmmcpm.h5ad`)
- 61,615 samples × 22,252 genes (log2 TMM-CPM normalized)
- Pre-computed `X_pca` and `X_umap` in `obsm` — used directly, no recomputation needed
- If UMAP absent: recomputes following reference notebook (top-2000 HVG, scale, 50 PCA, 100-neighbor graph, min_dist=0.5)

**UMAP styling**:
- No plot border (spines off), no axis labels, no tick marks
- Small inset axes in bottom-left corner with arrow labels: "↑ UMAP 2" (y-axis) and "UMAP 1 →" (x-axis)
  - Inset is purely decorative; drawn as annotation arrows + text, not a real axis

**Color palette for tissue**: batlow discrete (palette/batlow/DiscretePalettes/batlow100.txt)
- "Undetermined" → #9E9E9E (gray)
- "All anatomical structures" → #BDBDBD (light gray)
- All other categories: golden-ratio stride through batlow100

**Color palette for devstage**: discrete chronological gradient
- Chronological order: Embryo → Larval → Juvenile → Adult (4 steps, evenly spaced from batlow sequential)
- Multi-stage → distinct non-gradient color (batlow midpoint, visually separate)
- Undetermined → #9E9E9E (gray)

**Outputs** (no title on any):
- `Fig1_D_1_umap_devstage.svg` (devstage_curation_coarse, 6 categories)
- `Fig1_D_1_umap_tissue.svg` (tissue_curation_coarse, 21 categories)

**Legend step** (drawn as separate figures, placed to the left of UMAPs in layout):
- Shared legend+bar panel per UMAP view (devstage and tissue)
- Each legend entry: color swatch + category label (no n= count in label text)
- Devstage legend order: Embryo, Larval, Juvenile, Adult, Multi-stage, Undetermined (chronological)
- Tissue legend order: by count descending
- Alongside legend: vertical stacked bar showing proportion each category comprises of total dataset
  - n= sample counts used as bar segment labels (not in legend text)
- Outputs:
  - `Fig1_D_1_umap_devstage_legend.svg`
  - `Fig1_D_1_umap_tissue_legend.svg`


## Panel C: UMAPs — data quality / technology views (Fig1_D_2)
status: NEW
new script path: `Fig1_D_2_umap.py`

Reuses pre-computed `X_umap` from D.1 anndata; no recomputation.

**View 1**: colored by seq-detective filtering outcome
- Outcome column: grade1+grade2 for PE (BB/BT/TB/TT), grade1 only for SE (B/T)
- Missing = "Unknown" → gray `#9E9E9E`
- Colors: same fixed palette as `Fig1_B_2 OUTCOME_COLORS` / `Fig1_B_3 BATLOW_COLORS`

**View 2**: colored by sequencing technology
- Same category mapping as Fig1_B_2 (`map_technology_category`)
- "Other" → `#BDBDBD`, "Unknown" → `#9E9E9E`; named techs use batlow by count rank

**Outputs**:
- `Fig1_D_2_umap_sd_outcomes.svg`
- `Fig1_D_2_umap_technology.svg`
