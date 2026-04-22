# RNAquarium "75k" Figure 1 — 2026-04 layout

Reference draft: `Figure1draft-2026-04.png` (200 dpi export of the current
assembly). Panel letters follow the manuscript legend (Figure 1 A/B/C); 
the teaser row between A and B is unlettered.

## Panels

### Panel A: Pipeline overview (top, full width)
- Source: external drawio file post-processed into SVG.
- Artifacts: `figures/Fig1_A_1_pipeline-sized.drawio.svg` (+ `.png` preview,
  `.svg.bak` last edit backup).
- drawio tweak helper: `scripts/pipeline_diagram/modify_drawio.py` (regex
  patcher for HTML content inside `mxCell` values).
- Run-count totals and CPU summaries are embedded directly in the diagram
  SVG, not rendered as separate figure files.

### Teaser row (between A and B, unlettered, 16.1 × 2.46 cm)
Four forward-reference thumbnails for analyses developed in later figures:
  1. ML embedding / foundation model output (Figure 3)
  2. Tiny transcriptome UMAP (Figure 2 full UMAP)
  3. Gene–taxa correlation heatmap (Figure 5)
  4. Virus phylogeny tree (Figure 4)

Only thumbnail 2 is scripted here (`Fig1_teaser_umap.py` →
`figures/Fig1_teaser_umap_{devstage,tissue}.svg`). The other three come from
their respective figure directories and are dropped in at layout time; they
are not regenerated from Fig1/.

### Panel B: Seq-Detective (bottom left, 8.1 × 5.35 cm)
Shows that Seq-Detective is needed to triage technical vs. biological mates
before RNA-seq processing. Three sub-elements, all drawn by two scripts:

**B.1 — Technology and Seq-Detective outcome pies**
- Script: `Fig1_B_1_pies.py` (produces both pies)
- Outputs:
  - `figures/Fig1_B_1_tech_pie.svg` — simple technology pie covering every
    accession (including "Other"); no outer ring.
  - `figures/Fig1_B_1_seqdetective_pie.svg` — multi-level pie on the
    named-tech subset ("Other" excluded). Inner ring is the same technology
    category as the simple pie; outer ring is the Seq-Detective per-accession
    outcome (BB / TB / BT / TT / B / T).
- Tech colors are computed once from the full batlow palette over all
  categories sorted by count, so the two pies share colors for matching slices.

**B.2 — Mate 1 vs mate 2 mapping scatter**
- Script: `Fig1_B_2_seqdetective_scatter.py`
- Output: `figures/Fig1_B_2_seqdetective_scatter.svg`
- PE accessions only; each point colored by per-mate Seq-Detective
  classification (B-B, T-B, B-T, T-T).
- Draws all ~48k PE samples; hexbin gridsize 50, equal aspect.
- Supplemental companions (drawn by the same script):
  `figures/supplemental/Fig1_seqdetective_density_grid.svg` (2×2 per-category
  hexbin) and `figures/supplemental/Fig1_seqdetective_technical_reasons.svg`
  (technical-class reason breakdown).

**Color palette (shared by B.1 outer ring and B.2)**
```
BB → #99882c  (olive, both biological)
TB → #426f52  (forest green, M1 tech / M2 bio — common issue)
BT → #0b2c5c  (dark navy,    M1 bio  / M2 tech — minor issue)
TT → #f29d6c  (coral/peach,  both technical)
B  → #c4b030  (lighter olive, close to BB but distinguishable)
T  → #d4603c  (darker coral,  close to TT but distinguishable)
```

### Panel C: Pipeline processing outcomes (bottom right, 8.1 × 5.34 cm)
Two `great_tables` HTML tables, screenshotted or SVG-converted for layout.

- Script: `Fig1_C_1_tables.py`
- Outputs:
  - `figures/Fig1_C_1_run_read_filtering.html` (full breakdown by SE / PE /
    T-filt)
  - `figures/Fig1_C_1_run_read_totals.html` (totals only — rendered in the
    main figure)
  - `figures/Fig1_C_1_contigs_breakdown.html` (metatranscriptome transcript
    counts — rendered in the main figure)

Columns use OUTPUT-of-stage read counts (reads that survived the filter).
For aligners (kallisto, hisat2, star, bowtie2), `*_unaligned` counts are
read-through-to-next-stage values; this matches
`host-filtering.summary.after-recovery.txt`.

**Stage → data column mapping**
```
starting      → starting_reads
seq_detective → fastp_reads_before   (post-download)
fastp         → fastp_reads_after    (post-quality)
kb_negative   → kallisto_unaligned   (post host filter)
hisat2        → hisat2_unaligned
star          → star_unaligned
bowtie2       → bowtie2_unaligned
dedup         → dedup_reads_after
final         → final_reads          (input to metatranscriptome)
```

## Layout (cm)

```
[ A 16.1 × 5.22                                            ]
[ teaser 16.1 × 2.46 — ML | tiny UMAP | gene-taxa | virus ]
[ B (pie + scatter) 8.1 × 5.35 ][ C (tables) 8.1 × 5.34   ]
```

- Font: Arial; min 5 pt, max 7 pt.
- Table row padding ≤ 2 pt.
- Panel dimensions are exact — no added margins outside the stated bounds.
- Pie charts are drawn as polar bar plots (lets the outer ring carry a
  second categorical variable without stacking).

## Data sources

Current (post-recovery) files. See `.claude/data_sources.md` for discrepancy
counts and the full column schemas.

| Use | Path |
|---|---|
| Per-run pipeline stats (authoritative) | `data/75k_unstable/stats-with-dropouts-enhanced.csv` |
| Pipeline aggregate summary | `data/75k_unstable/host-filtering.summary.after-recovery.txt` |
| Seq-Detective per-run judgements | `data/75k_unstable/seq-detective-judgement-summary-augmented.txt` |
| Seq-Detective per-mate metrics | `data/75k_unstable/seqdetective_metrics.parquet` |
| SRA/GEO metadata (curated) | `/hpc/projects/balla_group/sra_experiments/SRA_metadata/dec2025_75k_submitteradded/all_zf_dates_devstage_tissue_tech_curated.tsv` |
| Transcriptome anndata | `data/75k_unstable/75k_unstable_anndata_zfin_aliases_metadata.log2tmmcpm.h5ad` (symlink to versioned output) |
| Manual technology annotations | `data/zf-core-v2-74K_problematic_with_bulk.kmers.allcols.csv` |
| SRA accession universe | `data/75k_unstable/ZF_SraEsearch-2025-06-22.csv` |
| Nextflow process trace | `data/75k_unstable/trace-merged-dangerously.txt` |

The pipeline recovery scripts that produced `stats-with-dropouts-enhanced.csv`
and the augmented Seq-Detective summary live under `scripts/recovery/` — see
`scripts/recovery/RECOVERY_README.md`.

## Supplemental

Kept in `supplemental/` (scripts) and `figures/supplemental/` (outputs):

- `supplemental/treemaps_metadata.py` — developmental-stage and tissue
  treemaps (was Panel C before 2026-03).
- `supplemental/resource_usage.py` — horizontal CPU-hour bar chart per
  pipeline step (the main-figure version is inline in the Panel A diagram).
- `figures/supplemental/Fig1_seqdetective_density_grid.svg` and
  `Fig1_seqdetective_technical_reasons.svg` from
  `Fig1_B_2_seqdetective_scatter.py`.

## Moved out

- Full UMAPs (devstage, tissue, Seq-Detective outcome, technology) moved to
  Figure 2 — scripts at `../Fig2/umap_devstage_tissue.py` and
  `../Fig2/umap_quality_technology.py`, renders at `figures/fig2/*.svg`.

## Archived / superseded

- `scripts/old/filtering.ipynb`, `technologies.ipynb` — predecessor notebooks
  for Panels C and B.
- `scripts/old/create_contour_plots.py`, `create_technology_plots.py`,
  `create_panel_c.py` — earlier pandas/matplotlib drafts.
- `scripts/old/Fig1_B_3_distribution_mapping_clustering.py` — abandoned
  Wasserstein-clustering view of the mate1/mate2 scatter.
- `scripts/old/Fig1_B_4_distribution_sparsity.py` — abandoned sparsity
  vs. mapping-ratio panel.
- `figures/old/` — stale output SVGs from the above.

## Conventions

- `polars` for tabular wrangling (avoid mixing with pandas in new scripts).
- Tables via `great_tables`, exported as HTML.
- Plots export as SVG; set `plt.rcParams['svg.fonttype'] = 'none'` so all
  text stays as `<text>` elements rather than paths.
- Pie charts drawn as polar bar plots so the outer ring can carry a second
  categorical variable.
- Color palette: `batlow` (main sequential), `batlow/CategoricalPalettes/`,
  `batlow/DiscretePalettes/`. Palette files in `palette/`.
- Distribution plots draw every sample — a clustering algorithm used for
  panel design has to scale to ~thousand distributions × 77k samples. We
  expected ~8 meaningful clusters when that approach was live.
- Main-figure scripts read like notebooks (block comments per section, no
  excessive print noise).
- Annotation text uses `horizontalalignment='center'`.
- Contour plots use a single level per cluster (`levels=1`).
- Panel dimensions are exact — draw to the boundary, no added margins.
- File-naming: main-figure outputs go `Fig1_<panel>_<n>_<description>.svg`;
  supplemental outputs drop the panel letter (`Fig1_<description>.svg`).
  Scripts may be run from the Fig1/ root (they resolve `data/` and `figures/`
  relative to the CWD).

## Open issues

Flagged during the 2026-04 reorg, to revisit:

- `Fig1_B_1_seqdetective_pie.py` still references the old kmers-annotation
  CSV at `data/zf-core-v2-74K_problematic_with_bulk.kmers.allcols.csv`
  (Sep 2025). Worth auditing whether a fresher annotation source exists.
- `supplemental/treemaps_metadata.py` reads the curated metadata tsv via
  a hard-coded absolute HPC path. A local symlink under `data/` would make
  the script portable.
- `scripts/recovery/extract_seqdetective_metrics.py` was built against
  `seq-detective-judgement-summary-all.txt`; confirm it should now read the
  `-augmented` file (or decide that one-shot extraction from the original
  is the intended behaviour).
- `modify_drawio.py` contains an absolute path to the drawio SVG; relocate
  to `scripts/pipeline_diagram/` kept it tidy but the path should be made
  relative.
