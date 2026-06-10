# Figure 1 — RNAquarium pipeline

Pipeline-overview figure for the 75k zebrafish-SRA processing paper.

Panels:

- **A** — pipeline diagram (full width, drawio-sourced SVG)
- **teaser** — unlettered strip showing thumbnails that forward-reference
  later figures (ML, UMAP, gene–taxa heatmap, virus tree)
- **B** — Seq-Detective: technology × outcome pie plus mate1-vs-mate2
  scatter (8.1 × 5.35 cm)
- **C** — pipeline processing outcome tables (8.1 × 5.34 cm)

## Directory layout

```
Fig1/
├── Fig1_B_1_pies.py                    # Panel B, tech pie + tech×outcome pie
├── Fig1_B_2_seqdetective_scatter.py    # Panel B, mate1 vs mate2 scatter
├── Fig1_C_1_tables.py                  # Panel C, run/read/contig tables
├── Fig1_teaser_umap.py                 # teaser strip, tiny UMAP
├── supplemental/
│   ├── treemaps_metadata.py            # devstage & tissue treemaps
│   └── resource_usage.py               # CPU-hour bar chart
├── scripts/
│   ├── recovery/                       # one-time dropout-recovery workflow
│   ├── pipeline_diagram/               # drawio SVG post-processors
│   ├── metadata/                       # metadata-annotation outputs
│   ├── metadata-annotate-v2.py         # current annotation driver
│   ├── metadata-annotate.py            # v1 driver (kept for reference)
│   └── old/                            # superseded notebooks + drafts
├── figures/
│   ├── Fig1_A_1_pipeline-sized.drawio.svg (+ .png, .svg.bak)
│   ├── Fig1_B_1_tech_pie.svg
│   ├── Fig1_B_1_seqdetective_pie.svg
│   ├── Fig1_B_2_seqdetective_scatter.svg
│   ├── Fig1_C_1_{run_read_filtering,run_read_totals,contigs_breakdown}.html
│   ├── Fig1_teaser_umap_{devstage,tissue}.svg
│   ├── supplemental/                   # treemaps, resource, extra pies
│   ├── fig2/                           # full UMAPs (moved to Figure 2)
│   └── old/                            # abandoned/superseded outputs
├── data/                               # local copies + symlinks (gitignored)
├── palette/batlow*                     # Crameri batlow palette data
└── RNAquarium-Pipeline-Counts.svg      # earlier diagram source
```

All main-figure scripts expect to be run from the `Fig1/` root — they use
`Path("data/...")` and write into `figures/` (or `figures/supplemental/`)
relative to the current working directory.

Full UMAP scripts live at `../Fig2/umap_devstage_tissue.py` and
`../Fig2/umap_quality_technology.py`.

## Data sources

Authoritative (post-recovery) files:

| Use | Path |
|---|---|
| Per-run pipeline stats | `data/75k_unstable/stats-with-dropouts-enhanced.csv` |
| Pipeline aggregate summary | `data/75k_unstable/host-filtering.summary.after-recovery.txt` |
| Seq-Detective judgements | `data/75k_unstable/seq-detective-judgement-summary-augmented.txt` |
| Seq-Detective per-mate metrics | `data/75k_unstable/seqdetective_metrics.parquet` |
| Transcriptome anndata | `data/75k_unstable/75k_unstable_anndata_zfin_aliases_metadata.log2tmmcpm.h5ad` |
| Curated metadata (devstage, tissue, tech) | `data/metadata/all_zf_dates_devstage_tissue_tech_curated.tsv` |
| SRA accession list | `data/75k_unstable/ZF_SraEsearch-2025-06-22.csv` |
| Manual tech annotations | `data/zf-core-v2-74K_problematic_with_bulk.kmers.allcols.csv` |
| Nextflow process trace | `data/75k_unstable/trace-merged-dangerously.txt` |

## Regenerate panels

```bash
# from Fig1/
python Fig1_B_1_pies.py                  # → figures/Fig1_B_1_tech_pie.svg
                                         #   figures/Fig1_B_1_seqdetective_pie.svg
python Fig1_B_2_seqdetective_scatter.py  # → figures/Fig1_B_2_seqdetective_scatter.svg
python Fig1_C_1_tables.py                # → figures/Fig1_C_1_*.html
python Fig1_teaser_umap.py               # → figures/Fig1_teaser_umap_*.svg

# supplemental
python supplemental/treemaps_metadata.py # → figures/supplemental/Fig1_treemap_*.svg
python supplemental/resource_usage.py    # → figures/supplemental/Fig1_resource_usage*.svg
```

Panel A is assembled by editing the drawio file and re-exporting; the
`scripts/pipeline_diagram/modify_drawio.py` helper patches the HTML-in-XML
content.

## Recovery workflow (one-time)

The post-recovery data files above were produced by
`scripts/recovery/run_recovery_finalization.sh` which chains five steps:
find recoverable files → copy download stats → augment Seq-Detective →
recompute per-stage stats → regenerate aggregate summary. Details in
`scripts/recovery/RECOVERY_README.md`. These scripts are idempotent but the
file-copy step touches NFS and takes hours.

## Conventions

- `polars` for tabular processing.
- Tables via `great_tables`, HTML output.
- SVG plots with `svg.fonttype = 'none'` so text stays selectable.
- Pie charts drawn as polar bar plots to allow a second outer ring.
- Colors from Crameri batlow (`palette/`); seq-detective outcome hexes are
  locked so pie and scatter stay consistent.
- Main-figure outputs: `Fig1_<panel>_<n>_<description>.svg`. Supplemental
  outputs drop the panel letter.
