# Figure 2 - RNAquarium reveals transcriptomic structure across development and tissue, as well as diverse non-host and microbial diversity

Panels A–D are produced by the Jupyter notebook in this directory
(`Fig2ABCD_umap_and_marker_genes.ipynb`). Panels E–G (non-host / microbial
diversity) are generated directly by the RNAquarium pipeline run, except for 2 heat trees generated via `Fig2_FG_heattrees.ibynb`.

## Overview

The notebook builds a UMAP embedding of the ~75k-sample zebrafish RNA-seq atlas
and characterizes its transcriptomic structure across tissues and developmental
stages:

- **UMAP embeddings** of all samples, colored by coarse anatomical tissue
  (`tissue_curation_coarse`) and by developmental stage (`devstage_curation`).
- **Marker-gene discovery** — Wilcoxon rank-sum, one-vs-rest — for each tissue
  and each developmental stage.
- **Literature validation** of the computed markers against curated zebrafish
  tissue- and stage-specific marker dictionaries, cross-referenced with ZFIN
  expression annotations.
- **Marker-gene panels and per-gene expression UMAPs**, including an
  "Adult vs. rest" panel annotated by immune-gene category (MHC, ISG,
  complement, …).

## What's in this directory

| Path | Description |
| --- | --- |
| `Fig2ABCD_umap_and_marker_genes.ipynb` | Main analysis notebook (panels A–D). |
| `tmm.sh` | Wrapper invoking the TMM normalization step on the raw counts matrix. |
| `scripts/tmm.py` | TMM normalization + explicit filtering; emits `log2(TMM-CPM + prior_count)`. |
| `scripts/zebrafish_tissue_markers.py` | Literature-curated tissue marker dictionary (transgenic lines + PMIDs). |
| `scripts/zebrafish_stage_markers.py` | Literature-curated developmental-stage marker dictionary (with PMIDs). |
| `umap_devstage_tissue.py` | Fig 2A (devstage) and 2C (tissue) UMAPs with legends |
| `umap_quality_technology.py` | UMAP colored by Seq-Detective filtering outcome and by sequencing technology |

## Inputs

The notebook and `tmm.sh` are run **from the repository root**, so all paths
below are relative to it. Input data lives under `data/` — use `make data` in
the root `Makefile` to symlink a data source directory.

- `data/75k_anndata_zfin_aliases_metadata.h5ad` — AnnData of the ~75k-sample
  atlas with curated sample metadata (`tissue_curation_coarse`,
  `devstage_curation`, …).
- `output/counts_tmm-norm.parquet` — TMM-normalized counts; intermediate
  produced by the TMM step below.
- `data/drerio_z12-geneid_to_z11-ENSDARG.csv` — gene-ID ↔ Ensembl-ID map.
- `data/immune_genes.csv` — immune-gene categories for the Adult-marker panel.
- `data/marker_gene_expression_manual_data.csv` — manual marker-expression
  annotations (optional; used only to validate the automated ZFIN lookup).
- ZFIN `xpat_fish.txt` — downloaded automatically from zfin.org and cached in
  `output/`.

## Outputs

Written to `output/`:

- UMAP figures (notebook): `umap_tissue_curation_coarse.{pdf,svg}`,
  `umap_devstage_curation.{pdf,svg}`.
- UMAP figures (standalone scripts): `Fig2_A_umap_devstage_curation.svg`,
  `Fig2_A_umap_devstage_curation_legend.svg`,
  `Fig2_C_umap_tissue_curation_coarse.svg`,
  `Fig2_C_umap_tissue_curation_coarse_legend.svg`,
  `umap_sd_outcomes.svg`, `umap_technology.svg`.
- Marker tables: `{tissue_curation_coarse,devstage_curation}_markers.csv`
  (with companion `*_rank_genes_groups.pkl`).
- Marker-gene panels: `*_marker_genes.pdf`,
  `{tissue,devstage}_marker_panels_highlighted.{pdf,png}`,
  `adult_markers_annotated_categories.{pdf,png}`.
- Per-gene expression UMAPs: `umap_expression_*_{lit,novel}_markers.{pdf,png}`.
- Literature-reference tables: `highlighted_*_literature_refs.csv`.
- Cached AnnData: `75k_anndata_zfin_aliases_metadata.log2tmmcpm.h5ad`,
  `adata_umap_processed.h5ad`.

## Usage

Run everything from the repository root.

1. **TMM-normalize the raw counts matrix.** `scripts/tmm.py` filters the matrix
   and writes `log2(TMM-CPM + prior_count)` to a parquet file. See `tmm.sh` for
   the invocation; set the `$PATH_TO_*` variables to point at your data first.

2. **Run the standalone UMAP scripts** (optional; produce panel-labelled SVGs).
   Run from the `Fig2/` directory:
   ```sh
   python umap_devstage_tissue.py
   python umap_quality_technology.py
   ```

3. **Run the notebook.**
   ```sh
   jupyter lab Fig2/Fig2ABCD_umap_and_marker_genes.ipynb
   ```
   Two control flags at the top of the notebook govern recomputation:
   - `FORCE_RECOMPUTE` (default `False`) — recompute from scratch vs. load cached
     artifacts from `output/`.
   - `SKIP_UMAP` (default `True`) — skip the (slow) UMAP recomputation and reuse
     cached coordinates.
