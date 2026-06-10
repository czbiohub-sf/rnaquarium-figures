# Figure 2 — Pan-archive atlas of zebrafish transcriptomes + metatranscriptomes

Per `manuscript_main.tex` Figure 2 legend: UMAP of curated developmental
stage + tissue, with accompanying metatranscriptome taxonomic breakdown.

## Scripts moved from Fig1/ (2026-04)

- `umap_devstage_tissue.py` — full UMAP colored by developmental stage and
  by tissue category (formerly `Fig1/Fig1_D_1_umap.py`). Draws the
  large-format transcriptome UMAP plus legend+stacked-bar panels.
- `umap_quality_technology.py` — UMAP colored by Seq-Detective filtering
  outcome and by sequencing technology (formerly `Fig1/Fig1_D_2_umap.py`).

Both scripts reuse the pre-computed `X_umap` in the anndata object, so they
don't recompute the embedding. Both were written expecting to run from
`Fig1/` root (their `data/` and `figures/` paths are still relative to
that). Adjust paths (or `cd ../Fig1 && python ../Fig2/...`) when wiring
them into a Figure-2 regeneration routine.

## Outputs currently in Fig1/

Rendered SVGs from these scripts currently live at
`../Fig1/figures/fig2/umap_*.svg` pending a proper Figure 2 output dir.

## TODO

- Design panels for metatranscriptome taxonomic breakdown (taxa bar plot,
  treemap, heat tree).
- Move UMAP scripts' data/output paths so they resolve relative to `Fig2/`
  rather than `Fig1/`.

## Earlier scope (now superseded)

The original Figure 2 plan described virus/microbe trees — that content has
moved to Figure 4 per the current manuscript draft.
