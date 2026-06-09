# Figure 3 - AI modeling of archive-wide transcriptomes and analysis of gene representations

This figure trains a transformer model (a [fork](https://huggingface.co/pengxunduo/Geneformer-RNAquarium) of [Geneformer](https://huggingface.co/ctheodoris/Geneformer)) on the RNAquarium archive of zebrafish bulk transcriptomes, then probes what the model learned about gene biology from its learned gene representations (embeddings):

- **A, B** — classification performance: Geneformer pretrained on RNAquarium and fine-tuned to predict developmental stage and tissue, benchmarked against scrambled/randomized-data controls, a Zebrahub-pretrained model, and simpler baselines (PCA + nearest neighbor, XGBoost).
- **C, D** — gene-representation structure: gene embeddings extracted from the fine-tuned model, turned into a similarity network, Leiden-clustered, and annotated by GO / Reactome enrichment (with a co-expression network as baseline).
- **E, F** — "transcription-factor vector arithmetic": composing gene embedding vectors (à la word2vec analogies) to trace developmental trajectories, assessed against random null models.
- **G** — comparison of the RNAquarium-Geneformer gene neighborhoods against ESM2 protein-language-model neighborhoods.

> **Requires the RNAquarium Geneformer fork.** The modeling in this figure is built on a fork of Geneformer, now public at [`pengxunduo/Geneformer-RNAquarium`](https://huggingface.co/pengxunduo/Geneformer-RNAquarium) (referred to as `Geneformer_RQfork` in the scripts and the conda environment). The scripts below import and run that fork's code, so you need it to reproduce the modeling work. They are an archival copy of the workflow as it was run inside the fork, and the subdirectories mirror its `tokenize/`, `pretrain/`, `finetune/`, and `embeddings/` layout. They document the analysis rather than run turnkey from this repo — cluster filesystem paths are redacted to the `/path/to/` placeholder, and most steps are SLURM/GPU jobs. The two master shell scripts (`3AB_.../tokenization_pretraining_fintuning.sh` and `3CD_.../shell_commands_part3_embedding_analysis.sh`) are the best entry points: each lists every command in run order.

## Pipeline overview

The panels run as a dependency chain:

```
3AB  tokenize → pretrain → fine-tune  ─────────────►  fine-tuned Geneformer model
3CD  extract gene embeddings → cosine similarity + co-expression baseline
        → Leiden clustering + GO/Reactome enrichment + Cytoscape networks
        → full cosine-similarity matrices  ──────────┐
3EF  TF vector arithmetic on embeddings + null models │ (consume 3CD embeddings /
3G   Geneformer vs ESM2 neighbor comparison  ◄────────┘  cosine matrices)
```

| Subdirectory | Panels | What it produces |
| --- | --- | --- |
| `3AB_tokenize_pretrain_finetune/` | A, B | Tokenized datasets, pretrained + fine-tuned model checkpoints, hyperopt result tables, classification metrics & confusion matrices. |
| `3CD_embedding_network_GO/` | C, D | Gene embeddings, cosine-similarity matrices, co-expression correlation matrices, Leiden clusters, GO/Reactome enrichment tables & plots, Cytoscape network files. |
| `3EF_TF_vector_arithmetic/` | E, F | Composite/pairwise embedding null distributions, vector-arithmetic similarity rankings, cosine-distribution plots. |
| `3G_ESM2_comparison/` | G | Gene-neighbor overlap (Jaccard) tables between Geneformer and ESM2 embeddings. |

## Shared inputs

- **RNAquarium bulk counts** — the archive-wide zebrafish RNA-seq count matrix (`counts.parquet` / `counts.csv.gz`), used for tokenization/pretraining and for the co-expression baseline.
- **Zebrahub single-cell atlas** — `zf_atlas_full_v4_release.developmentalStage_anatomyOntologyClass.100-3000.h5ad`, the labeled dataset used for fine-tuning and the baseline models (developmental stage + tissue/anatomy labels).
- **Genome annotation** — `GCF_049306965.1_GRCz12tu_genomic.gtf` (gene biotypes for tokenization).
- **Gene-ID map** — `drerio_z12-geneid_to_z11-ENSDARG.csv` (gene symbol ↔ Ensembl `ENSDARG` ID).
- **ESM2 embeddings** (panel G only) — protein-language-model embeddings, precomputed externally.

## Panels in detail

### `3AB_tokenize_pretrain_finetune/` — model training & benchmarking

Master script: `tokenization_pretraining_fintuning.sh`. Steps:

1. **Tokenize** expression matrices into Geneformer rank-value tokens — `tokenize_RNAquarium/` (real data) plus `tokenize_RNAquarium_randomized/` and `tokenize_RNAquarium_scrambled/` (controls), and `tokenize_zebrahub/`. Helpers live in `tokenize_scripts/` (`tokenization.py`, `table_to_loom.py`, `h5ad_to_loom.py`, `nonzero_median_digests.py`, `token_dicts.py`, `compare_tokens.py`, `utils.py`).
2. **Pretrain** Geneformer for 30 epochs with DeepSpeed — `pretrain_RNAquarium/` (real + `RQrandom` + `RQscramble` variants) and `pretrain_Zebrahub/` (`pretrain_geneformer_w_deepspeed.py`).
3. **Hyperparameter-optimize** the fine-tuning — `finetune_hyperparm_optimization/` and `finetune_zebrahub_model_hyperparam_opt/`; `parse_hyperopt.py` aggregates the SLURM trial logs into `hyperopt_results_{RQ,zebrahub}-trained_model.csv` (committed here).
4. **Fine-tune** the final models with the best hyperparameters — `finetune_final_model/` (`cmd_devtissue_prepdata.sh` → `cmd_devtissue_finetuneOnly.sh`, driving `finetune.py`).
5. **Baselines** — `simple_model_xgb/` runs PCA + XGBoost (Ray Tune / Optuna sweep) on the same labels.
6. **Plot** — `plot_scan_results.ipynb` produces the macro-F1-vs-training-size panels; `replot_conf_mat.py` regenerates confusion matrices from saved `*_test_metrics_dict.pkl`.

### `3CD_embedding_network_GO/` — embedding network & enrichment

Master script: `shell_commands_part3_embedding_analysis.sh`. Steps:

1. **Extract gene embeddings** from the fine-tuned checkpoints — `sbatch_job_extract_gene_embeddings.sh` → `extract_gene_embeddings.py` (one CSV of gene × embedding-dim per model/layer, e.g. `second_to_last_layer`).
2. **Co-expression baseline** — `coexpression/tmm.py` normalizes counts to `log2(TMM-CPM + prior)`; `coexpression/leiden.py` builds Pearson/Spearman/robust correlation matrices and Leiden clusters (`batch_leiden.sh` sweeps resolutions as a SLURM array).
3. **Cluster + enrich** — `sbatch_job_cluster_enrichment.sh` / `_reactome.sh` → `gp_enrich_cli.py` (uses `embedding_analysis.py`) Leiden-clusters the embeddings and runs g:Profiler GO/Reactome enrichment per cluster/resolution; `plot_enrichment_summary_from_csv.py` plots term discovery vs. resolution.
4. **Networks for Cytoscape** — `save_graph_for_cytoscape.py`, `save_thresholded_graph_for_cytoscape.py`, `save_multigene_graph_for_cytoscape.py` export top-k / thresholded gene networks (the paper highlights interferon genes, e.g. `ifnphi1`, `stat1b`, `usp18`, `rsad2`).
5. **Per-gene similarity rankings** — `get_similarity_ranking.py` ranks neighbors by embedding cosine and by co-expression; `merge_similarity_ranking_table.py` combines them; `cosine.py` is a small pairwise utility. `cluster_ranking.ipynb` slices the enrichment results interactively.

This subdir also writes the full cosine-similarity matrices (including the ESM2 one) consumed by panels E–G.

### `3EF_TF_vector_arithmetic/` — vector arithmetic & null models

Entry point: `embedding_analysis.ipynb`. The "NMP experiment" section at the bottom of the 3CD master script drives the supporting scripts: starting from a neuromesodermal-progenitor seed (`tbxta` + `sox2`), embedding vectors are composed step-by-step along a mesodermal track (`msgn1`, `meox1`, `myl1`) and a neural track (`pax6b`, `pax6a`) using `get_similarity_ranking_composite.py` (`--embed-combine mean_l2`). `composite_nulls.py` and `pairewise_nulls.py` generate random-gene null distributions, and `plot_cosine_distribution.py` overlays the observed composite similarities on those nulls.

### `3G_ESM2_comparison/` — Geneformer vs. ESM2

Single notebook: `compare_neighbors.ipynb`. It loads the precomputed cosine-similarity matrices for the RNAquarium-Geneformer embeddings, a randomized-embedding control, and ESM2 protein embeddings; for each shared gene it compares the top-50 neighbor sets (Jaccard / overlap coefficient). Most genes show low overlap while specific families (e.g. olfactory receptors, crystallins) agree strongly — quantifying how much the transcriptomic representation reflects protein-sequence similarity. (No GPU; matrices are precomputed in 3CD.)

## Outputs

Bulk intermediates (tokenized datasets, model checkpoints, embedding CSVs, correlation/cosine matrices, enrichment tables, Cytoscape graphs, per-gene rankings) are written under the working directories of the Geneformer fork and are **not** committed here (see `.gitignore`: `scratch`, `intermediate`). The figure-ready artifacts produced by the plotting notebooks/scripts are the panel PDFs/SVGs, plus the committed hyperopt summary tables in `3AB_tokenize_pretrain_finetune/`.

## Compute environment

Most steps are SLURM batch jobs (`sbatch` / `submit.sh`). Pretraining and embedding extraction use NVIDIA GPUs (A100/H100/H200-class, 80 GB) with multi-GPU DeepSpeed; tokenization, clustering, and enrichment run on CPU nodes. The conda/mamba environment is `Geneformer_RQfork` (Python 3.10.16); key libraries include the Geneformer fork + HuggingFace Transformers, scanpy, leidenalg/igraph, g:Profiler, XGBoost + Ray Tune, and (for panel G) precomputed ESM2 embeddings.
