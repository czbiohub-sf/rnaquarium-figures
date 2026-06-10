# Gene-expression SVM classification of viral infection in zebrafish RNA-seq

Code accompanying **Figure 5C–G** of the manuscript. Single-cell-style bulk
zebrafish RNA-seq samples are classified as virus-infected vs uninfected, and
as one of 13 specific zebrafish virus species, from per-sample
`log2(TMM-CPM + 1)` gene expression. Two SVM kernels are evaluated (RBF and
linear), in both **binary** (`infected` vs `no_infection`) and **multiclass**
(14 ground-truth labels) settings, and the gene-importance vectors of the
linear models are used to define small panels of top-N predictor genes that
are then re-evaluated with an RBF SVM.

This directory reproduces:

- **Figure 5C** – multiclass RBF confusion matrix (row-normalized %).
- **Figure 5D** – per-class PR / ROC curves for the top 5 viruses, overlaid
  with the binary `infected` vs `no_infection` curve as a baseline.
- **Figure 5E** – gene-importance heatmap combining the **top ±10 genes** from
  the binary LinearSVC with the **top ±5 genes per class** from the multiclass
  LinearSVC, with classes on the y-axis and genes (hierarchically clustered
  within each chunk) on the x-axis.
- **Figure 5F** – per-virus PR curves comparing whole-genome RBF SVM against
  top-N selected gene panels (N ∈ {1, 3, 5, 10, 20, 50, 100, 200, 500}).
- **Figure 5G** – normalized AUPRC (relative to the whole-genome RBF) as a
  function of N for the top 5 viruses, overlaid with a random-gene baseline
  (random panels drawn after excluding the top-50 SVM-selected genes).

---

## Repository contents

```
SVM_classification/
├── README.md
├── requirements.txt
│
├── figure5cd_analyze_multiclass_gene_rbf.ipynb          # Figure 5C, 5D
├── figure5e_plot_gene_importance_linear_final.ipynb     # Figure 5E
├── figure5fg_plot_gene_selection_top6.ipynb             # Figure 5F, 5G
│
├── train_multiclass_rbf_svm.py        # 14-class RBF SVM         -> multiclass_rbf_svm_*.csv
├── train_multiclass_linear_svm.py     # 14-class LinearSVC       -> multiclass_linear_svm_*.csv  (incl. gene_importance.csv)
├── train_binary_rbf_svm.py            # infected vs no_infection -> binary_rbf_svm_*.csv
├── train_binary_linear_svm.py         # infected vs no_infection -> binary_linear_svm_*.csv     (incl. gene_importance.csv)
│
├── train_top_n_selected_rbf_svm.py                  # top-N genes by importance, retrained as RBF
└── train_top_n_random_exclude_top50_rbf_svm.py      # matched random panels (after excluding top-50 SVM genes)
```

All training scripts write CSVs to a single `--result-dir`. The three plotting
notebooks read from that directory; bundled inputs (the AnnData file) live
under `data/`, and figures are written to `plots/`.

---

## Requirements & setup

Python 3.10 (tested with 3.10.16). Install the dependencies into a fresh
environment:

```bash
python -m venv .venv && source .venv/bin/activate     # or use conda
pip install -r requirements.txt
```

---

## Data

The notebooks resolve all paths **relative to this directory**, configured by a
small block of constants at the top of each notebook (`DATA_PATH`, `RESULT_DIR`,
`WHOLE_GENOME_SCORES_PATH`, `BINARY_SCORES_PATH`, `PLOT_DIR`, prefix
constants). The default layout is:

```
SVM_classification/
├── data/
│   └── 75k_unstable_anndata_log2tmmcpm_classification.h5ad   # input
├── results/                                                  # written by train_*.py
│   ├── multiclass_rbf_svm_*.csv
│   ├── multiclass_linear_svm_*.csv
│   ├── binary_rbf_svm_*.csv
│   ├── binary_linear_svm_*.csv
│   ├── top_n_selected_genes_*.csv
│   └── top_n_random_exclude_top50_*.csv
└── plots/                                                    # written by the notebooks
```

| File | Role |
|------|------|
| `data/75k_unstable_anndata_log2tmmcpm_classification.h5ad` | AnnData (~61k samples × 22,252 genes) of `log2(TMM-CPM + 1)` gene expression. `obs["virus_classification"]` is the ground-truth label (`no_infection` or one of 13 virus species). |
| `results/*.csv` | Outputs of the `train_*.py` scripts (see below). |

**Not included — the AnnData input.**
`75k_unstable_anndata_log2tmmcpm_classification.h5ad` is required to retrain
any model. It is available separately: **[TODO: add repository / DOI /
"available from the authors on request"]**. Place it at
`data/75k_unstable_anndata_log2tmmcpm_classification.h5ad`.

---

## How to run

Launch Jupyter **from this directory** so the relative paths resolve. The
typical reproduction order is:

```bash
cd SVM_classification

# 1. Train the four base classifiers (writes results/<prefix>_*.csv each).
python train_multiclass_rbf_svm.py    --adata data/75k_unstable_anndata_log2tmmcpm_classification.h5ad --result-dir results
python train_multiclass_linear_svm.py --adata data/75k_unstable_anndata_log2tmmcpm_classification.h5ad --result-dir results
python train_binary_rbf_svm.py        --adata data/75k_unstable_anndata_log2tmmcpm_classification.h5ad --result-dir results
python train_binary_linear_svm.py     --adata data/75k_unstable_anndata_log2tmmcpm_classification.h5ad --result-dir results

# 2. Train the top-N classifiers (one run, several N values).
python train_top_n_selected_rbf_svm.py \
    --adata data/75k_unstable_anndata_log2tmmcpm_classification.h5ad \
    --gene-importance-csv results/multiclass_linear_svm_gene_importance.csv \
    --result-dir results \
    --n-values "1,3,5,10,20,50,100,200,500"

# 3. (Optional, for the random-baseline curve in Fig 5G) train the random panels.
#    Each invocation processes one seed. Provide a random_shuffle_seedXX.csv that
#    enumerates the gene order to draw from, after excluding the top-50 SVM genes.
for SEED in 42 43 44 45 46 47 48 49 50 51; do
  python train_top_n_random_exclude_top50_rbf_svm.py \
      --adata data/75k_unstable_anndata_log2tmmcpm_classification.h5ad \
      --selection-result-dir results \
      --selection-prefix top_n_selected_genes \
      --shuffle-csv data/random_shuffle/exclude_top50/random_shuffle_seed${SEED}.csv \
      --result-dir results
done

# 4. Open the notebooks and run top to bottom.
jupyter lab            # or: jupyter notebook
```

Notebook ↔ script dependency map:

| Notebook | Reads (under `results/`) |
|----------|---------------------------|
| `figure5cd_analyze_multiclass_gene_rbf.ipynb` | `multiclass_rbf_svm_scores.csv`, `multiclass_rbf_svm_confusion_matrix.csv`, `binary_rbf_svm_scores.csv` |
| `figure5e_plot_gene_importance_linear_final.ipynb` | `multiclass_linear_svm_gene_importance.csv`, `multiclass_linear_svm_scores.csv`, `multiclass_linear_svm_confusion_matrix.csv`, `binary_linear_svm_gene_importance.csv` |
| `figure5fg_plot_gene_selection_top6.ipynb` | `top_n_selected_genes_rbf_n{N}_scores.csv` (and `*_n{N}_genes.csv`), `multiclass_rbf_svm_scores.csv`, `binary_rbf_svm_scores.csv`, `top_n_random_exclude_top50_random_shuffle_seed{S}_n{N}_rbf_scores.csv` (random overlay only) |

Each notebook starts with a small configuration cell that points to these
files via variables; if your folder layout differs, edit those constants and
everything downstream picks it up.

---

## Method summary

- **Features.** Per-cell `log2(TMM-CPM + 1)` gene expression matrix from
  `adata.X`. All 22,252 genes are used as features in the four base models;
  the two top-N models use only the indicated subset.
- **Labels.** `adata.obs["virus_classification"]` provides the 14-class label
  (`no_infection` plus 13 zebrafish virus species). Binary models collapse all
  virus species into a single `infected` label.
- **Train/test split.** Stratified, `test_size=0.2`, `random_state=42`. The
  same split (same seed) is used across all base models so that figures are
  directly comparable.
- **Models.** `sklearn.svm.SVC` (RBF) with `probability=True` and
  `class_weight="balanced"`; `sklearn.svm.LinearSVC` with the same class
  weighting. Features are standardized (`StandardScaler`); the scaler tolerates
  sparse input (`with_mean=False`).
- **Gene importance.** For LinearSVC, `coef_` is reported directly per class
  (multiclass: one-vs-rest with 14 rows; binary: a single row labelled
  `infected_vs_no_infection`).
- **Top-N panels.** From the multiclass linear importance matrix restricted to
  the top-6 most frequent classes, the union of top ±N genes per class defines
  a panel for each `N ∈ {1, 3, 5, 10, 20, 50, 100, 200, 500}`. The panel is
  then retrained as a multiclass RBF SVM and evaluated with the same
  train/test split. The random baseline draws panels of matching size after
  excluding the 50 highest-importance genes.

A full Materials & Methods description accompanies the manuscript.

---

## Software

NumPy, SciPy, pandas, scikit-learn, scanpy / anndata, seaborn, and Matplotlib
(Python 3.10). Exact versions are pinned in `requirements.txt`.
