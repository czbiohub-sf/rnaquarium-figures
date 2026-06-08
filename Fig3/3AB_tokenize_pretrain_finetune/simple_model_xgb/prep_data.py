import anndata as ad
import scanpy as sc
import numpy as np
import argparse, pickle, os
from pathlib import Path
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument("--input_h5ad", type=str, required=True)
parser.add_argument("--output_dir", type=str, default="splits")
parser.add_argument("--label_key", type=str, default="developmentalStage_anatomyOntologyClass")
parser.add_argument("--n_comps", type=int, default=20)
parser.add_argument("--reps", type=int, default=6)
parser.add_argument("--seed", type=int, default=123)
args = parser.parse_args()

label_key = args.label_key
n_comps = args.n_comps
reps = args.reps

# -----------------
# Load & sanity
# -----------------
adata = ad.read_h5ad(args.input_h5ad)
print(adata)
adata.obs_names_make_unique()
adata.var_names_make_unique()
print("after maing var and obs unique")
print(adata)
print("label key is set to: ", label_key)
print("label key values are: ", adata.obs[label_key].unique())

# If present, drop UMAP to save memory (no error if missing)
adata.obsm.pop("X_umap", None)


# -----------------
# PCA
# -----------------
sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor="cell_ranger")
sc.pp.pca(adata, n_comps=n_comps)
adata = adata[:, adata.var["highly_variable"]].copy()
print("after filtering highly variable genes")
print(adata)

# -----------------
# Label encoding ONCE
# -----------------
if label_key not in adata.obs:
    raise KeyError(f"Label key '{label_key}' not found in adata.obs")

cats = pd.Categorical(adata.obs[label_key])
y_all = cats.codes.astype(np.int32)  # -1 means missing label
int_to_label = list(cats.categories)  # index matches codes >= 0

# -----------------
# Split helper (no obs mutation)
# -----------------
def split_indices(n, test_frac, rng):
    test_size = int(round(n * test_frac))
    if test_size <= 0 or test_size >= n:
        raise ValueError(f"Invalid test size {test_size} for n={n}, frac={test_frac}")
    test_idx = rng.choice(n, size=test_size, replace=False)
    mask = np.zeros(n, dtype=bool)
    mask[test_idx] = True
    return ~mask, mask  # train_mask, test_mask

# -----------------
# Output dir
# -----------------
outdir = Path(args.output_dir)
outdir.mkdir(parents=True, exist_ok=True)

# -----------------
# Generate splits & dumps
# -----------------
X_pca = adata.obsm["X_pca"]  # (n_cells, n_comps)
n = X_pca.shape[0]

for rep in range(1, reps + 1):
    print(f"Replication #{rep}")
    rng = np.random.RandomState(args.seed + rep)  # deterministic but different per rep

    for k in range(5, 100, 5):  # 0.05, 0.10, ..., 0.95
        split = k / 100.0
        print(f"  Split fraction: {split:.2f}")

        train_mask, test_mask = split_indices(n, split, rng)

        X_train = X_pca[train_mask]
        X_test  = X_pca[test_mask]
        y_train = y_all[train_mask]
        y_test  = y_all[test_mask]

        train_data_dict = {
            "X": X_train,
            "y": y_train,
            "int_to_label_dict": {i: lbl for i, lbl in enumerate(int_to_label)}
        }
        eval_data_dict = {
            "X": X_test,
            "y": y_test,
            "int_to_label_dict": {i: lbl for i, lbl in enumerate(int_to_label)}
        }
        print("train data feature count: ", train_data_dict["X"].shape[1])
        print("train data length: ", len(train_data_dict["X"]))
        print("eval data length: ", len(eval_data_dict["X"]))

        with open(outdir / f"train_data_ncomps{n_comps}_testsize{split:.2f}_rep{rep}.pkl", "wb") as f:
            pickle.dump(train_data_dict, f, protocol=pickle.HIGHEST_PROTOCOL)
        with open(outdir / f"eval_data_ncomps{n_comps}_testsize{split:.2f}_rep{rep}.pkl", "wb") as f:
            pickle.dump(eval_data_dict, f, protocol=pickle.HIGHEST_PROTOCOL)


