#!/usr/bin/env python3
"""
Same as virus_svm_prediction_gene_multiclass_hard_random_top6.py, but intended for
gene order CSVs from random_shuffle/exclude_top50/ (shuffled genes after excluding
top-50 SVM-selected genes). Takes the first M genes from that CSV to match selection run M.

Default output directory: gene_selection_final/result_random_exclude_top50/
"""
import argparse
import os
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp

from sklearn.model_selection import train_test_split
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.svm import LinearSVC, SVC
from sklearn.metrics import confusion_matrix


def _default_result_dir() -> str:
    return str(Path(__file__).resolve().parent / "result_random_exclude_top50")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Multiclass SVMs on random gene order from exclude_top50 shuffles, "
        "matching gene count M from svm_gene_selection_top6_n{N}_genes.csv."
    )
    parser.add_argument(
        "--adata",
        required=True,
        help="Path to input .h5ad file.",
    )
    parser.add_argument(
        "--label-col",
        default="virus_classification",
        help="Column in adata.obs to use as ground-truth labels.",
    )
    parser.add_argument(
        "--selection-result-dir",
        required=True,
        help="Directory containing svm_gene_selection_top6_n{N}_genes.csv files.",
    )
    parser.add_argument(
        "--selection-prefix",
        default="svm_gene_selection_top6",
        help="Prefix used by the selection run for gene list files.",
    )
    parser.add_argument(
        "--shuffle-csv",
        required=True,
        help="Path to random_shuffle/exclude_top50/random_shuffle_seed*.csv (column 'gene').",
    )
    parser.add_argument(
        "--result-dir",
        default=_default_result_dir(),
        help="Directory to write outputs (default: .../gene_selection_final/result_random_exclude_top50).",
    )
    parser.add_argument(
        "--n-values",
        default="1,3,5,10,20,50,100,200,500",
        help="Comma-separated list of N values (matching the selection runs).",
    )
    parser.add_argument(
        "--test-size",
        type=float,
        default=0.2,
        help="Test split fraction.",
    )
    parser.add_argument(
        "--random-state",
        type=int,
        default=42,
        help="Random seed.",
    )
    parser.add_argument(
        "--metrics-prefix",
        default="top_n_random_exclude_top50",
        help="Prefix for saved outputs.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    os.makedirs(args.result_dir, exist_ok=True)

    n_values = [int(x) for x in args.n_values.split(",") if x.strip()]

    adata = sc.read_h5ad(args.adata)

    if args.label_col not in adata.obs.columns:
        raise KeyError(f"Label column not found in adata.obs: {args.label_col!r}")

    y = adata.obs[args.label_col].astype(str).fillna("no_infection")
    y = pd.Series(y.to_numpy(), index=adata.obs_names, name="label")

    class_names = sorted(pd.unique(y))
    if "no_infection" in class_names:
        class_names = ["no_infection"] + [c for c in class_names if c != "no_infection"]

    available_genes = set(adata.var_names)

    seed_name = os.path.splitext(os.path.basename(args.shuffle_csv))[0]
    shuffle_genes = pd.read_csv(args.shuffle_csv)["gene"].tolist()
    print(f"Shuffle file (exclude_top50): {args.shuffle_csv} ({seed_name}, {len(shuffle_genes)} genes)")

    for n in n_values:
        sel_genes_path = os.path.join(
            args.selection_result_dir,
            f"{args.selection_prefix}_n{n}_genes.csv",
        )
        if not os.path.exists(sel_genes_path):
            print(f"N={n} selection gene list not found: {sel_genes_path}, skipping.")
            continue

        sel_genes = pd.read_csv(sel_genes_path)["gene"].tolist()
        m = len(sel_genes)
        print(f"N={n}: selection used M={m} genes")

        genes_present = [g for g in shuffle_genes if g in available_genes][:m]
        if not genes_present:
            print(f"N={n} {seed_name}: no valid genes, skipping.")
            continue

        X = adata[:, genes_present].X
        if sp.issparse(X):
            X = X.tocsr()
            scaler = StandardScaler(with_mean=False)
        else:
            X = np.asarray(X)
            scaler = StandardScaler(with_mean=True)

        X_train, X_test, y_train, y_test, idx_train, idx_test = train_test_split(
            X,
            y,
            y.index.to_numpy(),
            test_size=args.test_size,
            random_state=args.random_state,
            stratify=y,
        )

        tag = f"{args.metrics_prefix}_{seed_name}_n{n}"

        genes_path = os.path.join(args.result_dir, f"{tag}_genes.csv")
        pd.Series(genes_present, name="gene").to_csv(genes_path, index=False)

        linear_pipeline = make_pipeline(
            scaler,
            LinearSVC(
                class_weight="balanced",
                random_state=args.random_state,
            ),
        )
        linear_pipeline.fit(X_train, y_train)
        y_pred_lin = linear_pipeline.predict(X_test)
        y_score = linear_pipeline.decision_function(X_test)

        cm = confusion_matrix(y_test, y_pred_lin, labels=class_names)
        cm_df = pd.DataFrame(cm, index=class_names, columns=class_names)
        cm_df.to_csv(os.path.join(args.result_dir, f"{tag}_linear_confusion_matrix.csv"))

        svc_lin = linear_pipeline.named_steps["linearsvc"]
        classes_lin = list(svc_lin.classes_)
        if y_score.ndim == 1:
            y_score = y_score[:, None]
        scores_df = pd.DataFrame(y_score, columns=classes_lin, index=idx_test)
        scores_df.insert(0, "y_true", y_test.to_numpy())
        scores_df.insert(1, "y_pred", y_pred_lin)
        scores_df.to_csv(os.path.join(args.result_dir, f"{tag}_linear_scores.csv"))

        rbf_pipeline = make_pipeline(
            scaler,
            SVC(
                kernel="rbf",
                probability=True,
                class_weight="balanced",
                random_state=args.random_state,
            ),
        )
        rbf_pipeline.fit(X_train, y_train)
        y_pred_rbf = rbf_pipeline.predict(X_test)
        y_proba = rbf_pipeline.predict_proba(X_test)

        cm = confusion_matrix(y_test, y_pred_rbf, labels=class_names)
        cm_df = pd.DataFrame(cm, index=class_names, columns=class_names)
        cm_df.to_csv(os.path.join(args.result_dir, f"{tag}_rbf_confusion_matrix.csv"))

        svc_rbf = rbf_pipeline.named_steps["svc"]
        classes_rbf = list(svc_rbf.classes_)
        proba_df = pd.DataFrame(y_proba, columns=classes_rbf, index=idx_test)
        proba_df.insert(0, "y_true", y_test.to_numpy())
        proba_df.insert(1, "y_pred", y_pred_rbf)
        proba_df.to_csv(os.path.join(args.result_dir, f"{tag}_rbf_scores.csv"))

        print(f"  {tag} ({len(genes_present)} genes) done")


if __name__ == "__main__":
    main()
