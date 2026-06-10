#!/usr/bin/env python3
import argparse
import os

import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp

from sklearn.model_selection import train_test_split
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.svm import LinearSVC, SVC
from sklearn.metrics import confusion_matrix


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Train multiclass SVMs on top-N gene sets from top-6 classes (linear + RBF)."
    )
    parser.add_argument(
        "--adata",
        required=True,
        help="Path to input .h5ad file.",
    )
    parser.add_argument(
        "--gene-importance-csv",
        required=True,
        help="Path to gene importance CSV (rows: classes, columns: genes).",
    )
    parser.add_argument(
        "--label-col",
        default="virus_classification",
        help="Column in adata.obs to use as ground-truth labels.",
    )
    parser.add_argument(
        "--n-top-classes",
        type=int,
        default=6,
        help="Number of most-frequent classes to keep from the importance matrix.",
    )
    parser.add_argument(
        "--result-dir",
        required=True,
        help="Directory to write outputs.",
    )
    parser.add_argument(
        "--n-values",
        default="10,20,50,100,200,500",
        help="Comma-separated list of N values.",
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
        default="top_n_selected_genes",
        help="Prefix for saved outputs.",
    )
    return parser.parse_args()


def build_gene_set(importance_df: pd.DataFrame, n: int) -> list[str]:
    gene_set: set[str] = set()
    for cls in importance_df.index:
        col = importance_df.loc[cls].sort_values(ascending=False)
        top_pos = col.head(n).index
        top_neg = col.tail(n).index
        gene_set.update(top_pos)
        gene_set.update(top_neg)
    return sorted(gene_set)


def main() -> None:
    args = parse_args()
    os.makedirs(args.result_dir, exist_ok=True)

    n_values = [int(x) for x in args.n_values.split(",") if x.strip()]

    adata = sc.read_h5ad(args.adata)

    # Labels from obs column
    if args.label_col not in adata.obs.columns:
        raise KeyError(f"Label column not found in adata.obs: {args.label_col!r}")

    y = adata.obs[args.label_col].astype(str).fillna("no_infection")
    y = pd.Series(y.to_numpy(), index=adata.obs_names, name="label")

    class_names = sorted(pd.unique(y))
    if "no_infection" in class_names:
        class_names = ["no_infection"] + [c for c in class_names if c != "no_infection"]

    # Load full gene importance (rows=classes, cols=genes) and keep top-K classes
    importance_full = pd.read_csv(args.gene_importance_csv, index_col=0)
    label_counts = y.value_counts()
    top_classes = label_counts.head(args.n_top_classes).index.tolist()
    available = [c for c in top_classes if c in importance_full.index]
    if not available:
        raise ValueError("None of the top classes found in gene importance rows.")
    importance_df = importance_full.loc[available]
    print(f"Using {len(available)} classes from importance matrix: {available}")

    for n in n_values:
        gene_set = build_gene_set(importance_df, n)
        genes_present = [g for g in gene_set if g in adata.var_names]
        missing = [g for g in gene_set if g not in adata.var_names]
        if missing:
            print(f"N={n} missing genes in adata.var_names: {len(missing)}")

        if not genes_present:
            print(f"N={n} no genes found in adata.var_names, skipping.")
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

        # Save selected genes
        genes_path = os.path.join(
            args.result_dir, f"{args.metrics_prefix}_n{n}_genes.csv"
        )
        pd.Series(genes_present, name="gene").to_csv(genes_path, index=False)

        # Linear SVM
        linear_pipeline = make_pipeline(
            scaler,
            LinearSVC(
                class_weight="balanced",
                random_state=args.random_state,
            ),
        )
        linear_pipeline.fit(X_train, y_train)
        y_pred = linear_pipeline.predict(X_test)
        y_score = linear_pipeline.decision_function(X_test)

        cm = confusion_matrix(y_test, y_pred, labels=class_names)
        cm_df = pd.DataFrame(cm, index=class_names, columns=class_names)
        cm_path = os.path.join(
            args.result_dir, f"{args.metrics_prefix}_linear_n{n}_confusion_matrix.csv"
        )
        cm_df.to_csv(cm_path)

        svc = linear_pipeline.named_steps["linearsvc"]
        classes = list(svc.classes_)
        if y_score.ndim == 1:
            y_score = y_score[:, None]
        scores_df = pd.DataFrame(y_score, columns=classes, index=idx_test)
        scores_df.insert(0, "y_true", y_test.to_numpy())
        scores_df.insert(1, "y_pred", y_pred)
        scores_path = os.path.join(
            args.result_dir, f"{args.metrics_prefix}_linear_n{n}_scores.csv"
        )
        scores_df.to_csv(scores_path)

        # RBF SVM
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
        y_pred = rbf_pipeline.predict(X_test)
        y_proba = rbf_pipeline.predict_proba(X_test)

        cm = confusion_matrix(y_test, y_pred, labels=class_names)
        cm_df = pd.DataFrame(cm, index=class_names, columns=class_names)
        cm_path = os.path.join(
            args.result_dir, f"{args.metrics_prefix}_rbf_n{n}_confusion_matrix.csv"
        )
        cm_df.to_csv(cm_path)

        svc = rbf_pipeline.named_steps["svc"]
        classes = list(svc.classes_)
        proba_df = pd.DataFrame(y_proba, columns=classes, index=idx_test)
        proba_df.insert(0, "y_true", y_test.to_numpy())
        proba_df.insert(1, "y_pred", y_pred)
        proba_path = os.path.join(
            args.result_dir, f"{args.metrics_prefix}_rbf_n{n}_scores.csv"
        )
        proba_df.to_csv(proba_path)

        print(f"N={n} ({len(genes_present)} genes) done")
        print(f"  linear cm:     {cm_path}")
        print(f"  linear scores: {scores_path}")
        print(f"  rbf cm:        {cm_path}")
        print(f"  rbf scores:    {proba_path}")


if __name__ == "__main__":
    main()
