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
from sklearn.svm import LinearSVC
from sklearn.metrics import confusion_matrix


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Train multiclass linear SVM on gene expression features with hard labels."
    )
    parser.add_argument(
        "--adata",
        required=True,
        help="Path to input .h5ad file.",
    )
    parser.add_argument(
        "--virus-counts-csv",
        required=False,
        default=None,
        help="(Unused) Path to virus_counts_df.csv. Kept for backward compatibility.",
    )
    parser.add_argument(
        "--label-col",
        default="virus_classification",
        help="Column in adata.obs to use as ground-truth labels.",
    )
    parser.add_argument(
        "--infected-label",
        default="infected",
        help="Label name to use for infected cells (non-no_infection).",
    )
    parser.add_argument(
        "--result-dir",
        required=True,
        help="Directory to write outputs.",
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
        default="binary_linear_svm",
        help="Prefix for saved outputs.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    os.makedirs(args.result_dir, exist_ok=True)

    adata = sc.read_h5ad(args.adata)

    if args.label_col not in adata.obs.columns:
        raise KeyError(f"Label column not found in adata.obs: {args.label_col!r}")

    y_raw = adata.obs[args.label_col].astype(str).fillna("no_infection")
    y_bin = np.where(y_raw.to_numpy() == "no_infection", "no_infection", args.infected_label)
    y = pd.Series(y_bin, index=adata.obs_names, name="label")
    class_names = ["no_infection", args.infected_label]

    X = adata.X
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

    # Linear SVM (no calibration)
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

    # Save overall confusion matrix
    cm = confusion_matrix(y_test, y_pred, labels=class_names)
    cm_df = pd.DataFrame(cm, index=class_names, columns=class_names)
    cm_path = os.path.join(args.result_dir, f"{args.metrics_prefix}_confusion_matrix.csv")
    cm_df.to_csv(cm_path)

    # Save score matrix (decision function)
    # Binary LinearSVC returns 1-D decision_function (positive → second class)
    svc = linear_pipeline.named_steps["linearsvc"]
    scores_df = pd.DataFrame({
        "y_true": y_test.to_numpy(),
        "y_pred": y_pred,
        "decision_score": y_score,
    }, index=idx_test)
    scores_path = os.path.join(args.result_dir, f"{args.metrics_prefix}_scores.csv")
    scores_df.to_csv(scores_path)

    # Save gene importance from LinearSVC coefficients
    # Binary case: coef_ is (1, n_features); rows labelled by the two classes
    coef = svc.coef_
    coef_df = pd.DataFrame(
        coef,
        index=[f"{svc.classes_[0]}_vs_{svc.classes_[1]}"],
        columns=adata.var_names,
    )
    coef_path = os.path.join(args.result_dir, f"{args.metrics_prefix}_gene_importance.csv")
    coef_df.to_csv(coef_path)

    print(f"Saved confusion matrix: {cm_path}")
    print(f"Saved scores: {scores_path}")
    print(f"Saved gene importance: {coef_path}")


if __name__ == "__main__":
    main()
