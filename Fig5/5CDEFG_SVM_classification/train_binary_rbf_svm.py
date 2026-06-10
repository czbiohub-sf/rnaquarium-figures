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
from sklearn.svm import SVC
from sklearn.metrics import confusion_matrix


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Train multiclass RBF SVM on gene expression features with hard labels."
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
        default="binary_rbf_svm",
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

    # Save overall confusion matrix
    cm = confusion_matrix(y_test, y_pred, labels=class_names)
    cm_df = pd.DataFrame(cm, index=class_names, columns=class_names)
    cm_path = os.path.join(args.result_dir, f"{args.metrics_prefix}_confusion_matrix.csv")
    cm_df.to_csv(cm_path)

    # Save probability matrix
    svc = rbf_pipeline.named_steps["svc"]
    classes = list(svc.classes_)
    proba_df = pd.DataFrame(y_proba, columns=classes, index=idx_test)
    proba_df.insert(0, "y_true", y_test.to_numpy())
    proba_df.insert(1, "y_pred", y_pred)
    # Convenience column for binary workflows: P(infected)
    if args.infected_label in classes:
        proba_df.insert(2, "p_infected", proba_df[args.infected_label].to_numpy())
    else:
        proba_df.insert(2, "p_infected", np.nan)
    proba_path = os.path.join(args.result_dir, f"{args.metrics_prefix}_scores.csv")
    proba_df.to_csv(proba_path)

    print(f"Saved confusion matrix: {cm_path}")
    print(f"Saved probabilities: {proba_path}")


if __name__ == "__main__":
    main()
