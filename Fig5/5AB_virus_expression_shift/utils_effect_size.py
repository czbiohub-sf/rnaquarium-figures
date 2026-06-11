#!/usr/bin/env python3
"""
rna_effect_utils.py — Jupyter-friendly effect-size utilities for RNA-seq
=======================================================================

Drop this file next to your notebook and import the functions you need:

>>> import numpy as np
>>> from rna_effect_utils import (
...     preprocess_log2_tmm, compare_groups,
...     median_diff, trimmed_mean_diff, glass_delta_star,
...     vargha_delaney_A, cliffs_delta_from_A, wasserstein_log2,
...     shift_function, plot_shift_function
... )
>>> ctrl = np.load('ctrl.npy')  # TMM-normalized nonnegative values
>>> pert = np.load('pert.npy')
>>> ctrl_z = preprocess_log2_tmm(ctrl, offset=1.0, winsorize_upper=0.995)
>>> pert_z = preprocess_log2_tmm(pert, offset=1.0, winsorize_upper=0.995)
>>> results = compare_groups(ctrl_z, pert_z, n_boot=400, seed=42,
...                          zero_threshold=None,
...                          shift_quantiles=np.linspace(0.1,0.9,9))
>>> results["median_diff"]  # {'point': ..., 'ci': (..., ...), 'units': 'log2'}
>>> plot_shift_function(results)  # matplotlib figure (decile shifts + 95% CI)

Core features
-------------
- Works directly with **log2(TMM + offset)** values (supply already-transformed arrays
  or call `preprocess_log2_tmm`).
- Robust to heavy tails & unequal variances.
- Bootstrap percentile CIs for all scalar metrics and for decile shift profiles.
- Optional zero-inflation summary for single-cell data.
- Optional balanced-subsampling sensitivity check for imbalanced n.

Dependencies: numpy, scipy, matplotlib (only if plotting).
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional, Tuple

import numpy as np
from numpy.typing import ArrayLike
from scipy.stats import rankdata, trim_mean, wasserstein_distance

# --------------------------- Preprocessing ---------------------------

def _as_1d_float(x: ArrayLike) -> np.ndarray:
    x = np.asarray(x, dtype=float).ravel()
    x = x[np.isfinite(x)]
    if x.size == 0:
        raise ValueError("Array has no finite values after filtering.")
    return x


def preprocess_log2_tmm(
    x: ArrayLike,
    *,
    offset: float = 1.0,
    apply_log2: bool = True,
    winsorize_lower: Optional[float] = None,
    winsorize_upper: Optional[float] = None,
) -> np.ndarray:
    """Return cleaned 1D array, optionally log2(TMM+offset) and winsorized.

    Parameters
    ----------
    x : array-like
        TMM-normalized nonnegative values (or already log-scale if apply_log2=False).
    offset : float, default 1.0
        Offset for log2(TMM + offset). Ignored if apply_log2=False.
    winsorize_lower, winsorize_upper : Optional[float]
        Quantiles in [0,1] used to clip tails (e.g., upper=0.995).
    """
    z = _as_1d_float(x)
    if apply_log2:
        if offset < 0:
            raise ValueError("offset must be >= 0")
        z = np.log2(z + offset)
    if winsorize_lower is not None or winsorize_upper is not None:
        if winsorize_lower is not None and winsorize_lower > 0:
            lo = np.quantile(z, winsorize_lower)
            z = np.maximum(z, lo)
        if winsorize_upper is not None and winsorize_upper < 1:
            hi = np.quantile(z, winsorize_upper)
            z = np.minimum(z, hi)
    return z

# --------------------------- Robust stats & effects ---------------------------

def mad(x: np.ndarray, c: float = 1.4826) -> float:
    med = np.median(x)
    return c * np.median(np.abs(x - med))


def median_diff(ctrl: np.ndarray, pert: np.ndarray) -> float:
    return float(np.median(pert) - np.median(ctrl))


def trimmed_mean_diff(ctrl: np.ndarray, pert: np.ndarray, proportion_to_cut: float = 0.2) -> float:
    return float(trim_mean(pert, proportion_to_cut) - trim_mean(ctrl, proportion_to_cut))


def glass_delta_star(ctrl: np.ndarray, pert: np.ndarray, eps: float = 1e-12) -> Tuple[float, float]:
    """Return (Δ*, MAD_ctrl) where Δ* = median_diff / MAD_ctrl (robust Glass's Δ)."""
    delta = median_diff(ctrl, pert)
    mad_c = mad(ctrl)
    if mad_c <= eps:
        mad_pooled = mad(np.concatenate([ctrl, pert]))
        denom = mad_pooled if mad_pooled > eps else 1.0
    else:
        denom = mad_c
    return float(delta / denom), float(denom)


def vargha_delaney_A(ctrl: np.ndarray, pert: np.ndarray) -> float:
    """Probability a random perturbed value exceeds a random control (AUC of MWU)."""
    x = _as_1d_float(ctrl)
    y = _as_1d_float(pert)
    n_x, n_y = x.size, y.size
    ranks = rankdata(np.concatenate([x, y]), method="average")
    R_y = np.sum(ranks[n_x:])
    U_y = R_y - n_y * (n_y + 1) / 2.0
    return float(U_y / (n_x * n_y))


def cliffs_delta_from_A(A: float) -> float:
    return float(2.0 * A - 1.0)


def wasserstein_log2(ctrl: np.ndarray, pert: np.ndarray) -> float:
    """1-Wasserstein distance (same units as log2 scale)."""
    return float(wasserstein_distance(ctrl, pert))


def shift_function(ctrl: np.ndarray, pert: np.ndarray, qs: Iterable[float]) -> Dict[float, float]:
    """Return {q: Q_q(pert) - Q_q(ctrl)} over quantiles q in (0,1)."""
    qs = list(qs)
    return {float(q): float(np.quantile(pert, q) - np.quantile(ctrl, q)) for q in qs}


def detection_rate(x: np.ndarray, threshold: float = 0.0) -> float:
    return float(np.mean(x > threshold))

# --------------------------- Bootstrapping ---------------------------

def bootstrap_ci_scalar(
    ctrl: np.ndarray,
    pert: np.ndarray,
    stat_fn,
    *,
    n_boot: int = 400,
    seed: Optional[int] = 42,
    alpha: float = 0.05,
) -> Tuple[float, Tuple[float, float]]:
    """Bootstrap percentile CI for a scalar stat_fn(ctrl, pert)."""
    rng = np.random.default_rng(seed)
    stat_hat = float(stat_fn(ctrl, pert))
    n_c, n_p = ctrl.size, pert.size
    boots = np.empty(n_boot, dtype=float)
    for b in range(n_boot):
        idx_c = rng.integers(0, n_c, size=n_c)
        idx_p = rng.integers(0, n_p, size=n_p)
        boots[b] = float(stat_fn(ctrl[idx_c], pert[idx_p]))
    lo, hi = np.quantile(boots, [alpha / 2, 1 - alpha / 2])
    return stat_hat, (float(lo), float(hi))


def bootstrap_ci_deciles(
    ctrl: np.ndarray,
    pert: np.ndarray,
    qs: Iterable[float],
    *,
    n_boot: int = 400,
    seed: Optional[int] = 42,
    alpha: float = 0.05,
) -> Dict[float, Dict[str, Tuple[float, float]]]:
    qs = list(qs)
    base = shift_function(ctrl, pert, qs)
    rng = np.random.default_rng(seed)
    n_c, n_p = ctrl.size, pert.size
    boot_mat = np.empty((n_boot, len(qs)), dtype=float)
    for b in range(n_boot):
        idx_c = rng.integers(0, n_c, size=n_c)
        idx_p = rng.integers(0, n_p, size=n_p)
        diffs = shift_function(ctrl[idx_c], pert[idx_p], qs)
        boot_mat[b, :] = [diffs[q] for q in qs]
    out: Dict[float, Dict[str, Tuple[float, float]]] = {}
    for j, q in enumerate(qs):
        lo, hi = np.quantile(boot_mat[:, j], [alpha / 2, 1 - alpha / 2])
        out[float(q)] = {"point": float(base[q]), "ci": (float(lo), float(hi))}
    return out


def bootstrap_glass_delta(
    ctrl: np.ndarray,
    pert: np.ndarray,
    *,
    n_boot: int = 400,
    seed: Optional[int] = 42,
    alpha: float = 0.05,
) -> Tuple[Tuple[float, Tuple[float, float]], Tuple[float, Tuple[float, float]]]:
    """Return (Δ* point, CI), (MAD_ctrl point, CI) from bootstraps."""
    dstar_point, mad_ctrl_point = glass_delta_star(ctrl, pert)
    rng = np.random.default_rng(seed)
    n_c, n_p = ctrl.size, pert.size
    vals_delta = np.empty(n_boot)
    vals_mad = np.empty(n_boot)
    for b in range(n_boot):
        idx_c = rng.integers(0, n_c, size=n_c)
        idx_p = rng.integers(0, n_p, size=n_p)
        dstar, mad_c = glass_delta_star(ctrl[idx_c], pert[idx_p])
        vals_delta[b] = dstar
        vals_mad[b] = mad_c
    d_lo, d_hi = np.quantile(vals_delta, [alpha / 2, 1 - alpha / 2])
    m_lo, m_hi = np.quantile(vals_mad, [alpha / 2, 1 - alpha / 2])
    return (float(dstar_point), (float(d_lo), float(d_hi))), (float(mad_ctrl_point), (float(m_lo), float(m_hi)))

# --------------------------- Sensitivity to imbalance ---------------------------

def sensitivity_balanced_subsampling(
    ctrl: np.ndarray,
    pert: np.ndarray,
    *,
    reps: int = 200,
    seed: Optional[int] = 123,
) -> Dict[str, Dict[str, float]]:
    """Repeatedly draw min(n_ctrl,n_pert) from each group and recompute metrics.
    Returns mean and std for median Δ, A, and Wasserstein.
    """
    if reps <= 0:
        return {}
    rng = np.random.default_rng(seed)
    n_c, n_p = ctrl.size, pert.size
    m = min(n_c, n_p)
    med_vals = np.empty(reps)
    A_vals = np.empty(reps)
    W_vals = np.empty(reps)
    for r in range(reps):
        idx_c = rng.choice(n_c, size=m, replace=False)
        idx_p = rng.choice(n_p, size=m, replace=False)
        x = ctrl[idx_c]
        y = pert[idx_p]
        med_vals[r] = median_diff(x, y)
        A_vals[r] = vargha_delaney_A(x, y)
        W_vals[r] = wasserstein_log2(x, y)
    return {
        "median_diff": {"mean": float(np.mean(med_vals)), "std": float(np.std(med_vals, ddof=1))},
        "A": {"mean": float(np.mean(A_vals)), "std": float(np.std(A_vals, ddof=1))},
        "wasserstein": {"mean": float(np.mean(W_vals)), "std": float(np.std(W_vals, ddof=1))},
    }

# --------------------------- Atomic analysis functions ---------------------------

def compute_median_shift(ctrl: ArrayLike, pert: ArrayLike, *, n_boot: int = 400, alpha: float = 0.05, seed: Optional[int] = 42) -> Dict:
    x = _as_1d_float(ctrl); y = _as_1d_float(pert)
    point, ci = bootstrap_ci_scalar(x, y, median_diff, n_boot=n_boot, seed=seed, alpha=alpha)
    return {"point": point, "ci": ci, "units": "log2"}


def compute_trimmed_shift(ctrl: ArrayLike, pert: ArrayLike, *, trimmed_prop: float = 0.2, n_boot: int = 400, alpha: float = 0.05, seed: Optional[int] = 42) -> Dict:
    x = _as_1d_float(ctrl); y = _as_1d_float(pert)
    tm_fn = lambda c, p: trimmed_mean_diff(c, p, proportion_to_cut=trimmed_prop)
    point, ci = bootstrap_ci_scalar(x, y, tm_fn, n_boot=n_boot, seed=seed, alpha=alpha)
    return {"point": point, "ci": ci, "trimmed_prop": float(trimmed_prop), "units": "log2"}


def compute_glass_delta(ctrl: ArrayLike, pert: ArrayLike, *, n_boot: int = 400, alpha: float = 0.05, seed: Optional[int] = 42) -> Dict:
    x = _as_1d_float(ctrl); y = _as_1d_float(pert)
    (dstar_point, dstar_ci), (mad_point, mad_ci) = bootstrap_glass_delta(x, y, n_boot=n_boot, seed=seed, alpha=alpha)
    return {
        "glass_delta_star": {"point": dstar_point, "ci": dstar_ci},
        "mad_control": {"point": mad_point, "ci": mad_ci, "units": "log2"},
    }


def compute_stochastic_dominance(ctrl: ArrayLike, pert: ArrayLike, *, n_boot: int = 400, alpha: float = 0.05, seed: Optional[int] = 42) -> Dict:
    x = _as_1d_float(ctrl); y = _as_1d_float(pert)
    A_point, A_ci = bootstrap_ci_scalar(x, y, vargha_delaney_A, n_boot=n_boot, seed=seed, alpha=alpha)
    return {
        "A_vargha_delaney": {"point": A_point, "ci": A_ci},
        "cliffs_delta": {"point": cliffs_delta_from_A(A_point), "ci": (2*A_ci[0]-1, 2*A_ci[1]-1)},
    }


def compute_wasserstein(ctrl: ArrayLike, pert: ArrayLike, *, n_boot: int = 400, alpha: float = 0.05, seed: Optional[int] = 42) -> Dict:
    x = _as_1d_float(ctrl); y = _as_1d_float(pert)
    point, ci = bootstrap_ci_scalar(x, y, wasserstein_log2, n_boot=n_boot, seed=seed, alpha=alpha)
    return {"point": point, "ci": ci, "units": "log2"}


def compute_shift_profile(ctrl: ArrayLike, pert: ArrayLike, *, qs: Iterable[float] = (0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9), n_boot: int = 400, alpha: float = 0.05, seed: Optional[int] = 42) -> Dict[float, Dict[str, Tuple[float, float]]]:
    x = _as_1d_float(ctrl); y = _as_1d_float(pert)
    return bootstrap_ci_deciles(x, y, qs=qs, n_boot=n_boot, seed=seed, alpha=alpha)


def compute_detection_summaries(ctrl: ArrayLike, pert: ArrayLike, *, threshold: float = 0.0, n_boot: int = 400, alpha: float = 0.05, seed: Optional[int] = 42) -> Dict:
    x = _as_1d_float(ctrl); y = _as_1d_float(pert)
    det_fn = lambda c, p: detection_rate(p, threshold) - detection_rate(c, threshold)
    det_point, det_ci = bootstrap_ci_scalar(x, y, det_fn, n_boot=n_boot, seed=seed, alpha=alpha)
    def med_pos(c, p):
        c_pos = c[c > threshold]; p_pos = p[p > threshold]
        if c_pos.size == 0 or p_pos.size == 0:
            return np.nan
        return median_diff(c_pos, p_pos)
    medpos_point, medpos_ci = bootstrap_ci_scalar(x, y, med_pos, n_boot=n_boot, seed=seed, alpha=alpha)
    return {
        "detection_rate_diff": {"point": det_point, "ci": det_ci, "threshold": float(threshold)},
        "median_diff_among_positives": {"point": medpos_point, "ci": medpos_ci, "units": "log2"},
    }


def compute_balanced_subsampling_sensitivity(ctrl: ArrayLike, pert: ArrayLike, *, reps: int = 200, seed: Optional[int] = 123) -> Dict:
    x = _as_1d_float(ctrl); y = _as_1d_float(pert)
    return sensitivity_balanced_subsampling(x, y, reps=reps, seed=seed)

# --------------------------- Aggregated orchestrator ---------------------------

def compare_groups(
    ctrl: ArrayLike,
    pert: ArrayLike,
    *,
    n_boot: int = 400,
    alpha: float = 0.05,
    seed: Optional[int] = 42,
    compute_trimmed: bool = True,
    trimmed_prop: float = 0.2,
    shift_quantiles: Iterable[float] = (0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9),
    zero_threshold: Optional[float] = None,
    balanced_subsample_reps: int = 0,
) -> Dict:
    """Aggregate all analyses; inputs should already be on the log2 scale."""
    x = _as_1d_float(ctrl); y = _as_1d_float(pert)
    results: Dict = {"n_ctrl": int(x.size), "n_pert": int(y.size), "alpha": float(alpha), "n_boot": int(n_boot)}

    results["median_diff"] = compute_median_shift(x, y, n_boot=n_boot, alpha=alpha, seed=seed)
    if compute_trimmed:
        results["trimmed_mean_diff"] = compute_trimmed_shift(x, y, trimmed_prop=trimmed_prop, n_boot=n_boot, alpha=alpha, seed=seed)
    results.update(compute_glass_delta(x, y, n_boot=n_boot, alpha=alpha, seed=seed))
    results.update(compute_stochastic_dominance(x, y, n_boot=n_boot, alpha=alpha, seed=seed))
    results["wasserstein"] = compute_wasserstein(x, y, n_boot=n_boot, alpha=alpha, seed=seed)
    results["shift_function"] = compute_shift_profile(x, y, qs=shift_quantiles, n_boot=n_boot, alpha=alpha, seed=seed)

    if zero_threshold is not None:
        results.update(compute_detection_summaries(x, y, threshold=zero_threshold, n_boot=n_boot, alpha=alpha, seed=seed))

    sens = compute_balanced_subsampling_sensitivity(x, y, reps=balanced_subsample_reps, seed=seed)
    if sens:
        results["balanced_subsampling_sensitivity"] = sens

    return results

# --------------------------- Plotting helper ---------------------------

def plot_shift_function(results: Dict, *, title: Optional[str] = None):
    """Plot decile (or custom) shift function with 95% bootstrap CIs.

    Returns the matplotlib figure object. Requires matplotlib installed.
    """
    import matplotlib.pyplot as plt

    sf = results.get("shift_function", {})
    if not sf:
        raise ValueError("No 'shift_function' found in results. Run compare_groups with shift_quantiles.")

    qs = sorted(sf.keys())
    pts = [sf[q]["point"] for q in qs]
    los = [sf[q]["ci"][0] for q in qs]
    his = [sf[q]["ci"][1] for q in qs]

    fig, ax = plt.subplots()
    ax.plot(qs, pts, marker='o')
    ax.fill_between(qs, los, his, alpha=0.2)
    ax.axhline(0.0, linewidth=1.0)
    ax.set_xlabel("Quantile q")
    ax.set_ylabel("Pert − Ctrl (log2 units)")
    if title:
        ax.set_title(title)
    return fig




import json
from pathlib import Path

def df_nonzero_rows_by_column_to_json(df, json_path, *, row_name=None, orient="columns", indent=2):
    """
    Build a dict mapping each column -> list of row labels (or values from `row_name`)
    where df[col] is non-zero, then write it to a JSON file.

    Parameters
    ----------
    df : pandas.DataFrame
    json_path : str | Path
        Output JSON filepath.
    row_name : str | None
        If None, use df.index labels. If provided, use values from df[row_name].
        (Row labels returned will be strings in JSON.)
    orient : {"columns"}
        Kept for future extensibility; currently only 'columns' behavior is implemented.
    indent : int
        JSON indentation.
    """
    if orient != "columns":
        raise ValueError("Only orient='columns' is supported.")

    # Decide what to use as the "row name" to return
    if row_name is None:
        row_labels = df.index
    else:
        if row_name not in df.columns:
            raise KeyError(f"row_name '{row_name}' not found in df.columns")
        row_labels = df[row_name]
    
    out = {}
    for col in df.columns:
        if col == row_name:
            continue  # don't treat the label column as a data column
        
        # Non-zero mask (treat NaN as zero/False)
        s = df[col]
        mask = s.fillna(0).ne(0)
        
        # Collect corresponding row labels
        out[col] = [str(x) for x in row_labels[mask].tolist()]

    json_path = Path(json_path)
    with json_path.open("w", encoding="utf-8") as f:
        json.dump(out, f, indent=indent, sort_keys=True)

    return out

# --- Example usage ---
# out_dict = df_nonzero_rows_by_column_to_json(df, "nonzero_rows_by_col.json")
# or, if you want to use a column as the row identifier:
# out_dict = df_nonzero_rows_by_column_to_json(df, "nonzero_rows_by_col.json", row_name="gene_id")