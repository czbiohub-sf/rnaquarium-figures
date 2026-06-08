"""
Replot confusion matrices from saved test_metrics_dict.pkl files.

Usage:
    python replot_conf_mat.py --scan_dir /path/to/scan_devtissue/zebrahub_devtissue_classifier_freeze6layers_testsize20perc

This will find all *_test_metrics_dict.pkl files under the given directory
and regenerate the confusion matrix PDFs using the current plotting code.
"""

import argparse
import pickle
from pathlib import Path

import matplotlib.pyplot as plt
import seaborn as sns
from sklearn import preprocessing
from sklearn.metrics import ConfusionMatrixDisplay

def plot_conf_mat_no_text(conf_mat_df, title, output_dir, output_prefix,
                          custom_class_order=None, height=18, width=18,
                          tick_fontsize=None, title_fontsize=None,
                          xlabel="Predicted label", ylabel="True label",
                          xlabel_fontsize=None, ylabel_fontsize=None):
    sns.set(font_scale=1)
    sns.set_style("whitegrid", {"axes.grid": False})

    if custom_class_order is not None:
        conf_mat_df = conf_mat_df.reindex(index=custom_class_order, columns=custom_class_order)

    display_labels = [
        f"{label} n={conf_mat_df.iloc[i,:].sum():.0f}"
        for i, label in enumerate(conf_mat_df.index)
    ]
    conf_mat = preprocessing.normalize(conf_mat_df.to_numpy(), norm="l1")

    n_classes = len(display_labels)
    if tick_fontsize is None:
        tick_fontsize = max(4, min(10, 300 / n_classes))

    fig, ax = plt.subplots(figsize=(width, height))

    display = ConfusionMatrixDisplay(
        confusion_matrix=conf_mat, display_labels=display_labels
    )
    disp_obj = display.plot(cmap="Blues", values_format=".2f", ax=ax, colorbar=False)

    # Remove cell text
    for text_row in disp_obj.text_:
        for text in text_row:
            text.set_visible(False)

    fig.colorbar(disp_obj.im_, ax=ax, fraction=0.046, pad=0.04, shrink=0.6)
    disp_obj.im_.set_clim(0, 1)

    ax.tick_params(axis='both', which='major', labelsize=tick_fontsize)
    plt.xticks(rotation=45, ha="right")
    ax.set_xlabel(xlabel, fontsize=xlabel_fontsize)
    ax.set_ylabel(ylabel, fontsize=ylabel_fontsize)
    plt.title(title, fontsize=title_fontsize)

    output_file = (Path(output_dir) / f"{output_prefix}_conf_mat").with_suffix(".pdf")
    fig.savefig(output_file, bbox_inches="tight")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--scan_dir", type=str, required=True,
                        help="Directory containing replicate subdirs (e.g. 25Aug14t10pm1/)")
    parser.add_argument("--height", type=int, default=18)
    parser.add_argument("--width", type=int, default=18)
    parser.add_argument("--tick_fontsize", type=float, default=None,
                        help="Font size for x/y tick labels (default: auto-scaled by class count)")
    parser.add_argument("--title", type=str,
                        default="Confusion Matrix for RQ-GF on compound labels, 30/70 train/test split")
    parser.add_argument("--title_fontsize", type=float, default=None,
                        help="Font size for plot title (default: matplotlib default)")
    parser.add_argument("--xlabel", type=str, default="Predicted label")
    parser.add_argument("--ylabel", type=str, default="True label")
    parser.add_argument("--xlabel_fontsize", type=float, default=None,
                        help="Font size for x-axis label (default: matplotlib default)")
    parser.add_argument("--ylabel_fontsize", type=float, default=None,
                        help="Font size for y-axis label (default: matplotlib default)")
    args = parser.parse_args()

    scan_dir = Path(args.scan_dir)
    pkl_files = sorted(scan_dir.rglob("*_test_metrics_dict.pkl"))

    if not pkl_files:
        print(f"No *_test_metrics_dict.pkl files found under {scan_dir}")
        return

    print(f"Found {len(pkl_files)} metrics files")

    for pkl_path in pkl_files:
        output_dir = pkl_path.parent
        # derive output_prefix: strip "_test_metrics_dict.pkl"
        output_prefix = pkl_path.name.replace("_test_metrics_dict.pkl", "")

        print(f"Replotting: {pkl_path}")
        with open(pkl_path, "rb") as f:
            metrics = pickle.load(f)

        conf_mat_df = metrics["conf_matrix"]
        plot_conf_mat_no_text(
            conf_mat_df,
            title=args.title,
            output_dir=str(output_dir),
            output_prefix=output_prefix,
            height=args.height,
            width=args.width,
            tick_fontsize=args.tick_fontsize,
            title_fontsize=args.title_fontsize,
            xlabel=args.xlabel,
            ylabel=args.ylabel,
            xlabel_fontsize=args.xlabel_fontsize,
            ylabel_fontsize=args.ylabel_fontsize,
        )
        print(f"  -> {output_dir / f'{output_prefix}_conf_mat.pdf'}")

    print("Done")


if __name__ == "__main__":
    main()
