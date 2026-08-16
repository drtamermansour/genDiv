#!/usr/bin/env python3
"""
Replicate the top-right 'Relationship Comparison' subplot from
Robust_Matrix_Comparison_Enhanced.wholePop.png (Kinship_Std vs Kinship_ROH),
coloured once by pair-level Gait and once by pair-level Book Size, combining
the 1Mb / 5Mb / 10Mb ROH cutoffs into a 1x3 row per colour scheme.

Inputs (all relative to --output-dir):
  rep_ROHRM/roh_{1,5,10}Mb.Threshold_{sd}SD/Pairwise_Differences.wholePop.csv
  preprocess/USTA_Diversity_Study.bookSize   (tab-sep: idx, IID, bookSize)

Outputs:
  explore/relationship_comparison_byGait.wholePop.png
  explore/relationship_comparison_byBookSize.wholePop.png
  explore/relationship_comparison_byBookSize_HIGHgroups.wholePop.png
  explore/relationship_comparison_byBookSize_LOWandMEDIUMgroups.wholePop.png
"""

import argparse
import os
import sys
import pandas as pd
import matplotlib.pyplot as plt


def load_pairwise(out_dir, cutoff_mb, sd_label):
    subdir = f"roh_{int(float(cutoff_mb))}Mb.Threshold_{sd_label}"
    path = os.path.join(out_dir, "rep_ROHRM", subdir, "Pairwise_Differences.wholePop.csv")
    if not os.path.exists(path):
        sys.exit(f"[ERROR] Missing input: {path}")
    return pd.read_csv(path)


def pair_key_unordered(a, b):
    return "-".join(sorted([a, b]))


def add_gait_pair(df):
    df = df[(df["Pheno1"] != "undefined") & (df["Pheno2"] != "undefined")].copy()
    df["GaitPair"] = [pair_key_unordered(a, b) for a, b in zip(df["Pheno1"], df["Pheno2"])]
    return df


def add_booksize_pair(df, booksize_map):
    df = df.copy()
    df["BS1"] = df["ID1"].map(booksize_map)
    df["BS2"] = df["ID2"].map(booksize_map)
    df = df.dropna(subset=["BS1", "BS2"])
    df["BookSizePair"] = [pair_key_unordered(a, b) for a, b in zip(df["BS1"], df["BS2"])]
    return df


def scatter_panel(ax, df, color_col, category_order, palette, title):
    for cat in category_order:
        sub = df[df[color_col] == cat]
        if sub.empty:
            continue
        ax.scatter(
            sub["Kinship_Std"], sub["Kinship_ROH"],
            s=3, alpha=0.35, c=palette[cat], label=f"{cat} (n={len(sub)})",
            linewidths=0,
        )
    lo = min(df["Kinship_Std"].min(), df["Kinship_ROH"].min())
    hi = max(df["Kinship_Std"].max(), df["Kinship_ROH"].max())
    ax.plot([lo, hi], [lo, hi], "r--", linewidth=1.2, label="1:1")
    ax.set_xlabel(r"$G_{SNP}$")
    ax.set_ylabel(r"$G_{ROH}$")
    ax.set_title(title)
    ax.legend(markerscale=3, fontsize=8, loc="upper left")


def build_figure(dfs_by_cutoff, color_col, category_order, palette, out_path,
                 suptitle):
    fig, axes = plt.subplots(1, 3, figsize=(18, 6), sharex=True, sharey=True)
    for ax, (cutoff, df) in zip(axes, dfs_by_cutoff.items()):
        scatter_panel(ax, df, color_col, category_order, palette,
                      title=f"{cutoff} Mb cutoff")
    fig.suptitle(suptitle, fontsize=13)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(out_path, dpi=200)
    plt.close(fig)
    print(f"  -> {out_path}")


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--output-dir", required=True,
                   help="Pipeline run directory (results_<timestamp>)")
    p.add_argument("--cutoffs", nargs="+", default=["1.0", "5.0", "10.0"])
    p.add_argument("--sd-label", default="3SD")
    args = p.parse_args()

    out_dir = args.output_dir
    explore_dir = os.path.join(out_dir, "explore")
    os.makedirs(explore_dir, exist_ok=True)

    booksize_file = os.path.join(out_dir, "preprocess", "USTA_Diversity_Study.bookSize")
    if not os.path.exists(booksize_file):
        sys.exit(f"[ERROR] Missing bookSize metadata: {booksize_file}")
    bs = pd.read_csv(booksize_file, sep="\t", header=None,
                     names=["idx", "IID", "BookSize"])
    booksize_map = dict(zip(bs["IID"].astype(str), bs["BookSize"]))

    raw = {c: load_pairwise(out_dir, c, args.sd_label) for c in args.cutoffs}

    ## Gait pairs: 3 unordered categories (Pacer-Pacer, Pacer-Trotter, Trotter-Trotter)
    gait_dfs = {c: add_gait_pair(df) for c, df in raw.items()}
    gait_cats = ["Pacer-Pacer", "Pacer-Trotter", "Trotter-Trotter"]
    gait_palette = {
        "Pacer-Pacer":     "#1f77b4",
        "Pacer-Trotter":   "#7f7f7f",
        "Trotter-Trotter": "#d62728",
    }
    build_figure(gait_dfs, "GaitPair", gait_cats, gait_palette,
                 os.path.join(explore_dir, "relationship_comparison_byGait.wholePop.png"),
                 suptitle="Relationship Comparison (wholePop) — coloured by gait pair")

    ## Book-size pairs: 6 unordered categories
    bs_dfs = {c: add_booksize_pair(df, booksize_map) for c, df in raw.items()}
    bs_cats_all = ["HIGH-HIGH", "HIGH-LOW", "HIGH-MEDIUM",
                   "LOW-LOW",   "LOW-MEDIUM", "MEDIUM-MEDIUM"]
    bs_palette = {
        "HIGH-HIGH":     "#d62728",
        "HIGH-LOW":      "#9467bd",
        "HIGH-MEDIUM":   "#ff7f0e",
        "LOW-LOW":       "#1f77b4",
        "LOW-MEDIUM":    "#17becf",
        "MEDIUM-MEDIUM": "#2ca02c",
    }
    build_figure(bs_dfs, "BookSizePair", bs_cats_all, bs_palette,
                 os.path.join(explore_dir, "relationship_comparison_byBookSize.wholePop.png"),
                 suptitle="Relationship Comparison (wholePop) — coloured by book-size pair (all)")

    ## Focused subsets: keep category palette identical so the same colour
    ## means the same pair-category across all four book-size figures.
    ## Axes re-fit to each subset for better resolution of detail.
    high_cats = ["HIGH-HIGH", "HIGH-LOW", "HIGH-MEDIUM"]
    low_med_cats = ["LOW-LOW", "LOW-MEDIUM", "MEDIUM-MEDIUM"]

    high_dfs = {c: df[df["BookSizePair"].isin(high_cats)].copy()
                for c, df in bs_dfs.items()}
    low_med_dfs = {c: df[df["BookSizePair"].isin(low_med_cats)].copy()
                   for c, df in bs_dfs.items()}

    build_figure(high_dfs, "BookSizePair", high_cats, bs_palette,
                 os.path.join(explore_dir, "relationship_comparison_byBookSize_HIGHgroups.wholePop.png"),
                 suptitle="Relationship Comparison (wholePop) — pairs involving HIGH (HIGH-HIGH / HIGH-MEDIUM / HIGH-LOW)")
    build_figure(low_med_dfs, "BookSizePair", low_med_cats, bs_palette,
                 os.path.join(explore_dir, "relationship_comparison_byBookSize_LOWandMEDIUMgroups.wholePop.png"),
                 suptitle="Relationship Comparison (wholePop) — LOW and MEDIUM pairs only (LOW-LOW / LOW-MEDIUM / MEDIUM-MEDIUM)")


if __name__ == "__main__":
    main()
