"""Plotting for the ROH_common metric (Figures 2-4 of ROH_common.md).

Three modes share one CLI:

    --mode manhattan  : raw per-window f_w line plot along the genome,
                        alternating chromosome shading, no smoothing.
                        Subpanels for the 5 landscapes (genome-wide + 4
                        length classes) for ONE group.
    --mode scatter    : F_ROH vs ROH_common, one panel per gait,
                        colour = book-size category.
    --mode lengthbox  : grouped boxplots of class-specific ROH_common by
                        book-size, one panel per gait.

Style follows scripts/roh_plot.py (seaborn-v0_8-darkgrid, BOOK_SIZE_COLORS).
"""

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from utils import BOOK_SIZE_ORDER, BOOK_SIZE_COLORS

LENGTH_BIN_LABEL = ("1to3", "3to5", "5to10", "more10")
LENGTH_BIN_PRETTY = {
    "1to3":    "1-3 Mb",
    "3to5":    "3-5 Mb",
    "5to10":   "5-10 Mb",
    "more10":  ">10 Mb",
}


# ----------------------------- shared helpers -----------------------------

def load_landscape(path):
    return pd.read_csv(
        path, sep="\t",
        dtype={"chrom": str, "start": np.int64, "end": np.int64,
               "n_w": np.int64, "f_w": np.float64},
    )


def chrom_order_key(c):
    """Sort chromosomes numerically when possible, else lexically."""
    s = c.replace("chr", "")
    try:
        return (0, int(s))
    except ValueError:
        return (1, s)


def add_cumulative_position(land_df):
    """Add a cumulative x coordinate (bp) for whole-genome line plotting and
    return per-chromosome offsets for axis ticks."""
    chroms = sorted(land_df["chrom"].unique(), key=chrom_order_key)
    offset = {}
    cum = 0
    for c in chroms:
        offset[c] = cum
        cum += land_df.loc[land_df["chrom"] == c, "end"].max()
    land_df = land_df.copy()
    land_df["x"] = land_df["start"] + land_df["chrom"].map(offset)
    centers = {c: offset[c] + land_df.loc[land_df["chrom"] == c, "end"].max() / 2
               for c in chroms}
    return land_df, chroms, offset, centers


# ----------------------------- mode: manhattan -----------------------------

def plot_manhattan(args):
    """One PNG per group, five vertically-stacked panels (genome-wide + 4
    length classes). Raw per-window f_w connected by line; chromosomes
    drawn in alternating colour bands; no smoothing."""
    land_dir = Path(args.landscape_dir)
    panels = [("Genome-wide (>=1 Mb)", f"landscape.{args.rg}.tsv")] + \
             [(f"{LENGTH_BIN_PRETTY[l]} class",
               f"landscape.{args.rg}.{l}.tsv") for l in LENGTH_BIN_LABEL]

    fig, axes = plt.subplots(
        nrows=len(panels), ncols=1, figsize=(14, 2.4 * len(panels)),
        sharex=True,
    )
    if len(panels) == 1:
        axes = [axes]

    # Use the genome-wide landscape to define chromosome offsets so all
    # panels share x coords exactly.
    base = load_landscape(land_dir / f"landscape.{args.rg}.tsv")
    base, chroms, offset, centers = add_cumulative_position(base)

    band_colors = ("#3A66A0", "#7FA4D6")  # darker / lighter blue
    for ax, (title, fname) in zip(axes, panels):
        land = load_landscape(land_dir / fname)
        land = land.merge(base[["chrom", "start", "end", "x"]],
                          on=["chrom", "start", "end"], how="left")
        for i, c in enumerate(chroms):
            sub = land[land["chrom"] == c].sort_values("start")
            if sub.empty:
                continue
            ax.plot(sub["x"].values, sub["f_w"].values,
                    color=band_colors[i % 2], linewidth=0.6, antialiased=True)
        ax.set_ylabel("f_w", fontsize=11)
        ax.set_title(f"{args.rg} | {title}", fontsize=12, loc="left", pad=4)
        ax.set_ylim(0, max(0.05, land["f_w"].max() * 1.05))
        ax.grid(False)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
    axes[-1].set_xticks([centers[c] for c in chroms])
    axes[-1].set_xticklabels([c.replace("chr", "") for c in chroms], fontsize=8)
    axes[-1].set_xlabel("Chromosome", fontsize=11)

    fig.suptitle(f"Genome-wide population ROH frequency landscape — {args.rg}",
                 fontsize=14, y=0.995)
    fig.tight_layout(rect=(0, 0, 1, 0.985))
    fig.savefig(args.out, dpi=200, bbox_inches="tight")
    print(f"[plot manhattan] wrote {args.out}")


# ----------------------------- mode: scatter -----------------------------

def _join_with_factor(roh_common_df, froh_df):
    """Merge ROH_common and F_ROH tables on IID, splitting the gait_bookSize
    column into gait + book_size (matches roh_plot.py conventions)."""
    df = roh_common_df.merge(froh_df, on="IID", how="inner")
    df["gait_type"] = df["gait"].str.split("_").str[0]
    df["book_size"] = df["gait"].str.split("_").str[1]
    df["book_size"] = pd.Categorical(df["book_size"],
                                     categories=BOOK_SIZE_ORDER, ordered=True)
    return df


def plot_scatter(args):
    roh_common = pd.read_csv(args.roh_common, sep="\t", na_values=["NA"])
    froh = pd.read_csv(args.froh, sep=r"\s+")
    df = _join_with_factor(roh_common, froh)

    plt.style.use("seaborn-v0_8-darkgrid")
    fig, (ax_p, ax_t) = plt.subplots(2, 1, figsize=(12, 14))

    for ax, gait_name in ((ax_p, "Pacer"), (ax_t, "Trotter")):
        gait_df = df[df["gait_type"] == gait_name]
        for bs in BOOK_SIZE_ORDER:
            sub = gait_df[gait_df["book_size"] == bs]
            if sub.empty:
                continue
            ax.scatter(sub["F_ROH"], sub["ROH_common"],
                       color=BOOK_SIZE_COLORS[bs], s=110, alpha=0.75,
                       edgecolor="black", linewidth=0.6,
                       label=f"Book size: {bs}", zorder=5)
        if len(gait_df) > 1 and gait_df["ROH_common"].notna().sum() > 1:
            x = gait_df["F_ROH"].to_numpy()
            y = gait_df["ROH_common"].to_numpy()
            mask = ~np.isnan(y)
            if mask.sum() > 1:
                coef = np.polyfit(x[mask], y[mask], 1)
                xs = np.linspace(x[mask].min(), x[mask].max(), 100)
                ax.plot(xs, np.polyval(coef, xs), "k--",
                        linewidth=1.5, alpha=0.7, zorder=4,
                        label=f"slope={coef[0]:.3f}")
        ax.set_title(f"{gait_name} — F_ROH vs ROH_common",
                     fontsize=15, fontweight="bold")
        ax.set_xlabel("F_ROH", fontsize=12)
        ax.set_ylabel("ROH_common", fontsize=12)
        ax.legend(fontsize=10, loc="best")
        n_total = len(gait_df)
        n_scored = int(gait_df["ROH_common"].notna().sum())
        ax.text(0.02, 0.98,
                f"N (with ROH_common) = {n_scored} / {n_total}",
                transform=ax.transAxes, fontsize=10, va="top",
                bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.8))

    fig.tight_layout()
    fig.savefig(args.out, dpi=300, bbox_inches="tight")
    print(f"[plot scatter] wrote {args.out}")


# ----------------------------- mode: lengthbox -----------------------------

def _build_lengthbox_long(args):
    roh_common = pd.read_csv(args.roh_common, sep="\t", na_values=["NA"])
    froh = pd.read_csv(args.froh, sep=r"\s+")  # IID + gait (gait_bookSize)
    df = _join_with_factor(roh_common, froh)

    # Long form: one row per (IID, class). Drop the genome-wide ROH_common
    # column before melt so the value_name doesn't collide with an existing
    # column.
    class_cols = [f"ROH_common_{l}" for l in LENGTH_BIN_LABEL]
    keep = ["IID", "gait_type", "book_size"] + class_cols
    long = df[keep].melt(
        id_vars=["IID", "gait_type", "book_size"],
        value_vars=class_cols,
        var_name="class", value_name="ROH_common",
    )
    long["class"] = long["class"].str.replace("ROH_common_", "")
    long["class"] = pd.Categorical(
        long["class"], categories=list(LENGTH_BIN_LABEL), ordered=True
    )
    return long.dropna(subset=["ROH_common"])


def _tint(color, factor):
    """Multiply RGB channels by `factor` (<1 darkens, >1 lightens — clamped)."""
    r, g, b = mcolors.to_rgb(color)
    return tuple(min(1.0, max(0.0, c * factor)) for c in (r, g, b))


def _save_lengthbox_panels(long, out):
    plt.style.use("seaborn-v0_8-darkgrid")
    fig, (ax_p, ax_t) = plt.subplots(2, 1, figsize=(13, 11), sharey=False)

    for ax, gait_name in ((ax_p, "Pacer"), (ax_t, "Trotter")):
        sub = long[long["gait_type"] == gait_name]
        if sub.empty:
            ax.set_title(f"{gait_name} — (no data)", fontsize=14)
            continue
        sns.boxplot(
            data=sub, x="class", y="ROH_common", hue="book_size",
            order=list(LENGTH_BIN_LABEL),
            hue_order=BOOK_SIZE_ORDER,
            palette=BOOK_SIZE_COLORS, ax=ax,
        )
        ax.set_title(f"{gait_name} — class-specific ROH_common by book size",
                     fontsize=14, fontweight="bold")
        ax.set_xlabel("ROH length class (Mb)", fontsize=12)
        ax.set_xticklabels([LENGTH_BIN_PRETTY[l] for l in LENGTH_BIN_LABEL])
        ax.set_ylabel("ROH_common", fontsize=12)
        ax.set_ylim(0, 0.25)
        ax.set_yticks(np.arange(0, 0.2501, 0.025))
        ax.legend(title="Book size", fontsize=10)

    fig.tight_layout()
    fig.savefig(out, dpi=300, bbox_inches="tight")
    print(f"[plot lengthbox] wrote {out}")


def _save_lengthbox_nested(long, out):
    """Single-panel nested boxplot: x = length class, hue = gait × book_size
    (6 levels). Pacer keeps BOOK_SIZE_COLORS; Trotter uses a darker tint of
    the same hue so book size reads off the colour family and gait off the
    shade. Individual points overlaid (stripplot) for n / distribution shape;
    vertical separators mark the Pacer/Trotter boundary within each class."""
    long = long.copy()
    long["gait_book"] = (long["gait_type"].astype(str) + " "
                        + long["book_size"].astype(str))

    # hue_order arranges Pacers (first 3) then Trotters (last 3) within each
    # class — the vertical separator below sits between the two trios.
    hue_order = [f"{g} {bs}" for g in ("Pacer", "Trotter") for bs in BOOK_SIZE_ORDER]
    palette = {}
    for bs in BOOK_SIZE_ORDER:
        palette[f"Pacer {bs}"] = BOOK_SIZE_COLORS[bs]
        palette[f"Trotter {bs}"] = _tint(BOOK_SIZE_COLORS[bs], 0.55)

    plt.style.use("seaborn-v0_8-darkgrid")
    fig, ax = plt.subplots(figsize=(15, 7))
    sns.boxplot(
        data=long, x="class", y="ROH_common", hue="gait_book",
        order=list(LENGTH_BIN_LABEL),
        hue_order=hue_order, palette=palette, ax=ax,
        medianprops={"color": "white", "linewidth": 2.5},
        showfliers=False,
    )
    sns.stripplot(
        data=long, x="class", y="ROH_common", hue="gait_book",
        order=list(LENGTH_BIN_LABEL),
        hue_order=hue_order, dodge=True,
        palette="dark:black", alpha=0.25, size=2.2, jitter=0.18,
        ax=ax, legend=False,
    )
    # Vertical separator between the Pacer trio and Trotter trio at each
    # class centre. Subtle dashed line so it reads as "grouping aid" not
    # "data".
    for x in range(len(LENGTH_BIN_LABEL)):
        ax.axvline(x=x, color="dimgray", linewidth=1.0, alpha=0.55,
                   linestyle="--", zorder=1)

    ax.set_title("Class-specific ROH_common by gait × book size",
                 fontsize=14, fontweight="bold")
    ax.set_xlabel("ROH length class (Mb)", fontsize=12)
    ax.set_xticks(range(len(LENGTH_BIN_LABEL)))
    ax.set_xticklabels([LENGTH_BIN_PRETTY[l] for l in LENGTH_BIN_LABEL])
    ax.set_ylabel("ROH_common", fontsize=12)
    ax.set_ylim(0, 0.25)
    ax.set_yticks(np.arange(0, 0.2501, 0.025))
    ax.legend(title="Gait × Book size", fontsize=9, ncol=1,
              loc="upper left", bbox_to_anchor=(1.02, 1.0),
              borderaxespad=0)

    fig.tight_layout()
    fig.savefig(out, dpi=300, bbox_inches="tight")
    print(f"[plot lengthbox nested] wrote {out}")


def plot_lengthbox(args):
    long = _build_lengthbox_long(args)
    _save_lengthbox_panels(long, args.out)
    if args.nested_out:
        _save_lengthbox_nested(long, args.nested_out)


# ----------------------------- CLI -----------------------------

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="mode", required=True)

    p_m = sub.add_parser("manhattan", help="Figure 2: per-window f_w line plot.")
    p_m.add_argument("--landscape-dir", required=True)
    p_m.add_argument("--rg", required=True)
    p_m.add_argument("--out", required=True)

    p_s = sub.add_parser("scatter", help="Figure 3: F_ROH vs ROH_common.")
    p_s.add_argument("--roh-common", required=True,
                     help="TSV from roh_common_individual.py "
                          "(typically the twoGait concat).")
    p_s.add_argument("--froh", required=True,
                     help="F_ROH table with header IID + gait + F_ROH "
                          "(e.g., roh.L3_Froh_gait_bookSize.txt).")
    p_s.add_argument("--out", required=True)

    p_l = sub.add_parser("lengthbox", help="Figure 4: per-class boxplots.")
    p_l.add_argument("--roh-common", required=True)
    p_l.add_argument("--froh", required=True)
    p_l.add_argument("--out", required=True)
    p_l.add_argument("--nested-out", default=None,
                     help="Optional second output: single-panel nested "
                          "boxplot with gait × book_size as hue.")

    args = parser.parse_args()
    if args.mode == "manhattan":
        plot_manhattan(args)
    elif args.mode == "scatter":
        plot_scatter(args)
    elif args.mode == "lengthbox":
        plot_lengthbox(args)
    else:
        parser.error(f"Unknown mode: {args.mode}")


if __name__ == "__main__":
    main()
