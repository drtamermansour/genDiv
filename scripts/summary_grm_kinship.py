"""Per-subgroup pairwise GRM kinship summary.

Reads the PLINK2 --make-rel 'square' output (.rel + .rel.id) from the
wholePop run and the USTA_Diversity_Study.gait_bookSize factor file, and
writes a 9-row mean +/- SD table of within-subgroup pairwise GRM values
(off-diagonal only): wholePop + Pacer + Pacer_{LOW,MEDIUM,HIGH} +
Trotter + Trotter_{LOW,MEDIUM,HIGH}.

The pairwise GRM values are the same VanRaden additive-genetic
similarities that PLINK2's --pca uses internally, so this summary
quantifies the kinship structure visible in the PCA. By using the
wholePop GRM (single common AF basis across all 560 animals), the
cross-gait subgroup comparison is on a single yardstick.

Also writes a violin plot of the within-subgroup distributions.

Inputs:
  --rel:    PLINK2 .rel (square N x N) file
  --rel-id: corresponding .rel.id file (FID, IID); one row per matrix
            row in matching order
  --factor: USTA_Diversity_Study.gait_bookSize layout, no header:
            row_index<TAB>IID<TAB>gait_bookSize
            Animals missing from this file are labelled "undefined"
            and contribute to the wholePop row but not to the gait
            or gait x book-size rows.
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

SUBGROUP_ORDER = [
    "wholePop",
    "Pacer", "Pacer_LOW", "Pacer_MEDIUM", "Pacer_HIGH",
    "Trotter", "Trotter_LOW", "Trotter_MEDIUM", "Trotter_HIGH",
]


def _tint(color, factor):
    """Darken an RGB color by a multiplicative factor in [0, 1]."""
    r, g, b = mcolors.to_rgb(color)
    return (r * factor, g * factor, b * factor)


def fmt_mean_sd(values, decimals=4):
    if len(values) == 0:
        return "NA"
    return (f"{np.mean(values):.{decimals}f} +/- "
            f"{np.std(values, ddof=1):.{decimals}f}")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--rel", required=True,
                    help="PLINK2 .rel (square) file.")
    ap.add_argument("--rel-id", required=True,
                    help="PLINK2 .rel.id file (FID, IID).")
    ap.add_argument("--factor", required=True,
                    help="USTA_Diversity_Study.gait_bookSize factor file "
                         "(no header; row_index, IID, gait_bookSize).")
    ap.add_argument("--out-csv", required=True,
                    help="Output CSV path (subgroup-level summary).")
    ap.add_argument("--out-plot", required=True,
                    help="Output PNG (violin plot of within-subgroup distributions).")
    args = ap.parse_args()

    # Sample IDs in matrix row/column order. PLINK2 writes a "#FID IID"
    # comment header line; treat any leading "#" line as a comment.
    ids_df = pd.read_csv(args.rel_id, sep=r"\s+", header=None, comment="#",
                         names=["FID", "IID"], dtype=str)
    ids = ids_df["IID"].tolist()
    N = len(ids)

    grm = np.loadtxt(args.rel)
    if grm.shape != (N, N):
        raise SystemExit(
            f"GRM shape {grm.shape} does not match .rel.id row count {N}"
        )

    factor_df = pd.read_csv(args.factor, sep=r"\s+", header=None,
                            names=["row_idx", "IID", "gait_bookSize"],
                            dtype=str)
    iid_to_factor = dict(zip(factor_df["IID"], factor_df["gait_bookSize"]))

    subs = np.array([iid_to_factor.get(iid, "undefined") for iid in ids])

    # Off-diagonal upper-triangle pairs (i < j).
    i_idx, j_idx = np.triu_indices(N, k=1)
    pair_values = grm[i_idx, j_idx]
    pair_si = subs[i_idx]
    pair_sj = subs[j_idx]

    # --------------------------- summary CSV ---------------------------
    # Per-subgroup animal counts (N_animals), for readers who want to
    # compare against the N column in the other per-subgroup tables.
    n_animals = {
        "wholePop": N,
        "Pacer":   int(sum(s.startswith("Pacer_")   for s in subs)),
        "Trotter": int(sum(s.startswith("Trotter_") for s in subs)),
    }
    for gait in ("Pacer", "Trotter"):
        for bs in BOOK_SIZE_ORDER:
            label = f"{gait}_{bs}"
            n_animals[label] = int(sum(s == label for s in subs))

    rows = [{
        "Subgroup": "wholePop",
        "N_animals": n_animals["wholePop"],
        "N_pairs": int(len(pair_values)),
        "Mean_GRM_kinship": fmt_mean_sd(pair_values),
    }]
    for gait in ("Pacer", "Trotter"):
        gait_i = np.array([s.startswith(gait + "_") for s in pair_si])
        gait_j = np.array([s.startswith(gait + "_") for s in pair_sj])
        gait_mask = gait_i & gait_j
        rows.append({
            "Subgroup": gait,
            "N_animals": n_animals[gait],
            "N_pairs": int(gait_mask.sum()),
            "Mean_GRM_kinship": fmt_mean_sd(pair_values[gait_mask]),
        })
        for bs in BOOK_SIZE_ORDER:
            label = f"{gait}_{bs}"
            sub_mask = (pair_si == label) & (pair_sj == label)
            rows.append({
                "Subgroup": label,
                "N_animals": n_animals[label],
                "N_pairs": int(sub_mask.sum()),
                "Mean_GRM_kinship": fmt_mean_sd(pair_values[sub_mask]),
            })

    df = (pd.DataFrame(rows)
            .set_index("Subgroup")
            .reindex(SUBGROUP_ORDER)
            .reset_index())
    out_csv = Path(args.out_csv)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_csv, index=False)
    print(f"[grm kinship] wrote {out_csv} ({len(df)} rows)")

    # --------------------------- violin plot ---------------------------
    plot_rows = []

    def push(label, mask):
        for v in pair_values[mask]:
            plot_rows.append({"Subgroup": label, "GRM_kinship": float(v)})

    push("wholePop", np.ones_like(pair_values, dtype=bool))
    for gait in ("Pacer", "Trotter"):
        gait_i = np.array([s.startswith(gait + "_") for s in pair_si])
        gait_j = np.array([s.startswith(gait + "_") for s in pair_sj])
        push(gait, gait_i & gait_j)
        for bs in BOOK_SIZE_ORDER:
            push(f"{gait}_{bs}",
                 (pair_si == f"{gait}_{bs}") & (pair_sj == f"{gait}_{bs}"))

    plot_df = pd.DataFrame(plot_rows)
    plot_df["Subgroup"] = pd.Categorical(plot_df["Subgroup"],
                                         categories=SUBGROUP_ORDER, ordered=True)

    # Palette: wholePop = neutral gray; aggregate-gait rows take a soft
    # blue/red base tone; the six gait x book-size strata reuse the
    # BOOK_SIZE_COLORS palette (Pacer keeps base hues; Trotter uses
    # darker tints of the same hues, mirroring Fig4b).
    palette = {
        "wholePop":       "#808080",
        "Pacer":          "#7297C4",
        "Pacer_LOW":      BOOK_SIZE_COLORS["LOW"],
        "Pacer_MEDIUM":   BOOK_SIZE_COLORS["MEDIUM"],
        "Pacer_HIGH":     BOOK_SIZE_COLORS["HIGH"],
        "Trotter":        "#A05D5D",
        "Trotter_LOW":    _tint(BOOK_SIZE_COLORS["LOW"], 0.55),
        "Trotter_MEDIUM": _tint(BOOK_SIZE_COLORS["MEDIUM"], 0.55),
        "Trotter_HIGH":   _tint(BOOK_SIZE_COLORS["HIGH"], 0.55),
    }

    plt.style.use("seaborn-v0_8-darkgrid")
    fig, ax = plt.subplots(figsize=(14, 6))

    # Violin bodies (no inner — boxplot overlay handles the median).
    sns.violinplot(
        data=plot_df, x="Subgroup", y="GRM_kinship",
        order=SUBGROUP_ORDER, ax=ax,
        inner=None, cut=0,
        hue="Subgroup", palette=palette, legend=False,
    )
    # Narrow boxplot overlay with a thick dark median bar for visibility
    # against the white box fill (Fig4b uses white median on colored
    # boxes; here the overlay frame is white so a dark median contrasts).
    sns.boxplot(
        data=plot_df, x="Subgroup", y="GRM_kinship",
        order=SUBGROUP_ORDER, ax=ax,
        width=0.16, showfliers=False, color="white",
        boxprops={"facecolor": "white", "edgecolor": "black", "linewidth": 0.8},
        whiskerprops={"color": "black", "linewidth": 0.8},
        capprops={"color": "black", "linewidth": 0.8},
        medianprops={"color": "black", "linewidth": 2.5},
    )
    # Vertical separators between the three blocks:
    # wholePop | Pacer (aggregate + 3 book-size) | Trotter (aggregate + 3 book-size).
    for x_sep in (0.5, 4.5):
        ax.axvline(x=x_sep, color="dimgray", linewidth=1.0,
                   alpha=0.55, linestyle="--", zorder=1)

    ax.set_title("Within-subgroup pairwise GRM kinship (off-diagonal)",
                 fontsize=14, fontweight="bold")
    ax.set_xlabel("Subgroup", fontsize=12)
    ax.set_ylabel("GRM kinship", fontsize=12)
    plt.setp(ax.get_xticklabels(), rotation=30, ha="right")
    fig.tight_layout()
    out_plot = Path(args.out_plot)
    out_plot.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_plot, dpi=300, bbox_inches="tight")
    print(f"[grm kinship] wrote {out_plot}")


if __name__ == "__main__":
    main()
