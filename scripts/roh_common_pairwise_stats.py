"""Pairwise Mann-Whitney U + Cohen's d for ROH_common across the six
Trotter/Pacer x {LOW, MEDIUM, HIGH} subgroups.

Five metrics: ROH_common (>=1 Mb genome-wide) + ROH_common_{1to3,3to5,
5to10,more10} length classes. 15 unordered pairs x 5 metrics = 75 rows.

Output TSV columns:
    group_A group_B metric n_A n_B mean_A mean_B median_A median_B
    stat p_value p_adj_bonferroni cohens_d

p_adj_bonferroni is the Bonferroni FWER-adjusted p-value computed within
each metric (15 tests per family): p_adj = min(p * m, 1) with m the number
of tested pairs in the family. stat = Mann-Whitney U for sample A.
"""

import argparse
import itertools
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

from utils import BOOK_SIZE_ORDER

METRICS = (
    "ROH_common",
    "ROH_common_1to3",
    "ROH_common_3to5",
    "ROH_common_5to10",
    "ROH_common_more10",
)
GROUPS = [f"{g}_{bs}" for g in ("Pacer", "Trotter") for bs in BOOK_SIZE_ORDER]


def cohens_d(a, b):
    """Pooled-SD Cohen's d. NaN when either group has <2 finite values or
    pooled variance collapses."""
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    na, nb = len(a), len(b)
    if na < 2 or nb < 2:
        return np.nan
    va = float(np.var(a, ddof=1))
    vb = float(np.var(b, ddof=1))
    pooled = ((na - 1) * va + (nb - 1) * vb) / (na + nb - 2)
    if pooled <= 0:
        return np.nan
    return (float(np.mean(a)) - float(np.mean(b))) / np.sqrt(pooled)


def bonferroni_adjust(pvals):
    """Bonferroni FWER adjustment within a family: p_adj = min(p * m, 1),
    where m is the number of non-NaN tests in the family. Returns adjusted
    p-values in the original input order; NaN inputs stay NaN."""
    p = np.asarray(pvals, dtype=float)
    out = np.full_like(p, np.nan, dtype=float)
    mask = ~np.isnan(p)
    m = int(mask.sum())
    if m == 0:
        return out
    out[mask] = np.clip(p[mask] * m, 0, 1)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--roh-common", required=True,
                    help="roh_common.twoGait.tsv (per-IID ROH_common + 4 "
                         "length-class columns).")
    ap.add_argument("--froh", required=True,
                    help="roh.L3_Froh_gait_bookSize.txt; supplies the "
                         "gait_bookSize factor per IID in column 'gait'.")
    ap.add_argument("--out", required=True, help="Output TSV path.")
    args = ap.parse_args()

    rc = pd.read_csv(args.roh_common, sep="\t", na_values=["NA"])
    froh = pd.read_csv(args.froh, sep=r"\s+")
    if "gait" not in froh.columns:
        raise SystemExit("--froh file must contain a 'gait' column "
                         "(gait_bookSize labels).")
    df = rc.merge(froh[["IID", "gait"]], on="IID", how="inner")
    df = df.rename(columns={"gait": "subgroup"})

    missing_metrics = [m for m in METRICS if m not in df.columns]
    if missing_metrics:
        raise SystemExit(f"Missing metric columns: {missing_metrics}")

    rows = []
    for metric in METRICS:
        for a, b in itertools.combinations(GROUPS, 2):
            x = df.loc[df["subgroup"] == a, metric].dropna().to_numpy(dtype=float)
            y = df.loc[df["subgroup"] == b, metric].dropna().to_numpy(dtype=float)
            na, nb = len(x), len(y)
            row = {
                "group_A": a, "group_B": b, "metric": metric,
                "n_A": na, "n_B": nb,
                "mean_A":   float(np.mean(x))   if na else np.nan,
                "mean_B":   float(np.mean(y))   if nb else np.nan,
                "median_A": float(np.median(x)) if na else np.nan,
                "median_B": float(np.median(y)) if nb else np.nan,
                "stat": np.nan, "p_value": np.nan,
                "cohens_d": cohens_d(x, y),
            }
            if na >= 2 and nb >= 2:
                u, pv = mannwhitneyu(x, y, alternative="two-sided")
                row["stat"] = float(u)
                row["p_value"] = float(pv)
            rows.append(row)

    out_df = pd.DataFrame(rows)
    out_df["p_adj_bonferroni"] = np.nan
    for metric in METRICS:
        mask = out_df["metric"] == metric
        out_df.loc[mask, "p_adj_bonferroni"] = bonferroni_adjust(
            out_df.loc[mask, "p_value"].to_numpy()
        )

    cols = ["group_A", "group_B", "metric",
            "n_A", "n_B",
            "mean_A", "mean_B", "median_A", "median_B",
            "stat", "p_value", "p_adj_bonferroni", "cohens_d"]
    out_df = out_df[cols]

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_df.to_csv(out_path, sep="\t", index=False,
                  float_format="%.6g", na_rep="NA")
    print(f"[pairwise stats] wrote {out_path} ({len(out_df)} rows, "
          f"{out_df['p_value'].notna().sum()} tested)")


if __name__ == "__main__":
    main()
