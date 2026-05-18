"""Per-subgroup ROH_common summary table (mean +/- SD).

Produces a 9-row CSV with one row per analysis subgroup:

    wholePop, Pacer, Pacer_LOW, Pacer_MEDIUM, Pacer_HIGH,
    Trotter, Trotter_LOW, Trotter_MEDIUM, Trotter_HIGH

and five metric columns (genome-wide ROH_common + the four length classes:
1to3, 3to5, 5to10, more10).

NOTE on reference frames. The wholePop row uses scores computed against
the cohort-wide landscape (N = 560 in our cohort); per-gait and
per-booksize rows use scores computed against the within-gait landscape
(N = 271 per gait). The wholePop value is therefore on a different
reference scale than the other eight rows and should not be directly
compared to them numerically -- footnote this in the manuscript table
caption.
"""

import argparse
from pathlib import Path

import pandas as pd

METRICS = (
    "ROH_common",
    "ROH_common_1to3",
    "ROH_common_3to5",
    "ROH_common_5to10",
    "ROH_common_more10",
)


def fmt_mean_sd(series, decimals=4):
    """'mean +/- sd' with NA handling. n in the table is the non-NA count."""
    s = series.dropna()
    if len(s) == 0:
        return "NA"
    return f"{s.mean():.{decimals}f} +/- {s.std():.{decimals}f}"


def row_for(df, label):
    """One row: Subgroup, N (total rows in the slice), then a mean +/- SD
    cell per metric (each cell uses its own non-NA n; per-metric NA counts
    can be audited via n_windows_* in the per-individual TSVs)."""
    rec = {"Subgroup": label, "N": len(df)}
    for m in METRICS:
        rec[m] = fmt_mean_sd(df[m])
    return rec


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--wholepop", required=True,
                    help="roh_common.wholePop.tsv (cohort scored against the "
                         "cohort-wide landscape).")
    ap.add_argument("--threebooksize", required=True,
                    help="roh_common.threeBooksize.tsv (gait-level scoring "
                         "with a gait_bookSize column).")
    ap.add_argument("--out", required=True, help="Output CSV path.")
    args = ap.parse_args()

    wp = pd.read_csv(args.wholepop, sep="\t", na_values=["NA"])
    tb = pd.read_csv(args.threebooksize, sep="\t", na_values=["NA"])
    if "gait_bookSize" not in tb.columns:
        raise SystemExit("--threebooksize file must contain a 'gait_bookSize' column.")
    tb["gait"] = tb["gait_bookSize"].str.split("_").str[0]

    rows = [
        row_for(wp,                                            "wholePop"),
        row_for(tb[tb["gait"] == "Pacer"],                     "Pacer"),
        row_for(tb[tb["gait_bookSize"] == "Pacer_LOW"],        "Pacer_LOW"),
        row_for(tb[tb["gait_bookSize"] == "Pacer_MEDIUM"],     "Pacer_MEDIUM"),
        row_for(tb[tb["gait_bookSize"] == "Pacer_HIGH"],       "Pacer_HIGH"),
        row_for(tb[tb["gait"] == "Trotter"],                   "Trotter"),
        row_for(tb[tb["gait_bookSize"] == "Trotter_LOW"],      "Trotter_LOW"),
        row_for(tb[tb["gait_bookSize"] == "Trotter_MEDIUM"],   "Trotter_MEDIUM"),
        row_for(tb[tb["gait_bookSize"] == "Trotter_HIGH"],     "Trotter_HIGH"),
    ]
    out_df = pd.DataFrame(rows, columns=["Subgroup", "N", *METRICS])

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_df.to_csv(out_path, index=False)
    print(f"[subgroup summary] wrote {out_path} ({len(out_df)} rows)")


if __name__ == "__main__":
    main()
