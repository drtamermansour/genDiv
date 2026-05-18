"""Per-individual ROH_common score for one group.

For each individual i and each ROH length class c (plus genome-wide), compute

    ROH_common_{i,c} = mean over w in W_{i,c} of (n_w - 1) / (N - 1)

where W_{i,c} is the set of 100 kb windows overlapping any of individual i's
class-c ROH segments and n_w is the population count of distinct individuals
whose ROH overlaps window w. The (n_w - 1)/(N - 1) form is the exact
leave-one-out estimator: since i is by construction one of the n_w contributors
at every window in W_{i,c}, the proposal's f_w - 1/N approximation is biased
upward (~1/N relative error). The exact form costs nothing extra.

Inputs:
  * landscape directory produced by roh_common_landscape.py (windows.${rg}.bed,
    landscape.${rg}.tsv, landscape.${rg}.${class}.tsv, segments.${rg}.${class}.bed,
    N.${rg}.txt)
  * the per-segment BED used to build it (passed in for IID enumeration and
    the genome-wide pass)
"""

import argparse
import subprocess
from io import StringIO
from pathlib import Path

import numpy as np
import pandas as pd

LENGTH_BIN_LABEL = ("1to3", "3to5", "5to10", "more10")


def run_bedtools(args):
    return subprocess.run(args, capture_output=True, text=True, check=True).stdout


def intersect_window_iid_pairs(windows_bed, subset_bed):
    """Return distinct (chrom, start, end, IID) pairs where the IID's subset_bed
    intervals overlap the window. Empty DataFrame if no overlap."""
    out = run_bedtools(
        ["bedtools", "intersect", "-a", str(windows_bed), "-b", str(subset_bed),
         "-wa", "-wb"]
    )
    if not out.strip():
        return pd.DataFrame(columns=["chrom", "start", "end", "IID"])
    cols = ["chrom", "start", "end", "s_chrom", "s_start", "s_end", "IID"]
    df = pd.read_csv(StringIO(out), sep="\t", header=None, names=cols,
                     dtype={"chrom": str, "s_chrom": str, "IID": str})
    return df[["chrom", "start", "end", "IID"]].drop_duplicates()


def score_one_pass(windows_bed, subset_bed, landscape_df, N):
    """Return per-IID DataFrame with columns IID, ROH_common, n_windows for one
    BED subset (genome-wide or one class). IIDs absent from the subset BED do
    not appear in the result."""
    pairs = intersect_window_iid_pairs(windows_bed, subset_bed)
    if pairs.empty:
        return pd.DataFrame(columns=["IID", "ROH_common", "n_windows"])
    merged = pairs.merge(
        landscape_df[["chrom", "start", "end", "n_w"]],
        on=["chrom", "start", "end"], how="left",
    )
    # Exact LOO: every window in W_{i,c} has i as a contributor, so n_w >= 1.
    # If a join miss yields NaN n_w we treat it as a contract violation.
    if merged["n_w"].isna().any():
        missing = merged[merged["n_w"].isna()]
        raise RuntimeError(
            "Landscape lookup miss for windows: "
            f"{missing[['chrom', 'start', 'end']].head().to_dict('records')}"
        )
    merged["f_corr"] = (merged["n_w"] - 1) / (N - 1)
    agg = (merged.groupby("IID")
                 .agg(ROH_common=("f_corr", "mean"),
                      n_windows=("f_corr", "size"))
                 .reset_index())
    return agg


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--landscape-dir", required=True,
                        help="Directory produced by roh_common_landscape.py.")
    parser.add_argument("--segments-bed", required=True,
                        help="Per-segment L3 BED (same file used by the "
                             "landscape script).")
    parser.add_argument("--rg", required=True,
                        help="Group label.")
    parser.add_argument("--out", required=True,
                        help="Output TSV path.")
    args = parser.parse_args()

    land_dir = Path(args.landscape_dir)
    windows_bed = land_dir / f"windows.{args.rg}.bed"
    N = int((land_dir / f"N.{args.rg}.txt").read_text().strip())

    # Full IID universe from the per-segment BED (matches the N used to build
    # the landscape).
    seg_df = pd.read_csv(
        args.segments_bed, sep="\t", header=None,
        names=["chrom", "start", "end", "IID"],
        dtype={"chrom": str, "start": np.int64, "end": np.int64, "IID": str},
    )
    all_iids = pd.DataFrame({"IID": sorted(seg_df["IID"].unique())})
    assert len(all_iids) == N, f"IID count mismatch: bed has {len(all_iids)}, N file says {N}"

    out_df = all_iids.copy()

    # Genome-wide pass.
    landscape_gw = pd.read_csv(
        land_dir / f"landscape.{args.rg}.tsv", sep="\t",
        dtype={"chrom": str, "start": np.int64, "end": np.int64,
               "n_w": np.int64, "f_w": np.float64},
    )
    gw_scores = score_one_pass(windows_bed, args.segments_bed, landscape_gw, N)
    gw_scores = gw_scores.rename(columns={"ROH_common": "ROH_common",
                                          "n_windows": "n_windows_total"})
    out_df = out_df.merge(gw_scores, on="IID", how="left")

    # Per-class passes.
    for label in LENGTH_BIN_LABEL:
        class_subset_bed = land_dir / f"segments.{args.rg}.{label}.bed"
        class_landscape = pd.read_csv(
            land_dir / f"landscape.{args.rg}.{label}.tsv", sep="\t",
            dtype={"chrom": str, "start": np.int64, "end": np.int64,
                   "n_w": np.int64, "f_w": np.float64},
        )
        scores = score_one_pass(windows_bed, class_subset_bed, class_landscape, N)
        scores = scores.rename(columns={
            "ROH_common": f"ROH_common_{label}",
            "n_windows":  f"n_windows_{label}",
        })
        out_df = out_df.merge(scores, on="IID", how="left")

    # Zero-fill window counts (no segments => 0 windows touched); leave
    # ROH_common columns as NA for samples with empty class membership.
    for label in LENGTH_BIN_LABEL:
        out_df[f"n_windows_{label}"] = out_df[f"n_windows_{label}"].fillna(0).astype(np.int64)
    out_df["n_windows_total"] = out_df["n_windows_total"].fillna(0).astype(np.int64)

    col_order = ["IID", "ROH_common"] + [f"ROH_common_{l}" for l in LENGTH_BIN_LABEL] + \
                ["n_windows_total"] + [f"n_windows_{l}" for l in LENGTH_BIN_LABEL]
    out_df = out_df[col_order]
    out_df.to_csv(args.out, sep="\t", index=False, float_format="%.6f", na_rep="NA")
    print(f"[roh_common_individual] {args.rg}: scored {len(out_df)} samples; "
          f"wrote {args.out}")


if __name__ == "__main__":
    main()
