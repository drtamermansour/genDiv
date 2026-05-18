"""Build the per-window population ROH frequency landscape (f_w) for one group.

Implements the landscape step of the ROH_common metric (see manuscript/ROH_common.md):

    f_w = n_w / N

where n_w is the number of *distinct individuals* whose L3 ROH segments overlap
window w (any-overlap semantics) and N is the group sample size.

Two passes are run:
  * genome-wide: all L3 segments, no length filtering.
  * length-stratified: one landscape per ROH length class (default
    1-3, 3-5, 5-10, >10 Mb; configurable via --length-bins).

All interval arithmetic is delegated to `bedtools` via subprocess to keep
half-open BED semantics consistent with the rest of the pipeline. Pandas is
used only for distinct-IID counting and TSV emission.
"""

import argparse
import subprocess
import sys
from io import StringIO
from pathlib import Path

import numpy as np
import pandas as pd

LENGTH_BIN_LABEL = ("1to3", "3to5", "5to10", "more10")


def run_bedtools(args, stdout_path=None):
    """Run a bedtools command, returning stdout as text. Raises on non-zero exit."""
    if stdout_path is not None:
        with open(stdout_path, "w") as fh:
            subprocess.run(args, stdout=fh, check=True)
        return None
    result = subprocess.run(args, capture_output=True, text=True, check=True)
    return result.stdout


def parse_length_bins(spec):
    """'1,3,5,10' -> [1.0, 3.0, 5.0, 10.0, inf] in Mb, with labels."""
    edges = [float(x) for x in spec.split(",")]
    if len(edges) < 2:
        raise ValueError(f"Need at least two length-bin edges, got: {spec!r}")
    edges_extended = edges + [float("inf")]
    if len(edges_extended) - 1 != len(LENGTH_BIN_LABEL):
        raise ValueError(
            f"This script ships with {len(LENGTH_BIN_LABEL)} class labels "
            f"{LENGTH_BIN_LABEL}; got {len(edges_extended) - 1} bins from {spec!r}"
        )
    return list(zip(edges_extended[:-1], edges_extended[1:], LENGTH_BIN_LABEL))


def load_segments_bed(segments_path):
    """Load the per-segment L3 BED: chrom, start, end, IID. Compute size_mb."""
    df = pd.read_csv(
        segments_path, sep="\t", header=None,
        names=["chrom", "start", "end", "IID"],
        dtype={"chrom": str, "start": np.int64, "end": np.int64, "IID": str},
    )
    df["size_mb"] = (df["end"] - df["start"]) / 1e6
    return df


def write_bed(df, path):
    """Write chrom/start/end/IID BED, sorted by chrom/start."""
    df_sorted = df.sort_values(["chrom", "start", "end"])
    df_sorted[["chrom", "start", "end", "IID"]].to_csv(
        path, sep="\t", header=False, index=False
    )


def count_distinct_iids_per_window(windows_bed, subset_bed):
    """For each window, count the number of distinct IIDs whose intervals overlap it.

    Uses `bedtools intersect -wa -wb` so we get window + segment records, then
    dedups (window, IID) pairs and groups by window to get n_w. Windows with
    zero overlapping segments are not in the output and must be zero-filled
    by the caller against the full windows list.
    """
    out = run_bedtools(
        ["bedtools", "intersect", "-a", str(windows_bed), "-b", str(subset_bed),
         "-wa", "-wb"]
    )
    if not out.strip():
        return pd.DataFrame(columns=["chrom", "start", "end", "n_w"])
    cols = ["w_chrom", "w_start", "w_end", "s_chrom", "s_start", "s_end", "IID"]
    df = pd.read_csv(StringIO(out), sep="\t", header=None, names=cols,
                     dtype={"w_chrom": str, "s_chrom": str, "IID": str})
    df = df.drop_duplicates(["w_chrom", "w_start", "w_end", "IID"])
    n_w = (
        df.groupby(["w_chrom", "w_start", "w_end"]).size()
          .reset_index(name="n_w")
          .rename(columns={"w_chrom": "chrom", "w_start": "start", "w_end": "end"})
    )
    return n_w


def build_landscape(windows_df, subset_bed, N, windows_bed_path):
    """Run the count + zero-fill + f_w computation for one BED subset."""
    n_w = count_distinct_iids_per_window(windows_bed_path, subset_bed)
    landscape = windows_df.merge(n_w, on=["chrom", "start", "end"], how="left")
    landscape["n_w"] = landscape["n_w"].fillna(0).astype(np.int64)
    landscape["f_w"] = landscape["n_w"] / N
    return landscape[["chrom", "start", "end", "n_w", "f_w"]]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--segments-bed", required=True,
                        help="Per-segment L3 BED: chrom/start/end/IID "
                             "(e.g., divStats/roh.L3.${rg}.bed).")
    parser.add_argument("--autosomes", required=True,
                        help="Tab-separated chrom-size file used by bedtools "
                             "(e.g., divStats/autosomes.genome).")
    parser.add_argument("--window-kb", type=int, default=100,
                        help="Window size in kb (default: 100).")
    parser.add_argument("--length-bins", default="1,3,5,10",
                        help="Comma-separated Mb edges (default '1,3,5,10' -> "
                             "[1,3), [3,5), [5,10), [10,inf)).")
    parser.add_argument("--out-dir", required=True,
                        help="Output directory; created if absent.")
    parser.add_argument("--rg", required=True,
                        help="Group label used in output filenames "
                             "(e.g. wholePop, Trotter, Pacer, Trotter_LOW).")
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    bins = parse_length_bins(args.length_bins)
    window_bp = args.window_kb * 1000

    # 1. Generate the fixed window tiling.
    windows_bed_path = out_dir / f"windows.{args.rg}.bed"
    run_bedtools(
        ["bedtools", "makewindows", "-g", args.autosomes, "-w", str(window_bp)],
        stdout_path=str(windows_bed_path),
    )
    windows_df = pd.read_csv(
        windows_bed_path, sep="\t", header=None,
        names=["chrom", "start", "end"],
        dtype={"chrom": str, "start": np.int64, "end": np.int64},
    )

    # 2. Load the per-segment BED and determine N.
    seg_df = load_segments_bed(args.segments_bed)
    N = seg_df["IID"].nunique()
    if N == 0:
        sys.exit(f"ERROR: zero distinct IIDs in {args.segments_bed}; cannot proceed.")
    (out_dir / f"N.{args.rg}.txt").write_text(f"{N}\n")

    # 3. Genome-wide landscape (all L3 segments).
    landscape_gw = build_landscape(
        windows_df, args.segments_bed, N, windows_bed_path
    )
    landscape_gw.to_csv(
        out_dir / f"landscape.{args.rg}.tsv", sep="\t", index=False,
        float_format="%.6f",
    )

    # 4. Length-stratified landscapes.
    for low_mb, high_mb, label in bins:
        if np.isinf(high_mb):
            mask = seg_df["size_mb"] >= low_mb
        else:
            mask = (seg_df["size_mb"] >= low_mb) & (seg_df["size_mb"] < high_mb)
        subset = seg_df[mask]
        subset_bed = out_dir / f"segments.{args.rg}.{label}.bed"
        write_bed(subset, subset_bed)

        if subset.empty:
            # Empty subset -> all n_w = 0, f_w = 0. Emit a zero-filled landscape
            # so downstream readers always find the file.
            landscape = windows_df.copy()
            landscape["n_w"] = 0
            landscape["f_w"] = 0.0
        else:
            landscape = build_landscape(
                windows_df, subset_bed, N, windows_bed_path
            )
        landscape.to_csv(
            out_dir / f"landscape.{args.rg}.{label}.tsv", sep="\t", index=False,
            float_format="%.6f",
        )

    print(f"[roh_common_landscape] {args.rg}: N={N}, "
          f"{len(windows_df)} windows, {len(bins)} length classes.")


if __name__ == "__main__":
    main()
