"""ROH-based effective population size (N_e) by length class.

For each group (typically wholePop, Trotter, Pacer) and each ROH length
class (1-3, 3-5, 5-10, >=10 Mb), compute the within-group autozygous
fraction contributed by class-c segments and translate it into an
effective population size estimate at the corresponding ancestral
generation depth.

Method.
- For each L3-filtered ROH segment, the length L_Mb is computed.
- Segments are partitioned into length classes by --length-bins.
- For each class c, the within-group autozygous fraction is
      F_ROH(c) = (sum of class-c segment lengths in bp across animals)
                 / (N_animals * effective_autosomal_genome_length)
- The mean segment length in the class, converted to Morgans via the
  EquCab3 genome-wide recombination rate (--cm-per-mb, default 1.16
  cM/Mb; Beeson et al. 2020), gives both the corresponding generation
  depth t = 50 / L_cM and the per-class recombination distance L
  (Morgans) used in the N_e formula.
- N_e at depth t is then estimated as
      N_e(t) ≈ 1 / (4 * L_Morgan * F_ROH(c))
  (Hayes et al. 2003; Marras et al. 2015).

Caveats. The formula assumes a panmictic population, that autozygous
tracts in class c derive from a common ancestor at generation
t = 50/L_cM, and that recombination has not yet substantially broken
those tracts. In a closed studbook with popular sires (the very
demography we are studying), both assumptions hold imperfectly.
Estimates should be read as order-of-magnitude indicators of recent
effective population size, not as precise point values.

Inputs.
  --group LABEL:SEGMENTS_BED:SAMPLES_FILE
      Repeat per group. SEGMENTS_BED has columns chrom, start, end, IID
      (the L3 BED). SAMPLES_FILE has one IID per line, or a PLINK-style
      two-column FID IID layout (one row per sample).
  --autosomal-length-file
      File containing a single integer: effective autosomal length in bp
      (effective_autosomal_genome_length.txt produced by shared.sh).
  --cm-per-mb (default 1.16; Beeson et al. 2020 EquCab3 mean)
  --length-bins (default "1,3,5,10" -> [1,3), [3,5), [5,10), [10,inf))
  --out CSV output path
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

LENGTH_BIN_LABEL = ("1to3", "3to5", "5to10", "more10")
DEFAULT_CM_PER_MB = 1.16  # Beeson et al. 2020, EquCab3 genome-wide mean


def parse_length_bins(spec):
    edges = [float(x) for x in spec.split(",")]
    if len(edges) < 2:
        raise ValueError(f"Need at least 2 edges, got: {spec!r}")
    edges_extended = edges + [float("inf")]
    if len(edges_extended) - 1 != len(LENGTH_BIN_LABEL):
        raise ValueError(
            f"This script ships with {len(LENGTH_BIN_LABEL)} class labels "
            f"{LENGTH_BIN_LABEL}; got {len(edges_extended) - 1} bins from {spec!r}"
        )
    return list(zip(edges_extended[:-1], edges_extended[1:], LENGTH_BIN_LABEL))


def count_samples(path):
    with open(path) as f:
        return sum(1 for line in f if line.strip())


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--group", action="append", required=True,
                    help="Repeatable: LABEL:SEGMENTS_BED:SAMPLES_FILE.")
    ap.add_argument("--autosomal-length-file", required=True,
                    help="File with a single integer: effective autosomal "
                         "length in bp.")
    ap.add_argument("--cm-per-mb", type=float, default=DEFAULT_CM_PER_MB,
                    help="Recombination rate in cM/Mb (default 1.16; "
                         "Beeson et al. 2020 EquCab3 genome-wide mean).")
    ap.add_argument("--length-bins", default="1,3,5,10",
                    help="Comma-separated Mb edges (default '1,3,5,10').")
    ap.add_argument("--out", required=True, help="Output CSV path.")
    args = ap.parse_args()

    aut_len_bp = int(open(args.autosomal_length_file).read().strip())
    cm_per_mb = args.cm_per_mb
    bins = parse_length_bins(args.length_bins)

    rows = []
    for spec in args.group:
        parts = spec.split(":")
        if len(parts) != 3:
            raise SystemExit(
                f"--group expects LABEL:SEGMENTS_BED:SAMPLES_FILE; got {spec!r}"
            )
        label, seg_bed, sample_file = parts
        n_animals = count_samples(sample_file)

        seg_df = pd.read_csv(
            seg_bed, sep=r"\s+", header=None,
            names=["chrom", "start", "end", "IID"],
            dtype={"chrom": str, "start": np.int64, "end": np.int64, "IID": str},
        )
        seg_df["size_mb"] = (seg_df["end"] - seg_df["start"]) / 1e6

        # Per-class rows.
        class_f_roh_sum = 0.0
        class_n_seg_sum = 0
        for low_mb, high_mb, bin_label in bins:
            if np.isinf(high_mb):
                mask = seg_df["size_mb"] >= low_mb
                pretty_class = f">{int(low_mb)} Mb"
            else:
                mask = (seg_df["size_mb"] >= low_mb) & (seg_df["size_mb"] < high_mb)
                pretty_class = f"{int(low_mb)}-{int(high_mb)} Mb"

            class_seg = seg_df[mask]
            n_seg = int(len(class_seg))
            sum_mb = float(class_seg["size_mb"].sum()) if n_seg else 0.0
            mean_mb = (sum_mb / n_seg) if n_seg else float("nan")

            f_roh_class = ((sum_mb * 1e6) / (n_animals * aut_len_bp)
                           if n_animals else float("nan"))

            if n_seg and mean_mb > 0 and f_roh_class and f_roh_class > 0:
                mean_cm = mean_mb * cm_per_mb
                mean_morgan = mean_cm / 100.0
                t_gen = 50.0 / mean_cm
                ne = 1.0 / (4.0 * mean_morgan * f_roh_class)
            else:
                mean_cm = float("nan")
                t_gen = float("nan")
                ne = float("nan")

            rows.append({
                "Group": label,
                "Length_class": pretty_class,
                "N_animals": n_animals,
                "N_segments": n_seg,
                "Mean_segment_length_Mb": f"{mean_mb:.2f}" if not np.isnan(mean_mb) else "NA",
                "Mean_segment_length_cM": f"{mean_cm:.2f}" if not np.isnan(mean_cm) else "NA",
                "Generation_depth_t": f"{t_gen:.1f}" if not np.isnan(t_gen) else "NA",
                "F_ROH_class": f"{f_roh_class:.4f}" if not np.isnan(f_roh_class) else "NA",
                "Ne": f"{ne:.0f}" if not np.isnan(ne) else "NA",
            })

            if not np.isnan(f_roh_class):
                class_f_roh_sum += f_roh_class
            class_n_seg_sum += n_seg

        # Group-level sanity-check row: F_ROH computed independently across
        # the full L3 BED versus the sum of per-class F_ROH(c). The two
        # should match by construction (the length classes partition all
        # L3 segments); any mismatch flags a binning or filter bug.
        n_seg_total_bed = int(len(seg_df))
        f_roh_total_bed = (float(seg_df["size_mb"].sum()) * 1e6
                           / (n_animals * aut_len_bp)
                           if n_animals else float("nan"))
        match_flag = (
            "OK"
            if (not np.isnan(f_roh_total_bed)
                and abs(class_f_roh_sum - f_roh_total_bed) < 1e-6
                and class_n_seg_sum == n_seg_total_bed)
            else "MISMATCH"
        )
        rows.append({
            "Group": label,
            "Length_class": (
                f"Total (>=1 Mb, sanity check, Sum_class vs BED-total: "
                f"{match_flag})"
            ),
            "N_animals": n_animals,
            "N_segments": n_seg_total_bed,
            "Mean_segment_length_Mb": "NA",
            "Mean_segment_length_cM": "NA",
            "Generation_depth_t": "NA",
            "F_ROH_class": (
                f"{f_roh_total_bed:.4f}"
                if not np.isnan(f_roh_total_bed) else "NA"
            ),
            "Ne": "NA",
        })

    out_df = pd.DataFrame(rows)
    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_df.to_csv(out_path, index=False)
    print(f"[ne roh] wrote {out_path} ({len(out_df)} rows)")


if __name__ == "__main__":
    main()
