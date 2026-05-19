"""Confirm that PC outliers correspond to genuine high-kinship clusters.

For each principal component in the user-selected range (default PC2-PC4),
identify the top-N individuals at each tail (most positive and most
negative PC loadings) and compute the mean off-diagonal pairwise GRM
kinship among those tail individuals. Compare against the cohort-wide
mean (≈ 0 by GRM construction).

The reasoning: PLINK2 --pca decomposes the same additive-genetic
relationship matrix that --make-rel 'square' writes out. If a low-order
PC reflects "hidden familial structure" (as is often the case for
PC2-PC4 in a closed studbook), the individuals at the tails of that PC
should be enriched for high pairwise kinship with each other relative
to a random subset of the cohort. This script makes that confirmation
numerical.

Inputs:
  --eigenvec: PLINK2 .eigenvec file (#FID IID PC1 PC2 ... header).
  --rel:      PLINK2 .rel (square N x N) file from --make-rel 'square'.
  --rel-id:   .rel.id file (FID, IID; matrix row/column order).
  --pcs:      Comma-separated PC names to analyse (default "PC2,PC3,PC4").
  --top-n:    Number of individuals at each tail per PC (default 10).
  --out:      Output CSV path.
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


def mean_pairwise(grm, indices):
    """Mean off-diagonal pairwise GRM value among the given row/col indices."""
    if len(indices) < 2:
        return float("nan")
    sub = grm[np.ix_(indices, indices)]
    # Upper-triangle off-diagonal
    iu = np.triu_indices(len(indices), k=1)
    return float(sub[iu].mean())


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--eigenvec", required=True)
    ap.add_argument("--rel", required=True)
    ap.add_argument("--rel-id", required=True)
    ap.add_argument("--pcs", default="PC2,PC3,PC4",
                    help="Comma-separated PC names (default PC2,PC3,PC4).")
    ap.add_argument("--top-n", type=int, default=10,
                    help="Number of individuals at each tail per PC "
                         "(default 10).")
    ap.add_argument("--out", required=True, help="Output CSV path.")
    args = ap.parse_args()

    # Load matrix sample order.
    ids_df = pd.read_csv(args.rel_id, sep=r"\s+", header=None, comment="#",
                         names=["FID", "IID"], dtype=str)
    rel_ids = ids_df["IID"].tolist()
    iid_to_row = {iid: i for i, iid in enumerate(rel_ids)}
    N = len(rel_ids)

    # Load GRM.
    grm = np.loadtxt(args.rel)
    if grm.shape != (N, N):
        raise SystemExit(
            f"GRM shape {grm.shape} does not match .rel.id row count {N}"
        )

    # Cohort-wide off-diagonal mean.
    iu = np.triu_indices(N, k=1)
    cohort_mean = float(grm[iu].mean())

    # Load eigenvec.
    eig = pd.read_csv(args.eigenvec, sep=r"\s+")
    if eig.columns[0].startswith("#"):
        eig = eig.rename(columns={eig.columns[0]: eig.columns[0].lstrip("#")})
    if "IID" not in eig.columns:
        raise SystemExit(
            f"--eigenvec must contain an IID column; got {list(eig.columns)}"
        )

    pcs = [p.strip() for p in args.pcs.split(",") if p.strip()]
    missing = [p for p in pcs if p not in eig.columns]
    if missing:
        raise SystemExit(f"--eigenvec is missing PC columns: {missing}")

    rows = []
    for pc in pcs:
        sorted_iids = eig.sort_values(pc)["IID"].tolist()
        bottom = sorted_iids[: args.top_n]
        top = sorted_iids[-args.top_n:]

        for tail_label, iids in (("bottom", bottom), ("top", top)):
            indices = [iid_to_row[i] for i in iids if i in iid_to_row]
            mean_k = mean_pairwise(grm, indices)
            excess = (mean_k - cohort_mean) if not np.isnan(mean_k) else float("nan")
            rows.append({
                "PC": pc,
                "tail": tail_label,
                "n_outliers": len(indices),
                "mean_within_cluster_kinship": (
                    f"{mean_k:.4f}" if not np.isnan(mean_k) else "NA"
                ),
                "cohort_mean_kinship": f"{cohort_mean:.4f}",
                "excess_over_cohort": (
                    f"{excess:+.4f}" if not np.isnan(excess) else "NA"
                ),
                "outlier_IIDs": ";".join(iids),
            })

    out_df = pd.DataFrame(rows)
    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_df.to_csv(out_path, index=False)
    print(f"[pc outlier kinship] wrote {out_path} ({len(out_df)} rows)")


if __name__ == "__main__":
    main()
