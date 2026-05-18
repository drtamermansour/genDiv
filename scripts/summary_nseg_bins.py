import argparse
import pandas as pd
import numpy as np
from utils import format_stats

# ================= CLI SETUP =================
parser = argparse.ArgumentParser(
    description="Per-subgroup summary of ROH segment counts binned by length (Mb)."
)
parser.add_argument("-i", "--roh-file", required=True,
                    help="Per-segment bcftools ROH file with RG-rows, e.g. roh.L3.${rg}.txt.")
parser.add_argument("-f", "--factor-file", required=True,
                    help="Per-IID factor TSV with header; IID in column 1, grouping factor in "
                         "the last column (same shape produced for summary_roh.py, e.g. "
                         "roh.L3_Froh_gait.txt).")
parser.add_argument("-o", "--output-file", required=True,
                    help="Path to write the output summary CSV.")
args = parser.parse_args()

BIN_EDGES = [1.0, 3.0, 5.0, 10.0, np.inf]
BIN_LABELS = ["NSEG_1to3", "NSEG_3to5", "NSEG_5to10", "NSEG_more10"]

# 1. Load factor file -> full IID universe + factor mapping
factor_df = pd.read_csv(args.factor_file, sep='\t')
iid_col = factor_df.columns[0]
factor_col = factor_df.columns[-1]
factor_df = factor_df[[iid_col, factor_col]].rename(columns={iid_col: "IID"})

# 2. Load per-segment ROH file (skip "# RG ..." header; data rows start with "RG\t")
seg_df = pd.read_csv(
    args.roh_file, sep='\t', comment='#', header=None,
    usecols=[1, 5], names=["IID", "SIZE_BP"],
)
seg_df["SIZE_MB"] = seg_df["SIZE_BP"] / 1e6
seg_df["bin"] = pd.cut(seg_df["SIZE_MB"], bins=BIN_EDGES, labels=BIN_LABELS, right=False)

# 3. Pivot to per-IID bin counts
counts = (
    seg_df.dropna(subset=["bin"])
          .groupby(["IID", "bin"], observed=False)
          .size()
          .unstack(fill_value=0)
          .reindex(columns=BIN_LABELS, fill_value=0)
)

# 4. Left-join onto the full IID universe so zero-segment samples appear as 0s
df = factor_df.merge(counts, left_on="IID", right_index=True, how="left")
df[BIN_LABELS] = df[BIN_LABELS].fillna(0).astype(int)

# 5. Whole-population + per-subgroup stats
pop_stats = format_stats(df, BIN_LABELS)
pop_stats.update({"Subgroup": "Whole Population", "N": len(df)})
results_list = [pop_stats]

for name, group_df in df.groupby(factor_col):
    sub_stats = format_stats(group_df, BIN_LABELS)
    sub_stats.update({"Subgroup": name, "N": len(group_df)})
    results_list.append(sub_stats)

# 6. Write CSV
output_df = pd.DataFrame(results_list)
final_cols = ["Subgroup", "N"] + BIN_LABELS
output_df = output_df[final_cols]
output_df.to_csv(args.output_file, index=False)
print(f"Summary saved to {args.output_file}")
