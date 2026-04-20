import argparse
import pandas as pd
import numpy as np
import sys
from utils import format_stats

# ================= CLI SETUP =================
parser = argparse.ArgumentParser(description="Generate a summary table of ROH analysis.")
parser.add_argument("-i", "--roh_file", type=str, required=True,
                    help="Path to input ROH analysis TSV with extra column for group ID")
parser.add_argument("-o", "--output_file", type=str, required=True,
                    help="Path to write the output summary stats")
parser.add_argument("-n", "--n-numeric-cols", type=int, default=3,
                    help="Number of numeric columns after IID to summarise (columns 2..1+N). "
                         "The next column (index N+1) is treated as the grouping factor. "
                         "Default: 3 (matches the legacy summary_roh). Use 4 to match legacy summary_roh_v2.")
args = parser.parse_args()

INPUT_ROH = args.roh_file
OUTPUT_FILE = args.output_file
N_COLS = args.n_numeric_cols

# 1. Read the tab-separated file
df = pd.read_csv(INPUT_ROH, sep='\t')

# Handle resulting infinities (if any) by turning them into NaN
df.replace([np.inf, -np.inf], np.nan, inplace=True)

# Numeric columns are 1..N_COLS (IID is column 0); grouping factor is column N_COLS+1.
cols_idx = list(range(1, 1 + N_COLS))
target_cols = [df.columns[i] for i in cols_idx]
factor_col = df.columns[1 + N_COLS]

# 2. Calculate Whole Population Statistics
pop_stats = format_stats(df, target_cols)
pop_stats.update({'Subgroup': 'Whole Population', 'N': len(df)})
results_list = [pop_stats]

# 3. Calculate Subpopulation Statistics
grouped = df.groupby(factor_col)
for name, group_df in grouped:
    sub_stats = format_stats(group_df, target_cols)
    sub_stats.update({'Subgroup': name, 'N': len(group_df)})
    results_list.append(sub_stats)

# 4. Create final DataFrame and save to CSV
output_df = pd.DataFrame(results_list)
final_cols = ['Subgroup', 'N'] + target_cols
output_df = output_df[final_cols]

output_df.to_csv(OUTPUT_FILE, index=False)
print(f"Summary saved to {OUTPUT_FILE}")
