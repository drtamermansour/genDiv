import argparse
import pandas as pd
import numpy as np
import os
import sys

# ================= CLI SETUP =================
parser = argparse.ArgumentParser(description="Generate a summary table of heterozygosity and inbreeding coefficient.")
parser.add_argument("-i", "--het_file", type=str, required=True, help="Path to input PLINK2 het TSV with extracolumn for group ID")
parser.add_argument("-o", "--output_file", type=str, required=True, help="Path to write the output summary stats")
args = parser.parse_args()

INPUT_HET = args.het_file
OUTPUT_FILE = args.output_file


# 1. Read the tab-separated file
# Column 5, 6, 8, and 9 correspond to indices 4, 5, 7, and 8
df = pd.read_csv(INPUT_HET, sep='\t')

# Get column names to avoid .iloc assignment issues
col_4_name = df.columns[4]
col_5_name = df.columns[5]
col_6_name = df.columns[6] # Denominator column

# Safe implementation of logic
# Use .replace(0, np.nan) to prevent division by zero errors (inf)
denominator = df[col_6_name].replace(0, np.nan)

# calculate observed and expected Heterozygosity as a fraction
df[col_4_name] = df[col_4_name] / denominator
df[col_5_name] = df[col_5_name] / denominator

# Handle resulting infinities (if any) by turning them into NaN
df.replace([np.inf, -np.inf], np.nan, inplace=True)

# Get actual column names from the header
cols_idx = [4, 5, 7]
target_cols = [df.columns[i] for i in cols_idx]
factor_col = df.columns[8]  # Column 9

# Helper function to format Mean +/- SD
def format_stats(group_df, columns):
    results = {}
    for col in columns:
        m = group_df[col].mean()
        s = group_df[col].std()
        results[col] = f"{m:.2f} +/- {s:.2f}"
    return results

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
# Reorder columns to: Subgroup, N, [Col 5 Name], [Col 6 Name], [Col 8 Name]
final_cols = ['Subgroup', 'N'] + target_cols
output_df = output_df[final_cols]

output_df.to_csv(OUTPUT_FILE, index=False)
print(f"Summary saved to {OUTPUT_FILE}")

