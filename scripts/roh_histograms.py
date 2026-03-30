"""
Unified ROH histogram script.

Usage:
    python roh_histograms.py --metric ratio  <roh_shared_file> <froh_file> <output_prefix>
    python roh_histograms.py --metric shared <roh_shared_file> <froh_file> <output_prefix>

--metric ratio:  plots ROH_shared / F_ROH ratio per individual (normalized inbreeding)
--metric shared: plots Percent_of_Consensus_ROH per individual (raw sharing)
"""

import argparse
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
from utils import BOOK_SIZE_ORDER, BOOK_SIZE_COLORS


def load_and_merge(file1_path, file2_path):
    """Load and merge the ROH-shared and F_ROH data files on IID."""
    try:
        df1 = pd.read_csv(file1_path, sep=r'\s+')
        df2 = pd.read_csv(file2_path, sep=r'\s+')
    except FileNotFoundError as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)

    df = pd.merge(df1, df2, on='IID', how='inner')
    df['gait_type'] = df['gait'].str.split('_').str[0]
    df['book_size'] = df['gait'].str.split('_').str[1]
    df['book_size'] = pd.Categorical(df['book_size'], categories=BOOK_SIZE_ORDER, ordered=True)
    return df


def add_metric(df, metric):
    """Add the target metric column to the dataframe."""
    if metric == 'ratio':
        epsilon = 1e-10
        df['_metric'] = (df['Percent_of_Consensus_ROH'] / 100) / (df['F_ROH'] + epsilon)
        col_label = 'ROH_shared / F_ROH'
        title_suffix = 'ROH_shared / F_ROH Ratio'
    else:  # 'shared'
        df['_metric'] = df['Percent_of_Consensus_ROH']
        col_label = 'ROH_shared (%)'
        title_suffix = 'ROH_shared (%)'
    return df, col_label, title_suffix


def _plot_gait(ax, gait_data, gait_name, col_label, title_suffix, plot_type):
    """Plot histogram or density for one gait group."""
    hist_data, labels = [], []
    for bs in BOOK_SIZE_ORDER:
        subset = gait_data[gait_data['book_size'] == bs]
        if len(subset) > 0:
            hist_data.append(subset['_metric'])
            labels.append(f'{bs} (n={len(subset)})')

    if not hist_data:
        return

    if plot_type == 'histogram':
        n_bins = min(20, int(np.sqrt(len(gait_data))))
        ax.hist(hist_data, bins=n_bins, stacked=False,
                color=[BOOK_SIZE_COLORS[bs] for bs in BOOK_SIZE_ORDER[:len(hist_data)]],
                alpha=0.7, edgecolor='black', linewidth=0.8, label=labels)
        # Mean lines
        for bs, color in zip(BOOK_SIZE_ORDER[:len(hist_data)],
                              [BOOK_SIZE_COLORS[bs] for bs in BOOK_SIZE_ORDER[:len(hist_data)]]):
            subset = gait_data[gait_data['book_size'] == bs]
            if len(subset) > 0:
                mean_val = subset['_metric'].mean()
                ax.axvline(mean_val, color=color, linestyle='--', linewidth=2,
                           alpha=0.8, label=f'{bs} Mean: {mean_val:.1f}')
        ax.set_ylabel('Frequency', fontsize=12)
        ax.legend(fontsize=10, title='Book Size', title_fontsize=11)

        if not gait_data.empty:
            stats_text = (f'Mean: {gait_data["_metric"].mean():.1f}\n'
                          f'Std:  {gait_data["_metric"].std():.1f}\n'
                          f'Median: {gait_data["_metric"].median():.1f}\n'
                          f'Min: {gait_data["_metric"].min():.1f}, '
                          f'Max: {gait_data["_metric"].max():.1f}')
            ax.text(0.02, 0.98, stats_text, transform=ax.transAxes,
                    fontsize=9, verticalalignment='top',
                    bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))

    else:  # density
        for bs in BOOK_SIZE_ORDER:
            subset = gait_data[gait_data['book_size'] == bs]
            if len(subset) > 0:
                ax.hist(subset['_metric'], bins=15, density=True, alpha=0.3,
                        color=BOOK_SIZE_COLORS[bs], edgecolor='black')
                kde = gaussian_kde(subset['_metric'])
                x_range = np.linspace(subset['_metric'].min() * 0.9,
                                      subset['_metric'].max() * 1.1, 200)
                ax.plot(x_range, kde(x_range), color=BOOK_SIZE_COLORS[bs],
                        linewidth=2, label=f'{bs} (n={len(subset)})')
        ax.set_ylabel('Density', fontsize=12)
        ax.legend(fontsize=10)

    ax.set_title(f'{gait_name} - Distribution of {title_suffix}',
                 fontsize=14, fontweight='bold', pad=15)
    ax.set_xlabel(col_label, fontsize=12)
    ax.grid(True, alpha=0.3, axis='y' if plot_type == 'histogram' else 'both')


def create_plots(df, col_label, title_suffix, output_path, plot_type):
    """Create 2-panel (Pacer/Trotter) histogram or density plot."""
    figsize = (12, 14) if plot_type == 'histogram' else (14, 12)
    fig, axes = plt.subplots(2, 1, figsize=figsize, sharex=False)
    for ax, gait in zip(axes, ['Pacer', 'Trotter']):
        _plot_gait(ax, df[df['gait_type'] == gait], gait, col_label, title_suffix, plot_type)
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"Plot saved as: {output_path}")
    plt.close()


def print_statistics(df, title_suffix):
    """Print comprehensive statistics to stdout."""
    print("\n" + "=" * 60)
    print(f"STATISTICS: {title_suffix.upper()}")
    print("=" * 60)
    for gait in ['Pacer', 'Trotter']:
        gd = df[df['gait_type'] == gait]
        if gd.empty:
            continue
        print(f"\n--- {gait.upper()} (n={len(gd)}) ---")
        print(f"  Mean:   {gd['_metric'].mean():.2f}")
        print(f"  Median: {gd['_metric'].median():.2f}")
        print(f"  Std:    {gd['_metric'].std():.2f}")
        print(f"  Range:  [{gd['_metric'].min():.2f}, {gd['_metric'].max():.2f}]")
        print("  By book size:")
        for bs in BOOK_SIZE_ORDER:
            subset = gd[gd['book_size'] == bs]
            if len(subset) > 0:
                print(f"    {bs} (n={len(subset)}): "
                      f"Mean={subset['_metric'].mean():.2f}, "
                      f"Median={subset['_metric'].median():.2f}")

    corr = df['Percent_of_Consensus_ROH'].corr(df['F_ROH'])
    print(f"\nCorrelation (ROH_shared vs F_ROH): {corr:.3f}")


def main():
    parser = argparse.ArgumentParser(
        description='Create ROH histograms (ratio or shared %) by gait and book size.')
    parser.add_argument('--metric', choices=['ratio', 'shared'], required=True,
                        help='"ratio": ROH_shared/F_ROH; "shared": Percent_of_Consensus_ROH')
    parser.add_argument('roh_shared_file', help='Per-sample consensus ROH intersection file')
    parser.add_argument('froh_file', help='Per-sample F_ROH file (with gait column)')
    parser.add_argument('output_prefix', help='Output file prefix (will add .histogram.png, .density.png)')
    args = parser.parse_args()

    df = load_and_merge(args.roh_shared_file, args.froh_file)
    df, col_label, title_suffix = add_metric(df, args.metric)

    create_plots(df, col_label, title_suffix, args.output_prefix + '.histogram.png', 'histogram')
    create_plots(df, col_label, title_suffix, args.output_prefix + '.density.png', 'density')
    print_statistics(df, title_suffix)


if __name__ == '__main__':
    main()
