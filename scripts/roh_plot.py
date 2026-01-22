import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from scipy import stats
import sys

def create_enhanced_roh_scatter_plots(file1_path, file2_path, output_path):
    """
    Create enhanced scatter diagrams with regression lines and statistics.
    """
    
    # Read and merge data
    df1 = pd.read_csv(file1_path, sep=r'\s+')
    df2 = pd.read_csv(file2_path, sep=r'\s+')
    df_merged = pd.merge(df1, df2, on='IID', how='inner')
    
    # Extract gait type and book size
    df_merged['gait_type'] = df_merged['gait'].str.split('_').str[0]
    df_merged['book_size'] = df_merged['gait'].str.split('_').str[1]
    
    # Order categories
    book_size_order = ['HIGH', 'MEDIUM', 'LOW']
    df_merged['book_size'] = pd.Categorical(df_merged['book_size'], 
                                           categories=book_size_order, 
                                           ordered=True)
    
    # Create figure with custom style
    plt.style.use('seaborn-v0_8-darkgrid')
    fig = plt.figure(figsize=(12, 14))
    
    # Create two subplots
    ax1 = plt.subplot(2, 1, 1)  # Pacer
    ax2 = plt.subplot(2, 1, 2)  # Trotter
    
    # Color palette
    palette = {'HIGH': '#D32F2F', 'MEDIUM': '#FF9800', 'LOW': '#388E3C'}
    
    # Plot Pacer data
    pacer_data = df_merged[df_merged['gait_type'] == 'Pacer']
    for book_size in book_size_order:
        subset = pacer_data[pacer_data['book_size'] == book_size]
        if not subset.empty:
            ax1.scatter(subset['F_ROH'], subset['Percent_of_Consensus_ROH'],
                       color=palette[book_size], s=120, alpha=0.7,
                       edgecolor='black', linewidth=0.8,
                       label=f'Book Size: {book_size}', zorder=5)
    
    # Add regression line for Pacer
    if len(pacer_data) > 1:
        x_pacer = pacer_data['F_ROH']
        y_pacer = pacer_data['Percent_of_Consensus_ROH']
        slope, intercept, r_value, p_value, std_err = stats.linregress(x_pacer, y_pacer)
        x_fit = np.linspace(x_pacer.min(), x_pacer.max(), 100)
        y_fit = slope * x_fit + intercept
        ax1.plot(x_fit, y_fit, 'k--', linewidth=2, alpha=0.7,
                label=f'R² = {r_value**2:.3f}', zorder=4)
    
    ax1.set_title('Pacer - Relationship between F_ROH and ROH Shared', 
                  fontsize=16, fontweight='bold', pad=15)
    ax1.set_ylabel('ROH_shared (%)', fontsize=14, labelpad=10)
    ax1.tick_params(axis='both', labelsize=12)
    ax1.legend(fontsize=11, title_fontsize=12, loc='best')
    ax1.grid(True, alpha=0.4)
    
    # Add text box with Pacer statistics
    if len(pacer_data) > 0:
        stats_text = f'N = {len(pacer_data)}\nMean F_ROH = {pacer_data["F_ROH"].mean():.4f}\nMean ROH_sh = {pacer_data["Percent_of_Consensus_ROH"].mean():.2f}%'
        ax1.text(0.02, 0.98, stats_text, transform=ax1.transAxes,
                fontsize=10, verticalalignment='top',
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
    
    # Plot Trotter data
    trotter_data = df_merged[df_merged['gait_type'] == 'Trotter']
    for book_size in book_size_order:
        subset = trotter_data[trotter_data['book_size'] == book_size]
        if not subset.empty:
            ax2.scatter(subset['F_ROH'], subset['Percent_of_Consensus_ROH'],
                       color=palette[book_size], s=120, alpha=0.7,
                       edgecolor='black', linewidth=0.8,
                       label=f'Book Size: {book_size}', zorder=5)
    
    # Add regression line for Trotter
    if len(trotter_data) > 1:
        x_trotter = trotter_data['F_ROH']
        y_trotter = trotter_data['Percent_of_Consensus_ROH']
        slope, intercept, r_value, p_value, std_err = stats.linregress(x_trotter, y_trotter)
        x_fit = np.linspace(x_trotter.min(), x_trotter.max(), 100)
        y_fit = slope * x_fit + intercept
        ax2.plot(x_fit, y_fit, 'k--', linewidth=2, alpha=0.7,
                label=f'R² = {r_value**2:.3f}', zorder=4)
    
    ax2.set_title('Trotter - Relationship between F_ROH and ROH Shared', 
                  fontsize=16, fontweight='bold', pad=15)
    ax2.set_xlabel('F_ROH', fontsize=14, labelpad=10)
    ax2.set_ylabel('ROH_shared (%)', fontsize=14, labelpad=10)
    ax2.tick_params(axis='both', labelsize=12)
    ax2.legend(fontsize=11, title_fontsize=12, loc='best')
    ax2.grid(True, alpha=0.4)
    
    # Add text box with Trotter statistics
    if len(trotter_data) > 0:
        stats_text = f'N = {len(trotter_data)}\nMean F_ROH = {trotter_data["F_ROH"].mean():.4f}\nMean ROH_sh = {trotter_data["Percent_of_Consensus_ROH"].mean():.2f}%'
        ax2.text(0.02, 0.98, stats_text, transform=ax2.transAxes,
                fontsize=10, verticalalignment='top',
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
    
    # Adjust layout and save
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"Enhanced plot saved as: {output_path}")
    
    return df_merged

# Usage
if __name__ == "__main__":
    if len(sys.argv) == 4:
        create_enhanced_roh_scatter_plots(sys.argv[1], sys.argv[2], sys.argv[3])
    else:
        print("Usage: python enhanced_script.py <Roh_shared_file> <froh_file> <output_file>")
