import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import sys
from scipy import stats

def create_roh_ratio_histograms(file1_path, file2_path, output_path="roh_ratio_histograms.png"):
    """
    Create histograms of Percent_of_Consensus_ROH / F_ROH ratio
    for Pacer and Trotter gaits, faceted by book size.
    """
    
    # Read the data files
    try:
        df1 = pd.read_csv(file1_path, sep=r'\s+')
        df2 = pd.read_csv(file2_path, sep=r'\s+')
    except FileNotFoundError as e:
        print(f"Error: {e}")
        sys.exit(1)
    
    # Merge the dataframes on IID
    df_merged = pd.merge(df1, df2, on='IID', how='inner')
    
    # Extract gait type and book size from the 'gait' column
    df_merged['gait_type'] = df_merged['gait'].str.split('_').str[0]
    df_merged['book_size'] = df_merged['gait'].str.split('_').str[1]
    
    # Calculate the ratio: Percent_of_Consensus_ROH / F_ROH
    # Adding a small epsilon to avoid division by zero
    epsilon = 1e-10
    df_merged['ROH_F_ratio'] = (df_merged['Percent_of_Consensus_ROH'] / 100) / (df_merged['F_ROH'] + epsilon)
    
    # Order for book size categories
    book_size_order = ['HIGH', 'MEDIUM', 'LOW']
    df_merged['book_size'] = pd.Categorical(df_merged['book_size'], 
                                           categories=book_size_order, 
                                           ordered=True)
    
    # Create subplots (2 rows, 1 column)
    fig, axes = plt.subplots(2, 1, figsize=(12, 14), sharex=False)
    
    # Color palette for book sizes
    colors = {'HIGH': '#E74C3C', 'MEDIUM': '#F39C12', 'LOW': '#2ECC71'}
    
    # Create histogram for Pacer
    ax1 = axes[0]
    pacer_data = df_merged[df_merged['gait_type'] == 'Pacer']
    
    # Prepare data for each book size
    hist_data = []
    labels = []
    for book_size in book_size_order:
        subset = pacer_data[pacer_data['book_size'] == book_size]
        if len(subset) > 0:
            hist_data.append(subset['ROH_F_ratio'])
            labels.append(f'{book_size} (n={len(subset)})')
    
    # Plot histogram
    if hist_data:
        n_bins = min(20, int(np.sqrt(len(pacer_data))))
        ax1.hist(hist_data, bins=n_bins, stacked=False, 
                color=[colors[bs] for bs in book_size_order[:len(hist_data)]],
                alpha=0.7, edgecolor='black', linewidth=0.8,
                label=labels)
        
        # Add vertical line for mean of each group
        for i, (book_size, color) in enumerate(zip(book_size_order[:len(hist_data)], [colors[bs] for bs in book_size_order[:len(hist_data)]])):
            subset = pacer_data[pacer_data['book_size'] == book_size]
            if len(subset) > 0:
                mean_val = subset['ROH_F_ratio'].mean()
                ax1.axvline(mean_val, color=color, linestyle='--', linewidth=2, 
                          alpha=0.8, label=f'{book_size} Mean: {mean_val:.1f}')
    
    ax1.set_title('Pacer - Distribution of ROH_shared / F_ROH Ratio', 
                  fontsize=14, fontweight='bold', pad=15)
    ax1.set_xlabel('ROH_shared / F_ROH', fontsize=12)
    ax1.set_ylabel('Frequency', fontsize=12)
    ax1.grid(True, alpha=0.3, axis='y')
    ax1.legend(fontsize=10, title='Book Size', title_fontsize=11)
    
    # Add summary statistics text for Pacer
    if not pacer_data.empty:
        stats_text = f'Overall Mean: {pacer_data["ROH_F_ratio"].mean():.1f}\n'
        stats_text += f'Overall Std: {pacer_data["ROH_F_ratio"].std():.1f}\n'
        stats_text += f'Overall Median: {pacer_data["ROH_F_ratio"].median():.1f}\n'
        stats_text += f'Min: {pacer_data["ROH_F_ratio"].min():.1f}, Max: {pacer_data["ROH_F_ratio"].max():.1f}'
        
        ax1.text(0.02, 0.98, stats_text, transform=ax1.transAxes,
                fontsize=9, verticalalignment='top',
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
    
    # Create histogram for Trotter
    ax2 = axes[1]
    trotter_data = df_merged[df_merged['gait_type'] == 'Trotter']
    
    # Prepare data for each book size
    hist_data = []
    labels = []
    for book_size in book_size_order:
        subset = trotter_data[trotter_data['book_size'] == book_size]
        if len(subset) > 0:
            hist_data.append(subset['ROH_F_ratio'])
            labels.append(f'{book_size} (n={len(subset)})')
    
    # Plot histogram
    if hist_data:
        n_bins = min(20, int(np.sqrt(len(trotter_data))))
        ax2.hist(hist_data, bins=n_bins, stacked=False,
                color=[colors[bs] for bs in book_size_order[:len(hist_data)]],
                alpha=0.7, edgecolor='black', linewidth=0.8,
                label=labels)
        
        # Add vertical line for mean of each group
        for i, (book_size, color) in enumerate(zip(book_size_order[:len(hist_data)], [colors[bs] for bs in book_size_order[:len(hist_data)]])):
            subset = trotter_data[trotter_data['book_size'] == book_size]
            if len(subset) > 0:
                mean_val = subset['ROH_F_ratio'].mean()
                ax2.axvline(mean_val, color=color, linestyle='--', linewidth=2,
                          alpha=0.8, label=f'{book_size} Mean: {mean_val:.1f}')
    
    ax2.set_title('Trotter - Distribution of ROH_shared / F_ROH Ratio', 
                  fontsize=14, fontweight='bold', pad=15)
    ax2.set_xlabel('ROH_shared / F_ROH', fontsize=12)
    ax2.set_ylabel('Frequency', fontsize=12)
    ax2.grid(True, alpha=0.3, axis='y')
    ax2.legend(fontsize=10, title='Book Size', title_fontsize=11)
    
    # Add summary statistics text for Trotter
    if not trotter_data.empty:
        stats_text = f'Overall Mean: {trotter_data["ROH_F_ratio"].mean():.1f}\n'
        stats_text += f'Overall Std: {trotter_data["ROH_F_ratio"].std():.1f}\n'
        stats_text += f'Overall Median: {trotter_data["ROH_F_ratio"].median():.1f}\n'
        stats_text += f'Min: {trotter_data["ROH_F_ratio"].min():.1f}, Max: {trotter_data["ROH_F_ratio"].max():.1f}'
        
        ax2.text(0.02, 0.98, stats_text, transform=ax2.transAxes,
                fontsize=9, verticalalignment='top',
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
    
    # Adjust layout
    plt.tight_layout()
    
    # Save the figure
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"Histograms saved as: {output_path}")
    
    # Show comprehensive statistics
    print("\n" + "="*60)
    print("COMPREHENSIVE STATISTICS: ROH_shared / F_ROH RATIO")
    print("="*60)
    
    for gait in ['Pacer', 'Trotter']:
        gait_data = df_merged[df_merged['gait_type'] == gait]
        if not gait_data.empty:
            print(f"\n--- {gait.upper()} ---")
            print(f"Total samples: {len(gait_data)}")
            print(f"Overall ratio statistics:")
            print(f"  Mean: {gait_data['ROH_F_ratio'].mean():.2f}")
            print(f"  Median: {gait_data['ROH_F_ratio'].median():.2f}")
            print(f"  Std: {gait_data['ROH_F_ratio'].std():.2f}")
            print(f"  Range: [{gait_data['ROH_F_ratio'].min():.2f}, {gait_data['ROH_F_ratio'].max():.2f}]")
            
            print("\n  By book size:")
            for book_size in book_size_order:
                subset = gait_data[gait_data['book_size'] == book_size]
                if len(subset) > 0:
                    print(f"    {book_size} (n={len(subset)}): "
                          f"Mean={subset['ROH_F_ratio'].mean():.2f}, "
                          f"Median={subset['ROH_F_ratio'].median():.2f}")
    
    return df_merged

# Alternative version with density plots and kernel density estimation
def create_roh_ratio_density_plots(file1_path, file2_path, output_path="roh_ratio_density.png"):
    """
    Create density plots of ROH/F_ROH ratio with kernel density estimation.
    """
    
    # Read and process data
    df1 = pd.read_csv(file1_path, sep=r'\s+')
    df2 = pd.read_csv(file2_path, sep=r'\s+')
    df_merged = pd.merge(df1, df2, on='IID', how='inner')
    
    # Extract columns and calculate ratio
    df_merged['gait_type'] = df_merged['gait'].str.split('_').str[0]
    df_merged['book_size'] = df_merged['gait'].str.split('_').str[1]
    epsilon = 1e-10
    df_merged['ROH_F_ratio'] = (df_merged['Percent_of_Consensus_ROH'] / 100) / (df_merged['F_ROH'] + epsilon)
    
    # Order categories
    book_size_order = ['HIGH', 'MEDIUM', 'LOW']
    df_merged['book_size'] = pd.Categorical(df_merged['book_size'], 
                                           categories=book_size_order, 
                                           ordered=True)
    
    # Create figure with subplots
    fig, axes = plt.subplots(2, 1, figsize=(14, 12))
    
    # Color palette
    colors = {'HIGH': '#E74C3C', 'MEDIUM': '#F39C12', 'LOW': '#2ECC71'}
    
    # Plot density for Pacer
    ax1 = axes[0]
    pacer_data = df_merged[df_merged['gait_type'] == 'Pacer']
    
    for book_size in book_size_order:
        subset = pacer_data[pacer_data['book_size'] == book_size]
        if len(subset) > 0:
            # Plot histogram
            ax1.hist(subset['ROH_F_ratio'], bins=15, density=True, alpha=0.3,
                    color=colors[book_size], edgecolor='black')
            
            # Plot kernel density estimate
            from scipy.stats import gaussian_kde
            kde = gaussian_kde(subset['ROH_F_ratio'])
            x_range = np.linspace(subset['ROH_F_ratio'].min() * 0.9, 
                                 subset['ROH_F_ratio'].max() * 1.1, 200)
            ax1.plot(x_range, kde(x_range), color=colors[book_size], 
                    linewidth=2, label=f'{book_size} (n={len(subset)})')
    
    ax1.set_title('Pacer - Density Distribution of ROH_shared / F_ROH Ratio',
                  fontsize=14, fontweight='bold')
    ax1.set_xlabel('ROH_shared / F_ROH', fontsize=12)
    ax1.set_ylabel('Density', fontsize=12)
    ax1.legend(fontsize=10)
    ax1.grid(True, alpha=0.3)
    
    # Plot density for Trotter
    ax2 = axes[1]
    trotter_data = df_merged[df_merged['gait_type'] == 'Trotter']
    
    for book_size in book_size_order:
        subset = trotter_data[trotter_data['book_size'] == book_size]
        if len(subset) > 0:
            # Plot histogram
            ax2.hist(subset['ROH_F_ratio'], bins=15, density=True, alpha=0.3,
                    color=colors[book_size], edgecolor='black')
            
            # Plot kernel density estimate
            kde = gaussian_kde(subset['ROH_F_ratio'])
            x_range = np.linspace(subset['ROH_F_ratio'].min() * 0.9,
                                 subset['ROH_F_ratio'].max() * 1.1, 200)
            ax2.plot(x_range, kde(x_range), color=colors[book_size],
                    linewidth=2, label=f'{book_size} (n={len(subset)})')
    
    ax2.set_title('Trotter - Density Distribution of ROH_shared / F_ROH Ratio',
                  fontsize=14, fontweight='bold')
    ax2.set_xlabel('ROH_shared / F_ROH', fontsize=12)
    ax2.set_ylabel('Density', fontsize=12)
    ax2.legend(fontsize=10)
    ax2.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"Density plots saved as: {output_path}")
    
    return df_merged

def main():
    if len(sys.argv) != 4:
        print("Usage: python script.py <Roh_shared_file> <froh_file> <output_prefix>")
        sys.exit(1)
    
    file1 = sys.argv[1]
    file2 = sys.argv[2]
    outHisto = sys.argv[3]+".histogram.png"
    outDensity = sys.argv[3]+".density.png"
    
    # Create histogram plots
    df = create_roh_ratio_histograms(file1, file2, outHisto)
    
    # Optional: Create density plots as well
    create_roh_ratio_density_plots(file1, file2, outDensity)
    
    # Show additional insights
    print("\n" + "="*60)
    print("BIOLOGICAL INTERPRETATION INSIGHTS:")
    print("="*60)
    print("""
The ROH_shared / F_ROH ratio represents:
- How much shared ROH (Percent_of_Consensus_ROH) is present per unit of individual inbreeding (F_ROH)
- Higher ratio: More shared segments relative to individual inbreeding
- Lower ratio: More individual-specific inbreeding relative to shared segments

Interpretation by gait and book size:
1. Higher ratios may indicate more recent shared ancestry within groups
2. Lower ratios may suggest more ancient or individual-specific inbreeding
3. Differences between Pacer/Trotter may reflect breeding strategies
4. Variation by book size could indicate population structure effects
""")
    
    # Calculate and display correlation
    correlation = df['Percent_of_Consensus_ROH'].corr(df['F_ROH'])
    print(f"\nCorrelation between ROH_shared and F_ROH: {correlation:.3f}")
    
    if correlation > 0.7:
        print("Strong positive correlation: Shared and individual ROH tend to increase together")
    elif correlation > 0.3:
        print("Moderate positive correlation")
    elif correlation > -0.3:
        print("Weak correlation")
    else:
        print("Negative correlation: As individual inbreeding increases, shared segments decrease")

if __name__ == "__main__":
    main()
