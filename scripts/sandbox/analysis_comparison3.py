import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import pearsonr, linregress
import sys
import os

class RobustMatrixComparator:
    def __init__(self, roh_prefix, std_prefix):
        self.roh_mat_file = f"{roh_prefix}.txt"
        self.roh_id_file  = f"{roh_prefix}.id"
        self.std_mat_file = f"{std_prefix}.rel"
        self.std_id_file  = f"{std_prefix}.rel.id"

        self.G_roh = None
        self.G_std = None
        self.common_ids = []
        self.pheno_map = {} 

    def load_ids(self, filename):
        print(f"  Reading IDs from {filename}...")
        try:
            with open(filename, 'r') as f:
                header_line = f.readline()
            
            if header_line.startswith('#'):
                df = pd.read_csv(filename, sep=r'\s+')
                df.columns = [c.replace('#', '') for c in df.columns]
                col = 'IID' if 'IID' in df.columns else df.columns[1]
                return df[col].astype(str).values
            else:
                df = pd.read_csv(filename, sep=r'\s+', header=None)
                col_idx = 1 if df.shape[1] >= 2 else 0
                return df.iloc[:, col_idx].astype(str).values

        except Exception as e:
            print(f"[ERROR] Failed to read ID file {filename}: {e}")
            sys.exit(1)

    def load_and_align(self):
        print("--- Loading and Aligning Matrices ---")
        roh_ids = self.load_ids(self.roh_id_file)
        std_ids = self.load_ids(self.std_id_file)
        
        common_samples = np.intersect1d(roh_ids, std_ids)
        print(f"  Overlapping Samples: {len(common_samples)}")
        
        if len(common_samples) == 0:
            print("[ERROR] No common IIDs found!")
            sys.exit(1)

        print("  Loading Matrices...")
        raw_roh = np.loadtxt(self.roh_mat_file)
        raw_std = np.loadtxt(self.std_mat_file)
        
        # Mapping
        roh_idx_map = {iid: i for i, iid in enumerate(roh_ids)}
        std_idx_map = {iid: i for i, iid in enumerate(std_ids)}
        
        roh_indices = [roh_idx_map[iid] for iid in common_samples]
        std_indices = [std_idx_map[iid] for iid in common_samples]
        
        print("  Re-ordering...")
        self.G_roh = raw_roh[np.ix_(roh_indices, roh_indices)]
        self.G_std = raw_std[np.ix_(std_indices, std_indices)]
        self.common_ids = common_samples

    def load_phenotypes(self, phenotype_file):
        if not os.path.exists(phenotype_file):
            return
        try:
            df = pd.read_csv(phenotype_file, sep=r'\s+', header=None)
            # Try to robustly detect ID vs Phenotype columns
            if df.shape[1] >= 3:
                self.pheno_map = dict(zip(df.iloc[:, 1].astype(str), df.iloc[:, 2]))
            elif df.shape[1] == 2:
                self.pheno_map = dict(zip(df.iloc[:, 0].astype(str), df.iloc[:, 1]))
        except Exception as e:
            print(f"Warning: Could not load phenotypes: {e}")

    def save_inbreeding_data(self, output_file="Inbreeding_Comparison.csv"):
        print(f"\n--- Saving Inbreeding Data to {output_file} ---")
        diag_std = np.diag(self.G_std) - 1.0
        diag_roh = np.diag(self.G_roh) - 1.0 
        
        data = []
        for i, iid in enumerate(self.common_ids):
            data.append({
                'IID': iid,
                'F_Standard': diag_std[i],
                'F_ROH': diag_roh[i],
                'Phenotype': self.pheno_map.get(str(iid), 'Unknown')
            })
        pd.DataFrame(data).to_csv(output_file, index=False)
        print("Done.")

    def compare_and_plot(self):
        print(f"\n--- Comparing Matrices ---")
        
        # 1. Scaling
        diag_std = np.diag(self.G_std) - 1.0
        diag_roh = np.diag(self.G_roh) - 1.0
        
        # Global Correlations
        corr_diag, _ = pearsonr(diag_std, diag_roh)
        
        mask = np.triu_indices_from(self.G_std, k=1)
        off_diag_std = self.G_std[mask]
        off_diag_roh = self.G_roh[mask]
        
        corr_rel, _ = pearsonr(off_diag_std, off_diag_roh)
        
        print(f"Global Correlation (Inbreeding): r = {corr_diag:.4f}")
        print(f"Global Correlation (Relationships): r = {corr_rel:.4f}")

        # Generate Plots
        self._generate_plots(diag_std, diag_roh, off_diag_std, off_diag_roh, corr_diag, corr_rel)

    def _generate_plots(self, d_std, d_roh, od_std, od_roh, r_diag, r_rel):
        plt.figure(figsize=(16, 12))
        sns.set_style("whitegrid")

        # --- PLOT 1: Inbreeding Scatter (Top Left) ---
        plt.subplot(2, 2, 1)
        phenos = [self.pheno_map.get(iid, 'Unknown') for iid in self.common_ids]
        df_diag = pd.DataFrame({'Std': d_std, 'ROH': d_roh, 'Pop': phenos})
        
        # Plot points
        sns.scatterplot(data=df_diag, x='Std', y='ROH', hue='Pop', style='Pop', alpha=0.7, s=30)
        
        # Calculate & Print Subpop Stats
        print("\n--- Subpopulation Inbreeding Correlations ---")
        unique_pops = [p for p in df_diag['Pop'].unique() if p != 'Unknown']
        
        annot_text = f"Global r = {r_diag:.3f}\n"
        colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
        
        for i, pop in enumerate(unique_pops):
            sub = df_diag[df_diag['Pop'] == pop]
            if len(sub) > 5:
                r_sub, _ = pearsonr(sub['Std'], sub['ROH'])
                print(f"  > {pop}: r = {r_sub:.4f}")
                annot_text += f"{pop} r = {r_sub:.3f}\n"
                
                # Add trend line
                m, b = np.polyfit(sub['Std'], sub['ROH'], 1)
                plt.plot(sub['Std'], m*sub['Std'] + b, linestyle='--', linewidth=2, color=colors[i])

        # Add text box to plot
        plt.text(0.05, 0.95, annot_text, transform=plt.gca().transAxes, 
                 fontsize=10, verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.9))

        plt.xlabel("Standard Inbreeding ($F_{GRM}$)")
        plt.ylabel("ROH Inbreeding ($F_{ROH}$)")
        plt.title(f"Inbreeding Comparison")

        # --- PLOT 2: Relationship Scatter (Top Right) ---
        plt.subplot(2, 2, 2)
        # Subsample for speed/clarity
        if len(od_std) > 50000:
            idx = np.random.choice(len(od_std), 50000, replace=False)
            x_vals, y_vals = od_std[idx], od_roh[idx]
            plt.scatter(x_vals, y_vals, alpha=0.2, s=2, c='royalblue', label='Pairs (Subsampled)')
        else:
            x_vals, y_vals = od_std, od_roh
            plt.scatter(x_vals, y_vals, alpha=0.2, s=2, c='royalblue', label='Pairs')

        # 1:1 Line (Red Dashed)
        plt.plot([min(od_std), max(od_std)], [min(od_std), max(od_std)], 'r--', linewidth=1.5, label="1:1 Line")
        
        # Regression Line (Green Solid) - Shows the actual shift
        slope, intercept, _, _, _ = linregress(x_vals, y_vals)
        plt.plot(x_vals, slope*x_vals + intercept, 'g-', linewidth=1.5, alpha=0.8, label=f"Fit (y={slope:.2f}x + {intercept:.2f})")

        plt.xlabel("Standard Kinship")
        plt.ylabel("ROH Kinship")
        plt.legend()
        plt.title(f"Relationship Comparison (Global r={r_rel:.3f})")

        # --- PLOT 3: Distributions (Bottom Left) ---
        plt.subplot(2, 2, 3)
        sns.kdeplot(od_std, fill=True, label='Standard GRM', color='blue', alpha=0.2)
        sns.kdeplot(od_roh, fill=True, label='ROH GRM', color='orange', alpha=0.2)
        plt.xlabel("Relationship Coefficient")
        plt.title("Distribution of Kinship Values")
        plt.legend()

        # --- PLOT 4: Heatmap (Bottom Right) ---
        plt.subplot(2, 2, 4)
        
        if unique_pops:
            # Sort by Phenotype to reveal structure
            sorted_indices = np.argsort(phenos)
            
            # Select subset for heatmap if too large (e.g., top 2000)
            limit = 2000
            if len(sorted_indices) > limit:
                print(f"Heatmap: Downsampling to {limit} samples for visibility.")
                sorted_indices = sorted_indices[:limit]

            sorted_G_roh = self.G_roh[np.ix_(sorted_indices, sorted_indices)]
            sorted_G_std = self.G_std[np.ix_(sorted_indices, sorted_indices)]
            
            # Calculate Difference
            diff_matrix = sorted_G_roh - sorted_G_std
            
            # Determine Center for Colorbar (Mean Difference)
            # This fixes the "All Red" issue by making the average difference "White"
            center_val = np.median(diff_matrix)
            
            sns.heatmap(diff_matrix, cmap="RdBu_r", center=center_val, 
                        cbar_kws={'label': f'Difference (Centered at {center_val:.2f})'},
                        xticklabels=False, yticklabels=False)
            
            plt.title(f"Difference Matrix ($G_{{ROH}} - G_{{STD}}$)\nSorted by Pop (Red = ROH relatively higher)")
        else:
            plt.text(0.5, 0.5, "Need Phenotypes for Sorted Heatmap", ha='center')
            plt.axis('off')

        plt.tight_layout()
        plt.savefig("Robust_Matrix_Comparison_Enhanced.png", dpi=300)
        print("\nSaved high-res plot to 'Robust_Matrix_Comparison_Enhanced.png'")

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python analysis_comparison.py <ROH_Prefix> <Std_Prefix> [Phenotypes_File]")
        sys.exit(1)

    ROH_PREFIX = sys.argv[1]
    STD_PREFIX = sys.argv[2]
    PHENO_FILE = sys.argv[3] if len(sys.argv) > 3 else "phenotypes.txt"

    comp = RobustMatrixComparator(ROH_PREFIX, STD_PREFIX)
    comp.load_phenotypes(PHENO_FILE) 
    comp.load_and_align()
    
    comp.compare_and_plot()
    comp.save_inbreeding_data()

