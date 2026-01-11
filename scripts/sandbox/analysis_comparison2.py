import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import pearsonr
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
        self.pheno_map = {} # Store phenotypes for labeling

    def load_ids(self, filename):
        """ Robustly load IIDs using raw strings for regex. """
        print(f"  Reading IDs from {filename}...")
        try:
            with open(filename, 'r') as f:
                header_line = f.readline()
            
            # Use raw string r'\s+' for separator to avoid warnings
            if header_line.startswith('#'):
                df = pd.read_csv(filename, sep=r'\s+')
                df.columns = [c.replace('#', '') for c in df.columns]
                if 'IID' in df.columns:
                    return df['IID'].astype(str).values
                else:
                    return df.iloc[:, 1].astype(str).values
            else:
                df = pd.read_csv(filename, sep=r'\s+', header=None)
                if df.shape[1] >= 2:
                    return df.iloc[:, 1].astype(str).values
                else:
                    return df.iloc[:, 0].astype(str).values

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
        """ Loads phenotypes for plotting and sorting """
        if not os.path.exists(phenotype_file):
            return
        try:
            # Try 3 columns (FID IID Pheno) or 2 (IID Pheno)
            df = pd.read_csv(phenotype_file, sep=r'\s+', header=None)
            if df.shape[1] == 3:
                self.pheno_map = dict(zip(df.iloc[:, 1].astype(str), df.iloc[:, 2]))
            else:
                self.pheno_map = dict(zip(df.iloc[:, 0].astype(str), df.iloc[:, 1]))
        except Exception as e:
            print(f"Warning: Could not load phenotypes: {e}")

    def save_inbreeding_data(self, output_file="Inbreeding_Comparison.csv"):
        print(f"\n--- Saving Inbreeding Data to {output_file} ---")
        
        # Subtract 1.0 from BOTH to get F on 0-1 scale
        diag_std = np.diag(self.G_std) - 1.0
        diag_roh = np.diag(self.G_roh) - 1.0  # FIX: Subtracted 1 here too
        
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
        
        # 1. Correct Scaling (Subtract 1 from both diagonals)
        diag_std = np.diag(self.G_std) - 1.0
        diag_roh = np.diag(self.G_roh) - 1.0
        
        corr_diag, _ = pearsonr(diag_std, diag_roh)
        print(f"Correlation (Inbreeding): r = {corr_diag:.4f}")

        # 2. Off-Diagonals
        mask = np.triu_indices_from(self.G_std, k=1)
        off_diag_std = self.G_std[mask]
        off_diag_roh = self.G_roh[mask]
        
        corr_rel, _ = pearsonr(off_diag_std, off_diag_roh)
        print(f"Correlation (Relationships): r = {corr_rel:.4f}")

        self._generate_plots(diag_std, diag_roh, off_diag_std, off_diag_roh, corr_diag, corr_rel)

    def _generate_plots(self, d_std, d_roh, od_std, od_roh, r_diag, r_rel):
        plt.figure(figsize=(14, 12))

        # Plot A: Inbreeding Scatter (Colored by Pop if avail)
        plt.subplot(2, 2, 1)
        phenos = [self.pheno_map.get(iid, 'Unknown') for iid in self.common_ids]
        
        # Create temp DF for seaborn
        df_diag = pd.DataFrame({'Std': d_std, 'ROH': d_roh, 'Pop': phenos})
        sns.scatterplot(data=df_diag, x='Std', y='ROH', hue='Pop', alpha=0.6, s=15)
        
        # Add trend lines per pop
        unique_pops = [p for p in df_diag['Pop'].unique() if p != 'Unknown']
        for pop in unique_pops:
            sub = df_diag[df_diag['Pop'] == pop]
            if len(sub) > 2:
                m, b = np.polyfit(sub['Std'], sub['ROH'], 1)
                plt.plot(sub['Std'], m*sub['Std'] + b, linestyle='--', linewidth=1)

        plt.xlabel("Standard Inbreeding ($F_{GRM}$)")
        plt.ylabel("ROH Inbreeding ($F_{ROH}$)")
        plt.title(f"Inbreeding Comparison (r={r_diag:.3f})")

        # Plot B: Relationship Scatter
        plt.subplot(2, 2, 2)
        # Subsample for speed
        if len(od_std) > 50000:
            idx = np.random.choice(len(od_std), 50000, replace=False)
            plt.scatter(od_std[idx], od_roh[idx], alpha=0.2, s=1, c='blue')
        else:
            plt.scatter(od_std, od_roh, alpha=0.2, s=1, c='blue')
        plt.plot([min(od_std), max(od_std)], [min(od_std), max(od_std)], 'r--')
        plt.xlabel("Standard Kinship")
        plt.ylabel("ROH Kinship")
        plt.title(f"Relationship Comparison (r={r_rel:.3f})")

        # Plot C: Density Distribution (The "Bimodal" check)
        plt.subplot(2, 2, 3)
        sns.kdeplot(od_std, fill=True, label='Standard GRM', color='blue', alpha=0.3)
        sns.kdeplot(od_roh, fill=True, label='ROH GRM', color='orange', alpha=0.3)
        plt.xlabel("Relationship Coefficient")
        plt.title("Distribution of Relationships")
        plt.legend()

        # Plot D: Difference Heatmap (New!)
        # Shows (G_ROH - G_Std)
        # We need to sort the matrix by Population to make sense of it
        plt.subplot(2, 2, 4)
        
        if unique_pops:
            # Sort indices by phenotype
            sorted_indices = np.argsort(phenos)
            sorted_G_roh = self.G_roh[np.ix_(sorted_indices, sorted_indices)]
            sorted_G_std = self.G_std[np.ix_(sorted_indices, sorted_indices)]
            diff_matrix = sorted_G_roh - sorted_G_std
            
            # Downsample if matrix is huge (e.g., > 1000x1000) for heatmap performance
            if diff_matrix.shape[0] > 1000:
                 # Just take top 1000 for visualization or block average (skipping for simplicity)
                 diff_matrix = diff_matrix[:1000, :1000]
                 title_suffix = "(Top 1000 sorted samples)"
            else:
                 title_suffix = "(All sorted samples)"

            sns.heatmap(diff_matrix, cmap="RdBu_r", center=0, cbar_kws={'label': 'ROH - Standard'})
            plt.title(f"Difference Matrix {title_suffix}\nRed = ROH Higher, Blue = Std Higher")
            plt.axis('off')
        else:
            plt.text(0.5, 0.5, "Need Phenotypes for Sorted Heatmap", ha='center')

        plt.tight_layout()
        plt.savefig("Robust_Matrix_Comparison.png")
        print("Saved plots to 'Robust_Matrix_Comparison.png'")

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python analysis_comparison.py <ROH_Prefix> <Std_Prefix> [Phenotypes_File]")
        sys.exit(1)

    ROH_PREFIX = sys.argv[1]
    STD_PREFIX = sys.argv[2]
    PHENO_FILE = sys.argv[3] if len(sys.argv) > 3 else "phenotypes.txt"

    comp = RobustMatrixComparator(ROH_PREFIX, STD_PREFIX)
    comp.load_phenotypes(PHENO_FILE) # Load first to help with alignment sorting if needed later
    comp.load_and_align()
    
    comp.compare_and_plot()
    comp.save_inbreeding_data()

