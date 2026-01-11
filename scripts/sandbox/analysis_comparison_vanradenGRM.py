import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import pearsonr
import allel
import sys

def compute_vanraden_grm(genotype_matrix):
    """
    Computes Standard GRM (VanRaden Method 1).
    Input: Genotype Matrix (0, 1, 2) of shape (N_ind, M_snps)
    """
    print("Computing Standard VanRaden GRM...")
    
    # 1. Calculate Allele Frequencies (p)
    # Mean of each column / 2
    p = np.mean(genotype_matrix, axis=0) / 2.0
    
    # 2. Center the Genotype Matrix (Z = M - 2p)
    Z = genotype_matrix - (2 * p)
    
    # 3. Calculate Scaling Factor 2*sum(p*(1-p))
    denominator = 2 * np.sum(p * (1 - p))
    
    # 4. Compute GRM: G = Z * Z' / denominator
    # Standard dot product of centered matrix
    G_std = np.dot(Z, Z.T) / denominator
    
    return G_std

# ==========================================
# MAIN COMPARISON ANALYSIS
# ==========================================
if len(sys.argv) != 4:
    print("Usage: python ROH_analysis.py <filtered.norm.phased.vcf.gz> <ROH_Relationship_Matrix> <output_dir>")
    sys.exit(1)

VCF_FILE = sys.argv[1]
ROH_Relationship_Matrix = sys.argv[2]
out_dir = sys.argv[3]

# 1. Load Data Again (We need raw genotypes 0/1/2 for Standard GRM)
print("Loading VCF for Standard GRM calculation...")
callset = allel.read_vcf(VCF_FILE, fields=['calldata/GT'])
gt = allel.GenotypeArray(callset['calldata/GT'])
# Sum alleles to get 0/1/2 dosage (N_snps x N_ind) -> Transpose to (N_ind x N_snps)
G_dosage = gt.to_n_alt().T 

# 2. Compute Standard GRM
G_std = compute_vanraden_grm(G_dosage)

# 3. Load your ROH Matrix
print("Loading ROH Relationship Matrix...")
G_roh = np.loadtxt(ROH_Relationship_Matrix)

# 4. Compare Diagonals (Inbreeding Estimates)
diag_std = np.diag(G_std) - 1.0  # VanRaden diagonal is 1+F
diag_roh = np.diag(G_roh)        # ROH diagonal is frequency of shared haplotypes

print(f"\nCorrelation between Inbreeding estimates: {pearsonr(diag_std, diag_roh)[0]:.4f}")

# 5. Compare Off-Diagonals (Relationships)
# We flatten the matrices and take the upper triangle to avoid duplicates
mask = np.triu_indices_from(G_std, k=1)
off_diag_std = G_std[mask]
off_diag_roh = G_roh[mask]

corr_rel = pearsonr(off_diag_std, off_diag_roh)[0]
print(f"Correlation between Relationship estimates: {corr_rel:.4f}")

# ==========================================
# VISUALIZATION
# ==========================================
plt.figure(figsize=(12, 5))

# Plot 1: Scatter of Off-Diagonal Elements
plt.subplot(1, 2, 1)
plt.scatter(off_diag_std, off_diag_roh, alpha=0.5, s=1)
plt.xlabel("Standard GRM Relationship")
plt.ylabel("ROH-Based Relationship")
plt.title(f"Relationship Comparison (r={corr_rel:.3f})")
plt.plot([0, max(off_diag_std)], [0, max(off_diag_std)], 'r--') # 1:1 line

# Plot 2: Distribution of Relationships
plt.subplot(1, 2, 2)
sns.kdeplot(off_diag_std, label='Standard GRM', fill=True)
sns.kdeplot(off_diag_roh, label='ROH Matrix', fill=True)
plt.xlabel("Relationship Coefficient")
plt.title("Distribution of Genomic Relationships")
plt.legend()

plt.tight_layout()
plt.savefig(f"{out_dir}/ROHRM_vs_VanradenGRM_Comparison.png")
print("\nPlot saved to 'ROHRM_vs_VanradenGRM_Comparison.png'. Check this image!")
