import allel
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys
from tqdm import tqdm

class ROHAnalysis:
    def __init__(self, vcf_path, phenotype_path=None, roh_mb_cutoff=1.0):
        self.vcf_path = vcf_path
        self.pheno_path = phenotype_path
        self.roh_cutoff = roh_mb_cutoff
        
        # Data
        self.map_df = None       # SNP Map
        self.hap1 = None         # Paternal
        self.hap2 = None         # Maternal
        self.sample_ids = None   # IDs from VCF
        
        # Phenotypes
        self.subpops = {}        # Dict: {'Whole': [indices], 'Case': [indices], ...}
        
        # Analysis Results
        self.roh_matrix = None   # Boolean Matrix (N_ind x N_snp): True if SNP is in ROH

    def load_data(self):
        print(f"--- Loading Data ---")
        # 1. Load VCF (Same as before)
        callset = allel.read_vcf(self.vcf_path, fields=['calldata/GT', 'variants/CHROM', 'variants/POS', 'samples'])
        
        # Map
        chroms = [str(c).replace('chr', '') for c in callset['variants/CHROM']]
        chroms = np.array([int(c) if c.isdigit() else 0 for c in chroms])
        self.map_df = pd.DataFrame({
            'chrom': chroms,
            'pos_mb': callset['variants/POS'] / 1e6,
            'pos_bp': callset['variants/POS']
        })
        
        # Haplotypes
        gt = allel.GenotypeArray(callset['calldata/GT'])
        self.hap1 = np.array(gt[:, :, 0].T, dtype=np.int8)
        self.hap2 = np.array(gt[:, :, 1].T, dtype=np.int8)
        self.sample_ids = callset['samples']
        
        print(f"Genotypes: {self.hap1.shape[0]} Inds x {self.hap1.shape[1]} SNPs")

        # 2. Load Phenotypes (If provided)
        if self.pheno_path:
            self._load_phenotypes()
        else:
            # Default: Everyone is one population
            self.subpops['Whole'] = np.arange(len(self.sample_ids))

    def _load_phenotypes(self):
        """
        Reads headerless tab-separated file: [unused, ID, Phenotype]
        """
        print(f"Loading Phenotypes: {self.pheno_path}")
        try:
            # Read file (assuming no header, columns 0, 1, 2)
            # Col 1 = ID, Col 2 = Binary Phenotype
            df = pd.read_csv(self.pheno_path, sep='\t', header=None, dtype=str)
            
            # map IDs to VCF indices
            vcf_id_map = {id_str: i for i, id_str in enumerate(self.sample_ids)}
            
            pop_indices = {} # Key: Phenotype Label, Value: List of VCF indices
            
            found_count = 0
            for _, row in df.iterrows():
                # Adjust column indices if your file format differs
                # Assuming: Col 0 (ignore), Col 1 (ID), Col 2 (Pheno)
                # If file is only 2 columns, change to row[0] and row[1]
                if len(row) >= 3:
                    sid = row[1]
                    pheno = row[2]
                else:
                    sid = row[0] # Fallback
                    pheno = row[1]

                if sid in vcf_id_map:
                    idx = vcf_id_map[sid]
                    if pheno not in pop_indices:
                        pop_indices[pheno] = []
                    pop_indices[pheno].append(idx)
                    found_count += 1
            
            # Save to self.subpops
            self.subpops['Whole'] = np.arange(len(self.sample_ids))
            for p, idxs in pop_indices.items():
                self.subpops[f"Subpop_{p}"] = np.array(idxs)
                
            print(f"Matched {found_count} phenotypes. Subgroups found: {list(pop_indices.keys())}")
            
        except Exception as e:
            print(f"Error loading phenotypes: {e}")
            print("Proceeding with 'Whole' population only.")
            self.subpops['Whole'] = np.arange(len(self.sample_ids))

    # ==========================================
    # 1. ROH Calling (The Foundation)
    # ==========================================
    def call_rohs_per_individual(self):
        """
        Identifies which SNPs are in an ROH for each individual.
        Logic: Find contiguous 'homozygous' stretches > roh_cutoff.
        """
        print("\n--- Calling ROHs per Individual ---")
        n_ind, n_snp = self.hap1.shape
        
        # Initialize boolean matrix (False = Not in ROH)
        self.roh_matrix = np.zeros((n_ind, n_snp), dtype=bool)
        
        # Calculate homozygous state (0=Het, 1=Hom)
        # Note: We check if hap1 == hap2. 
        # (Beagle phasing allows this; raw unphased data would be less accurate)
        is_hom = (self.hap1 == self.hap2)
        
        chrom_unique = self.map_df['chrom'].unique()
        chrom_unique = chrom_unique[chrom_unique != 0]

        # Iterate over chromosomes to prevent ROHs crossing chr boundaries
        for chrom in tqdm(chrom_unique, desc="Processing Chromosomes"):
            chr_mask = (self.map_df['chrom'] == chrom).values
            chr_indices = np.where(chr_mask)[0]
            positions = self.map_df.loc[chr_mask, 'pos_mb'].values
            
            # Extract hom status for this chromosome (N_ind x N_chr_snps)
            sub_hom = is_hom[:, chr_mask]
            
            # For each individual, find runs
            for i in range(n_ind):
                # Identify changes in state to find runs
                # Add False at ends to capture boundary runs
                row = sub_hom[i, :]
                bounded = np.concatenate(([False], row, [False]))
                
                # Find where values change
                diffs = np.diff(bounded.astype(int))
                run_starts = np.where(diffs == 1)[0]
                run_ends = np.where(diffs == -1)[0]
                
                # Check lengths
                for start, end in zip(run_starts, run_ends):
                    # -1 because end is exclusive index in diff, but inclusive in positions map
                    # Use positions to check physical length
                    length_mb = positions[end-1] - positions[start]
                    
                    if length_mb >= self.roh_cutoff:
                        # Mark these SNPs as ROH in the global matrix
                        global_start = chr_indices[start]
                        #global_end = chr_indices[end] # exclusive for slicing
                        global_end = global_start + (end - start)
                        self.roh_matrix[i, global_start:global_end] = True

        print("ROH Calling Complete.")

    # ==========================================
    # 2. Island Detection & 3. Genome Stats
    # ==========================================
    def analyze_subpopulations(self):
        """
        Performs the analysis requested:
        1. ROH Islands (Top 5%) -> SAVES TO CSV
        2. Frequency Plots -> SAVES TO PNG
        3. Genome Proportion Stats -> SAVES TO CSV
        """
        print("\n--- Analyzing Subpopulations & Filtering Islands ---")
        
        # 1. Calculate Global ROH Length Statistics (needed for the filter)
        # We need the length (in SNPs) of every single ROH detected in the population
        all_roh_lengths = []
        
        # Scan a subset or all individuals to get a stable mean/SD
        # Optimization: We don't need to check every single person if N is large, 
        # but for N=576 it's fast enough to check all.
        print("Calculating global ROH stats to define size cutoff...")
        for i in range(self.roh_matrix.shape[0]):
            row = self.roh_matrix[i, :]
            if not np.any(row): continue
            
            # Find lengths of True blocks
            bounded = np.concatenate(([False], row, [False]))
            diffs = np.diff(bounded.astype(int))
            run_lengths = np.where(diffs == -1)[0] - np.where(diffs == 1)[0]
            all_roh_lengths.extend(run_lengths)

        if not all_roh_lengths:
            print("No ROHs found in the entire population. Skipping analysis.")
            return

        mean_len = np.mean(all_roh_lengths)
        std_len = np.std(all_roh_lengths)
        
        # Manuscript Rule: Remove regions < Mean - 2*SD
        min_snp_cutoff = mean_len - (2 * std_len)
        
        # Safety: If SD is huge, cutoff might be negative. Ensure at least 3 SNPs.
        if min_snp_cutoff < 3:
            min_snp_cutoff = 3 
            
        print(f"Global ROH Stats: Mean={mean_len:.2f} SNPs, SD={std_len:.2f}")
        print(f"Island Size Cutoff (Mean - 2SD): {min_snp_cutoff:.2f} SNPs")

        # ---------------------------------------------------------

        results_stats = []
        all_islands = [] # Store all found islands here
        
        # Setup Plotting
        plt.figure(figsize=(15, 6))
        colors = ['black', 'red', 'blue', 'green']
        
        for i, (pop_name, indices) in enumerate(self.subpops.items()):
            print(f"\nAnalyzing: {pop_name} (N={len(indices)})")
            
            if len(indices) == 0: continue

            # --- A. Frequency & Islands ---
            # Extract sub-matrix
            sub_roh = self.roh_matrix[indices, :]
            
            # Frequency per SNP (0.0 to 1.0)
            freq_per_snp = np.mean(sub_roh, axis=0)
            
            # Detect Islands (Top 5%)
            # Logic: Top 5% of SNPs with *non-zero* frequency, or top 5% overall?
            # Manuscript implies "regions with high ROH frequency... keeping top 5% SNPs"
            # We take the 95th percentile of the frequency distribution
            threshold = np.percentile(freq_per_snp, 95)
            
            # Identify SNPs above threshold
            is_high = freq_per_snp >= threshold
            
            # Group into regions (Islands)
            # (Simple grouping of contiguous high-freq SNPs)
            #islands = []
            
            # Add False to edges to find boundaries
            bounded = np.concatenate(([False], is_high, [False]))
            diffs = np.diff(bounded.astype(int))
            starts = np.where(diffs == 1)[0]
            ends = np.where(diffs == -1)[0]
            
            print(f"  > Found {len(starts)} potential islands (Threshold freq={threshold:.3f})")
            
            # Filter Islands (Mean - 2SD rule from manuscript)
            # First, calculate stats for "standard" ROH windows to get the Mean/SD baseline
            # For simplicity here, we use the island lengths themselves or a standard heuristic
            # The manuscript uses "average number of SNPs included in a ROH window"
            # We'll skip the complex cross-referencing and just report the islands for now.


            # Store Island Details
            kept_count = 0
            for s, e in zip(starts, ends):
                # e is exclusive index
                n_snps = e - s
                # THE FILTER STEP
                if n_snps >= min_snp_cutoff:
                    chrom = self.map_df.iloc[s]['chrom']
                    start_mb = self.map_df.iloc[s]['pos_mb']
                    end_mb = self.map_df.iloc[e-1]['pos_mb']
                    mean_freq = np.mean(freq_per_snp[s:e])
                
                    all_islands.append({
                        'Population': pop_name,
                        'Chromosome': chrom,
                        'Start_Mb': start_mb,
                        'End_Mb': end_mb,
                        'Num_SNPs': n_snps,
                        'Mean_Frequency': mean_freq
                    })
                    kept_count += 1
            
            print(f"  > Raw clusters: {len(starts)} -> After Filter: {kept_count} Islands")

            # --- B. Calculate Genome Proportion (Froh) ---
            # Sum of ROH lengths / Total Genome Length
            # Total Genome Length (Autosomes) ~ 2500 Mb usually, but let's calc from map
            #total_genome_mb = self.map_df.groupby('chrom')['pos_mb'].max().sum()
            
            # Per individual in this subpop
            froh_vals = []
            for idx in range(sub_roh.shape[0]):
                # Indices where this individual has ROH
                # We need to sum physical distances, not just SNP counts
                # Approximation: Sum of (Num SNPs in ROH * Avg SNP spacing)
                # More Accurate: We already filtered by Mb, so we can sum the segments.
                # Fast Approximation: Count ROH SNPs * (TotalMb / TotalSNPs)
                # Only valid if density is uniform.
                
                # Let's use SNP count fraction as proxy for genome fraction
                # Froh = (Count of ROH SNPs) / (Total SNPs)
                # This is a standard proxy when exact lengths aren't stored per person
                froh = np.sum(sub_roh[idx, :]) / len(self.map_df)
                froh_vals.append(froh)
            
            mean_froh = np.mean(froh_vals)
            std_froh = np.std(froh_vals)
            results_stats.append({'Pop': pop_name, 'Mean_Froh': mean_froh, 'SD_Froh': std_froh})
            print(f"  > Genome Proportion (Froh): {mean_froh:.4f} +/- {std_froh:.4f}")

            # --- C. Add to Plot ---
            # To avoid plotting 50k points, we can bin or smooth, 
            # but usually genome-wide plots just plot all points or a moving average.
            # We'll plot a simple moving average for visualization clarity
            window_size = 100
            smoothed_freq = np.convolve(freq_per_snp, np.ones(window_size)/window_size, mode='same')
            
            # Create a cumulative position for X-axis (Manhattan style)
            # But for simple overlay, we just use index
            plt.plot(smoothed_freq, label=f"{pop_name} (N={len(indices)})", color=colors[i % len(colors)], alpha=0.7)

        # Finalize Plot
        plt.title("Genome-Wide ROH Frequency (Island Cutoff > {int(min_snp_cutoff)} SNPs)")
        plt.xlabel("SNP Index (Genome Order)")
        plt.ylabel("Frequency of ROH")
        plt.legend()
        plt.tight_layout()
        plt.savefig(f"Unfiltered_ROH_Frequency_Plot.{roh_mb_cutoff}.png")
        print(f"\nPlot saved to 'Unfiltered_ROH_Frequency_Plot.{roh_mb_cutoff}.png'")
        
        # Save Stats
        pd.DataFrame(results_stats).to_csv(f"Unfiltered_ROH_Subpop_Stats.{roh_mb_cutoff}.csv", index=False)
        print(f"Stats saved to 'Unfiltered_ROH_Subpop_Stats.{roh_mb_cutoff}.csv'")

        # Save ISLANDS Table (The detailed regions)
        islands_df = pd.DataFrame(all_islands)
        islands_df.to_csv(f"Filtered_ROH_Islands_Detailed.{roh_mb_cutoff}.csv", index=False)
        print(f"A List of detected genomic regions saved to 'Filtered_ROH_Islands_Detailed.{roh_mb_cutoff}.csv'")

# ==========================================
# EXECUTION
# ==========================================
import sys

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python ROH_analysis.py <filtered.norm.phased.vcf.gz> <phenotypes.txt> <roh_mb_cutoff>")
        sys.exit(1)

    # Settings
    VCF_FILE = sys.argv[1]
    PHENO_FILE = sys.argv[2]
    roh_mb_cutoff = float(sys.argv[3])
    
    # Run
    # 1. Init
    analysis = ROHAnalysis(VCF_FILE, PHENO_FILE, roh_mb_cutoff)
    
    # 2. Load
    analysis.load_data()
    
    # 3. Call ROHs (The heavy lifting)
    analysis.call_rohs_per_individual()
    
    # 4. Analyze Subpops (Islands, Plots, Stats)
    analysis.analyze_subpopulations()
