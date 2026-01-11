import allel
import numpy as np
import pandas as pd
import sys
from tqdm import tqdm  # Progress bar

class ROHAnalysis:
    def __init__(self, vcf_path, roh_mb_cutoff=1.0, roh_threshold=3):
        self.vcf_path = vcf_path
        self.roh_cutoff = roh_mb_cutoff
        self.roh_threshold = roh_threshold # Statistical threshold (Mean - k*SD)
        
        # Data
        self.map_df = None
        self.hap1 = None
        self.hap2 = None
        self.sample_ids = None
        self.windows = [] # Will store tuples: (chr, start_idx, end_idx, num_snps)
        self.G_matrix = None # Final Relationship Matrix

    # ==========================================
    # STAGE 1: Corrected Data Loading
    # ==========================================
    def load_data(self):
        print(f"Reading VCF: {self.vcf_path}")
        
        # 1. Read VCF
        # We explicitly assume phase is present because you ran Beagle.
        try:
            callset = allel.read_vcf(self.vcf_path, fields=['calldata/GT', 'variants/CHROM', 'variants/POS', 'samples'])
        except Exception as e:
            print(f"Error reading VCF: {e}")
            sys.exit(1)
            
        if callset is None:
            print("Error: VCF read failed.")
            sys.exit(1)

        # 2. Extract Sample IDs
        self.sample_ids = callset['samples']
        print(f"Loaded {len(self.sample_ids)} samples.")

        # 3. Extract Map & Convert to Mb
        chroms = callset['variants/CHROM']
        positions_bp = callset['variants/POS']
        
        # Clean Chromosomes (remove 'chr' if present)
        chroms = [str(c).replace('chr', '') for c in chroms]
        chroms = np.array([int(c) if c.isdigit() else 0 for c in chroms])

        self.map_df = pd.DataFrame({
            'chrom': chroms,
            'pos_mb': positions_bp / 1e6,  # Convert BP to Mb
            'original_idx': np.arange(len(chroms))
        })
        
        print(f"Loaded Map: {len(self.map_df)} SNPs.")

        # 4. Extract Haplotypes (Assuming Phased Input)
        gt = allel.GenotypeArray(callset['calldata/GT'])
        
        # Skip the .is_phased check. We trust the Beagle output.
        # hap1 = Paternal (index 0), hap2 = Maternal (index 1)
        self.hap1 = np.array(gt[:, :, 0], dtype=np.int8) # Shape: (Variants, Samples)
        self.hap2 = np.array(gt[:, :, 1], dtype=np.int8) # Shape: (Variants, Samples)
        
        # Transpose to (Samples, Variants) for faster row access in later stages
        self.hap1 = self.hap1.T
        self.hap2 = self.hap2.T
        
        print(f"Haplotypes Loaded: {self.hap1.shape[0]} Inds x {self.hap1.shape[1]} SNPs")

    # ==========================================
    # STAGE 2: Window Definition (Geometry)
    # ==========================================
    def define_windows(self):
        """
        Replicates the C++ logic:
        For every SNP 'i', find the first SNP 'j' such that distance(i, j) > cutoff.
        """
        print("\n--- Stage 2: Defining ROH Windows ---")
        
        # We process one chromosome at a time to avoid calculating cross-chr distances
        unique_chroms = self.map_df['chrom'].unique()
        unique_chroms = unique_chroms[unique_chroms != 0] # Exclude unmapped
        
        all_windows = []
        
        for chrom in unique_chroms:
            # Extract positions for this chromosome
            chr_mask = self.map_df['chrom'] == chrom
            chr_indices = self.map_df[chr_mask].index.values
            positions = self.map_df.loc[chr_mask, 'pos_mb'].values
            
            n_snps = len(positions)
            
            # Use a sliding window approach
            # The C++ code nested loop is O(N^2) which is slow.
            # We can use a 'two-pointer' approach for O(N) speed.
            
            right = 0
            for left in range(n_snps):
                # Move the right pointer forward until the distance exceeds cutoff
                # C++ Logic: "Determine where first one that is over threshold is reached"
                while right < n_snps:
                    dist = positions[right] - positions[left]
                    if dist > self.roh_cutoff:
                        break # Found the breaker
                    right += 1
                
                # If right reached the end, we might still have a window, 
                # but let's stick to C++ logic: "if(j < rows)" check
                if right < n_snps:
                    # Map back to original global indices
                    global_start = chr_indices[left]
                    global_end = chr_indices[right] # The C++ code includes the boundary breaker
                    
                    num_snp_in_window = (global_end - global_start) + 1
                    
                    # Store Window: (Chrom, Start_Global_Idx, End_Global_Idx, Num_SNPs)
                    all_windows.append((chrom, global_start, global_end, num_snp_in_window))
        
        self.windows = all_windows
        print(f"Total Raw Windows Found: {len(self.windows)}")


    # ==========================================
    # STAGE 3: Statistical Filtering
    # ==========================================
    def filter_windows(self):
        """
        Calculates the distribution of 'Number of SNPs' in ROH windows.
        Removes windows that fall below Mean - (Threshold * SD).
        """
        print("\n--- Stage 3: Filtering Windows ---")
        if not self.windows:
            print("No windows to filter.")
            return

        # Extract 'Num SNPs' (4th element in tuple)
        lengths = np.array([w[3] for w in self.windows])
        
        mean_len = np.mean(lengths)
        std_len = np.std(lengths)
        
        # Calculate Cutoff (rounding up as per C++: + 0.5 cast to int)
        # Cutoff = Mean - (Threshold * Stdev)
        cutoff = int(mean_len - (self.roh_threshold * std_len) + 0.5)
        
        print(f"Stats: Mean={mean_len:.2f}, SD={std_len:.2f}")
        print(f"Cutoff Calculated: {cutoff} SNPs")
        
        # Apply Filter
        original_count = len(self.windows)
        self.windows = [w for w in self.windows if w[3] >= cutoff]
        
        print(f"Windows filtered: {original_count} -> {len(self.windows)} remaining.")

    # ==========================================
    # STAGE 4: The Kernel (Haplotype Matching)
    # ==========================================
    def compute_relationship_matrix(self):
        """
        Builds the ROH Relationship Matrix.
        Optimized using NumPy broadcasting to avoid N^2 Python loops.
        """
        print("\n--- Stage 4: Computing Relationship Matrix ---")
        
        n_ind = self.hap1.shape[0]
        # Initialize G matrix (N x N)
        self.G_matrix = np.zeros((n_ind, n_ind), dtype=np.float32)
        
        # Pre-allocate a 2N x 2N buffer for the haplotype matches to save memory alloc time
        # This represents matches between all gametes (2 per individual)
        # But actually, we can do it smarter:
        # We need sum of: (Pat_i==Pat_j) + (Pat_i==Mat_j) + (Mat_i==Pat_j) + (Mat_i==Mat_j)
        
        print(f"Processing {len(self.windows)} windows...")
        
        for w_idx, (chrom, start_idx, end_idx, n_snps) in enumerate(tqdm(self.windows)):
            
            # 1. Extract Haplotypes for this window
            # Slice: all individuals, SNPs from start to end (inclusive in C++, so end+1 for Python slice)
            # Shapes: (N_ind, Window_Size)
            h1_block = self.hap1[:, start_idx : end_idx + 1]
            h2_block = self.hap2[:, start_idx : end_idx + 1]
            
            # 2. Assign Unique IDs to Haplotypes
            # We stack them to get shape (2*N_ind, Window_Size)
            # Rows 0..N-1 are Paternal, N..2N-1 are Maternal
            combined_haps = np.vstack([h1_block, h2_block])
            
            # Find unique rows and assign integer codes
            # axis=0 means we look for unique haplotypes (rows)
            # return_inverse gives us the ID for each row in combined_haps
            _, indices = np.unique(combined_haps, axis=0, return_inverse=True)
            
            # 'indices' is a vector of length 2*N_ind containing the Haplotype ID for each gamete
            # Split back into Paternal and Maternal IDs
            pat_ids = indices[:n_ind] # Shape (N,)
            mat_ids = indices[n_ind:] # Shape (N,)
            
            # 3. Compute Matches (Vectorized)
            # We need to add to G[i,j]:
            #   (Pat[i] == Pat[j]) + (Pat[i] == Mat[j]) + (Mat[i] == Pat[j]) + (Mat[i] == Mat[j])
            
            # Create boolean match matrices (N x N)
            # Broadcasting: pat_ids[:, None] is (N,1), pat_ids[None, :] is (1,N)
            # Result is (N,N) boolean matrix
            pp_match = (pat_ids[:, None] == pat_ids[None, :])
            pm_match = (pat_ids[:, None] == mat_ids[None, :])
            mp_match = (mat_ids[:, None] == pat_ids[None, :])
            mm_match = (mat_ids[:, None] == mat_ids[None, :])
            
            # Sum them up (converting bool to float implicitly or explicitly)
            # The C++ formula is: (Sum of 4 comparisons) / 2
            total_match = (pp_match.astype(np.float32) + 
                           pm_match.astype(np.float32) + 
                           mp_match.astype(np.float32) + 
                           mm_match.astype(np.float32)) / 2.0
            
            # Accumulate into Global Matrix
            self.G_matrix += total_match

        # ==========================================
        # STAGE 5: Normalization & Output
        # ==========================================
        print("\n--- Stage 5: Normalization ---")
        # Divide by number of windows
        self.G_matrix /= len(self.windows)
        print("Matrix construction complete.")

    def save_matrix(self, roh_mb_cutoff, roh_threshold):
        filename_prefix=f"{output_dir}/ROHRM.rohMinSize_{roh_mb_cutoff}.rohThreshold_{roh_threshold}"
        print(f"Saving Matrix and IDs to {filename_prefix}...")

        # 1. Save the Matrix
        np.savetxt(f"{filename_prefix}.txt", self.G_matrix, fmt='%.6f', delimiter=' ')

        # 2. Save the IDs (matching PLINK format: FamilyID IndividualID)
        # Since VCF usually just has one ID, we repeat it for FID and IID
        with open(f"{filename_prefix}.id", "w") as f:
            for sample_id in self.sample_ids:
                # PLINK .id files usually have 2 columns: FID IID.
                # We will write: ID ID
                f.write(f"{sample_id}\t{sample_id}\n")

        print("Done.")


# ==========================================
# FINAL EXECUTION SCRIPT
# ==========================================
import sys

if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: python ROHRM_Creator.py <filtered.norm.phased.vcf.gz> <roh_mb_cutoff> <roh_threshold> <output_dir>")
        sys.exit(1)

    # Update these paths to your specific files
    VCF_FILE = sys.argv[1]
    roh_mb_cutoff = float(sys.argv[2])
    roh_threshold = float(sys.argv[3])
    output_dir = sys.argv[4]
    
    # 1. Initialize
    # roh_threshold=3 corresponds to "Mean - 3*SD" from the C++ input
    analysis = ROHAnalysis(VCF_FILE, roh_mb_cutoff, roh_threshold)
    
    # 2. Load
    analysis.load_data()
    
    # 3. Define Windows
    analysis.define_windows()
    
    # 4. Filter Windows
    analysis.filter_windows()
    
    # 5. Compute Kernel
    if len(analysis.windows) > 0:
        analysis.compute_relationship_matrix()
        analysis.save_matrix(roh_mb_cutoff, roh_threshold)
    else:
        print("No windows survived filtering. Cannot compute matrix.")


