library(effsize)

# --- CONFIGURATION ---
gaits <- c("Trotter", "Pacer")
book_sizes <- c("LOW", "MEDIUM", "HIGH")
base_path <- "LD_pruned/USTA_Diversity_Study.remap.refAlleles.dedup.plink1.filtered.LD_prune.freq_stats" 

# --- ROBUST DATA LOADER ---
load_data <- function(filename, group_label) {
  if (file.exists(filename)) {
    # FIX 1: comment.char="" allows reading the line starting with #
    # FIX 2: header=TRUE allows us to use column names
    df <- read.table(filename, header=TRUE, comment.char="", stringsAsFactors=FALSE)
    
    # FIX 3: Select column by NAME ("A_e") instead of index (V7)
    # Note: R might sanitize "#CHROM" to "X.CHROM" or similar, but "A_e" will remain "A_e"
    if("A_e" %in% colnames(df)) {
      return(data.frame(Ae = df$A_e, Group = group_label))
    } else {
      stop(paste("Error: Column 'A_e' not found in file:", filename, 
                 "\nAvailable columns:", paste(colnames(df), collapse=", ")))
    }
  } else {
    warning(paste("File not found:", filename))
    return(NULL)
  }
}

# --- 1. GAIT COMPARISON ---
cat("\n=== ANALYSIS 1: Trotter vs Pacer ===\n")
data_gait <- rbind(
  load_data(paste0(base_path, ".Trotter.afreq.Ae"), "Trotter"),
  load_data(paste0(base_path, ".Pacer.afreq.Ae"), "Pacer")
)

# T-test (Pairwise comparison)
t_res <- t.test(Ae ~ Group, data = data_gait)
print(t_res)

# Cohen's d (Effect Size)
d_res <- cohen.d(data_gait$Ae, data_gait$Group)
cat("Cohen's d (Effect Size):", d_res$estimate, "\n")


# --- 2. BOOK SIZE COMPARISON ---
cat("\n=== ANALYSIS 2: Book Size Comparisons ===\n")
data_book <- data.frame()

# Load all 6 combinations
for (g in gaits) {
  for (b in book_sizes) {
    fname <- paste0(base_path, ".", g, "_", b, ".afreq.Ae")
    label <- paste0(g, "_", b) 
    temp <- load_data(fname, label)
    if(!is.null(temp)) data_book <- rbind(data_book, temp)
  }
}

data_book$Group <- factor(data_book$Group)

# Pairwise T-test with Bonferroni Correction
pairwise_res <- pairwise.t.test(data_book$Ae, data_book$Group, 
                                p.adjust.method = "bonferroni", 
                                pool.sd = FALSE) 
print(pairwise_res)
