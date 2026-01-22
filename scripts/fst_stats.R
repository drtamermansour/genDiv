library(dplyr)

# --- CONFIGURATION ---
file_sex      <- "divStats/filtered.LD_prune.fst_sex.fst.summary"
file_gait     <- "divStats/filtered.LD_prune.fst_gait.fst.summary"
file_book_all <- "divStats/filtered.LD_prune.fst_bookSize.fst.summary"
file_book_trotter <- "divStats/filtered.LD_prune.fst_bookSize.Trotter.fst.summary"
file_book_pacer   <- "divStats/filtered.LD_prune.fst_bookSize.Pacer.fst.summary"

# --- FIXED FUNCTION ---
load_fst <- function(filepath, category_label) {
  if (file.exists(filepath)) {
    # Check if file is empty before reading
    if (file.info(filepath)$size == 0) {
      warning(paste("File is empty:", filepath))
      return(NULL)
    }

    # FIX: Use colClasses to force the first two columns to be Character.
    # This prevents "F" (Female) from becoming FALSE.
    df <- read.table(filepath, 
                     comment.char="#", 
                     header=FALSE, 
                     stringsAsFactors=FALSE,
                     colClasses = c("character", "character", "numeric", "numeric"))
    
    colnames(df) <- c("POP1", "POP2", "Fst", "SE")
    df$Category <- category_label
    return(df)
  } else {
    warning(paste("File not found:", filepath))
    return(NULL)
  }
}

# --- 1. LOAD DATA ---
data_list <- list()

# Load and label
data_list[[1]] <- load_fst(file_sex,      "Control_Sex")
data_list[[2]] <- load_fst(file_gait,     "Gait_Comparison")
data_list[[3]] <- load_fst(file_book_all, "BookSize_Global")
data_list[[4]] <- load_fst(file_book_trotter, "BookSize_Within_Trotters")
data_list[[5]] <- load_fst(file_book_pacer,   "BookSize_Within_Pacers")

# Filter out any NULLs (in case a file was missing/empty) and bind
data_list <- data_list[!sapply(data_list, is.null)]
final_df <- bind_rows(data_list)

# --- 2. STATISTICAL ANALYSIS ---

# Z-Score
final_df$Z_score <- final_df$Fst / final_df$SE

# P-Value (One-tailed)
final_df$P_value <- pnorm(final_df$Z_score, lower.tail = FALSE)

# Bonferroni Correction
n_tests <- nrow(final_df) 
alpha_adj <- 0.05 / n_tests

final_df$Significant <- final_df$P_value < alpha_adj
final_df$Significance_Label <- ifelse(final_df$Significant, "***", "ns")

# --- 3. OUTPUT ---
output_table <- final_df %>%
  mutate(
    # Calculate Adjusted P-value (P_raw * n_tests)
    # We cap it at 1.0 because a probability cannot exceed 100%
    P_Value_Adj = pmin(P_value * n_tests, 1.0),
    
    # Format for display
    P_Raw_Scientific = formatC(P_value, format = "e", digits = 2),
    P_Adj_Scientific = formatC(P_Value_Adj, format = "e", digits = 2),
    
    # Significance Check (Compare Adjusted P to 0.05)
    Significance_Label = ifelse(P_Value_Adj < 0.05, "***", "ns")
  ) %>%
  select(Category, POP1, POP2, Fst, SE, Z_score, P_Raw_Scientific, P_Adj_Scientific, Significance_Label) %>%
  arrange(Category)

cat("\n=== FST STATISTICAL ANALYSIS (Adjusted P-Values) ===\n")
cat("Total Tests (Correction Factor):", n_tests, "\n")
print(as.data.frame(output_table))

write.csv(output_table, "Fst_Analysis_Results_Adjusted.csv", row.names = FALSE)
