library(ggplot2)
library(dplyr)
library(tidyr)

# --- CONFIGURATION ---
base_path <- "LD_pruned/USTA_Diversity_Study.remap.refAlleles.dedup.plink1.filtered.LD_prune.freq_stats" 
gaits <- c("Trotter", "Pacer")
book_sizes <- c("LOW", "MEDIUM", "HIGH")

# --- ROBUST DATA LOADER (Same as before) ---
load_data <- function(filename, group_label) {
  if (file.exists(filename)) {
    df <- read.table(filename, header=TRUE, comment.char="", stringsAsFactors=FALSE)
    if("A_e" %in% colnames(df)) {
      return(data.frame(Ae = df$A_e, Group = group_label))
    }
  }
  return(NULL)
}

# --- 1. LOAD AND PREPARE DATA ---
# We load the fine-scale data (Gait + Book Size) as it contains the most info
all_data <- data.frame()

for (g in gaits) {
  for (b in book_sizes) {
    fname <- paste0(base_path, ".", g, "_", b, ".afreq.Ae")
    label <- paste0(g, "_", b) 
    temp <- load_data(fname, label)
    if(!is.null(temp)) all_data <- rbind(all_data, temp)
  }
}

# Split the "Group" column into "Gait" and "BookSize" for plotting
plot_data <- all_data %>%
  separate(Group, into = c("Gait", "BookSize"), sep = "_")

# ENFORCE ORDER: Make sure "LOW" comes before "MEDIUM" (otherwise alphabetical puts HIGH first)
plot_data$BookSize <- factor(plot_data$BookSize, levels = c("LOW", "MEDIUM", "HIGH"))

# --- 2. GENERATE PLOT: VIOLIN PLOT WITH MEANS ---

# We use a violin plot to show density, and overlay the Mean (+/- SD)
# This matches your T-test logic (which compares Means, not Medians)

p <- ggplot(plot_data, aes(x = BookSize, y = Ae, fill = BookSize)) +
  # A. The Violin (Distribution Shape)
  geom_violin(trim = FALSE, alpha = 0.6, color = NA) +
  
  # B. The Boxplot (Inside, for Median/IQR reference)
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA, alpha = 0.8) +
  
  # C. The Mean (Red Dot - distinct from median line)
  stat_summary(fun = mean, geom = "point", shape = 18, size = 3, color = "red") +
  
  # D. Faceting (Split by Gait)
  facet_wrap(~Gait) +
  
  # E. Aesthetics
  scale_fill_manual(values = c("LOW"="#619CFF", "MEDIUM"="#00BA38", "HIGH"="#F8766D")) +
  labs(
    title = "Genetic Diversity (Ae) Stratified by Book Size",
    subtitle = "Red diamond indicates the Mean; Boxplot indicates Median",
    y = expression(paste("Effective Number of Alleles (", A[e], ")")),
    x = "Book Size (Breeding Intensity)"
  ) +
  theme_bw() +
  theme(
    legend.position = "none", # Hide legend since X-axis explains it
    text = element_text(size = 14),
    strip.background = element_rect(fill = "gray90"), # Gray header for facets
    panel.grid.minor = element_blank()
  )

# --- 3. SAVE AND VIEW ---
#print(p)

# Save for publication (TIFF/PDF are standard)
ggsave("divStats/Figure_Ae_BookSize.tiff", plot = p, width = 8, height = 6, dpi = 300)
ggsave("divStats/Figure_Ae_BookSize.pdf", plot = p, width = 8, height = 6)
