args = commandArgs(TRUE);
require(ggplot2);
require(reshape2);
require(GGally);

## Assign arguments
file1 <- args[1];
file2 <- args[2];
out_prefix <- args[3];

print("### --- 1. Load data ---");
homo <- read.table(file1, header=TRUE, comment.char="", stringsAsFactors=FALSE);
het  <- read.table(file2, header=TRUE, comment.char="", stringsAsFactors=FALSE);

print("### --- 2. Merge on IID ---");
dat <- merge(homo, het, by="IID");

print("### --- 3. Compute O(HET) - E(HET) ---");
dat$HET_diff <- dat$O.HET - dat$E.HET;

print("### --- 4. Select variables for correlation ---");
vars <- dat[, c("KB", "KBAVG", "HET_diff", "F")];

print("### --- 5. Correlation matrix ---");
cor_matrix <- cor(vars, use="pairwise.complete.obs");
melted_cor <- melt(cor_matrix);

print("### --- 6. Heatmap ---");
p_heatmap <- ggplot(melted_cor, aes(Var1, Var2, fill=value)) +
    geom_tile(color="white") +
    scale_fill_gradient2(low="blue", high="red", mid="white", midpoint=0) +
    geom_text(aes(label=sprintf("%.2f", value)), size=4) +
    theme_minimal(base_size=14) +
    labs(title="Correlation Matrix", x="", y="", fill="Correlation");
ggsave(paste0(out_prefix, ".correlation_heatmap.png"),
       p_heatmap, width=6, height=6, dpi=400);

print("### --- 7. Scatterplot matrix ---");
p_pairs <- ggpairs(vars,
           upper=list(continuous="cor"),
           diag=list(continuous="densityDiag"),
           lower=list(continuous="smooth"));
ggsave(paste0(out_prefix, ".pairplot.png"), p_pairs, width=10, height=10, dpi=300);
