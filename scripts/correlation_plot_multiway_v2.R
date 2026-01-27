args = commandArgs(TRUE);
require(ggplot2);
require(reshape2);
require(GGally);

## Assign arguments
file1 <- args[1];
file2 <- args[2];
file3 <- args[3];
out_prefix <- args[4];

print("### --- 1. Load data ---");
diag <- read.csv(file1, header=TRUE, comment.char="", stringsAsFactors=FALSE);
names(diag)[names(diag) == "F_Standard"] <- "D_STD"
names(diag)[names(diag) == "F_ROH"] <- "D_ROH"

het  <- read.table(file2, header=TRUE, comment.char="", stringsAsFactors=FALSE);
names(het)[names(het) == "F"] <- "F_SNP"

froh <- read.table(file3, header=TRUE, comment.char="", stringsAsFactors=FALSE);
names(froh)[names(froh) == "F_ROH"] <- "F_ROH"

print("### --- 2. Merge on IID ---");
list_of_dfs <- list(het,froh,diag)
dat <- Reduce(function(x, y) merge(x, y, by = "IID", all = TRUE), list_of_dfs)

#print("### --- 3. Compute O(HET) - E(HET) ---");
#dat$HET_diff <- dat$O.HET - dat$E.HET;

print("### --- 4. Select variables for correlation ---");
vars <- dat[, c("F_SNP", "F_ROH", "D_ROH", "D_STD")];

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
