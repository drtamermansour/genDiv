args = commandArgs(TRUE);
require(ggplot2);
require(reshape2);
require(GGally);

## Assign arguments
file1 <- args[1];
file2 <- args[2];
file3 <- args[3];
out_prefix <- args[4];

# Function to create a sorted, combined key for a row
create_sorted_key <- function(col1, col2) {
	  apply(cbind(as.character(col1), as.character(col2)), 1, function(x) paste(sort(x), collapse = "_"))
}

print("### --- 1. Load data ---");
diag <- read.csv(file1, header=TRUE, comment.char="", stringsAsFactors=FALSE);
diag$key <- create_sorted_key(diag$ID1, diag$ID2)
names(diag)[names(diag) == "Kinship_Std"] <- "G_STD"
names(diag)[names(diag) == "Kinship_ROH"] <- "G_ROH"

kingkin_wIBS  <- read.table(file2, header=TRUE, comment.char="", stringsAsFactors=FALSE);
kingkin_wIBS$key <- create_sorted_key(kingkin_wIBS$IID1, kingkin_wIBS$IID2)
names(kingkin_wIBS)[names(kingkin_wIBS) == "KINSHIP"] <- "King_kin"


euclDist <- read.table(file3, header=TRUE, comment.char="", stringsAsFactors=FALSE);
euclDist$key <- create_sorted_key(euclDist$IID1, euclDist$IID2)
names(euclDist)[names(euclDist) == "PCA_EUCLIDEAN_DIST"] <- "Eucl_dist"

print("### --- 2. Merge on IID ---");
list_of_dfs <- list(diag,kingkin_wIBS,euclDist)
dat <- Reduce(function(x, y) merge(x, y, by = "key", all = FALSE), list_of_dfs)

#print("### --- 3. Compute O(HET) - E(HET) ---");
#dat$HET_diff <- dat$O.HET - dat$E.HET;

print("### --- 4. Select variables for correlation ---");
#vars <- dat[, c("G_STD", "G_ROH", "King_kin","IBS", "Eucl_dist")];
vars <- dat[, c("G_STD", "G_ROH", "King_kin","IBS")];

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
       p_heatmap, width=8, height=8, dpi=400);

print("### --- 7. Scatterplot matrix ---");
p_pairs <- ggpairs(vars,
           upper=list(continuous="cor"),
           diag=list(continuous="densityDiag"),
           lower=list(continuous="smooth"));
ggsave(paste0(out_prefix, ".pairplot.png"), p_pairs, width=12, height=12, dpi=300);
