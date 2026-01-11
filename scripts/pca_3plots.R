args=(commandArgs(TRUE));
require(ggplot2);
require(gridExtra);

## Read inputs
eigenvec <- read.table(paste(args[1],"eigenvec",args[2],sep="."), header = TRUE, comment.char="");
eigenval <- scan(paste(args[1], "eigenval", sep="."), quiet = TRUE);

## prep labels
var_explained <- round(100 * eigenval / sum(eigenval), 2);
lab_PC1 <- paste0("PC1 (", var_explained[1], "%)");
lab_PC2 <- paste0("PC2 (", var_explained[2], "%)");
lab_PC3 <- paste0("PC3 (", var_explained[3], "%)");

# Plots
color_col <- args[3]
eigenvec[[color_col]] <-  as.factor(eigenvec[[color_col]]);
plot1 <- ggplot(eigenvec, aes(x = PC1, y = PC2, col = .data[[color_col]])) + geom_point() + labs(title = "PCA Plot", x = lab_PC1, y = lab_PC2);
plot2 <- ggplot(eigenvec, aes(x = PC1, y = PC3, col = .data[[color_col]])) + geom_point() + labs(title = "PCA Plot", x = lab_PC1, y = lab_PC3);
plot4 <- ggplot(eigenvec, aes(x = PC2, y = PC3, col = .data[[color_col]])) + geom_point() + labs(title = "PCA Plot", x = lab_PC2, y = lab_PC3);
combined_plot <- grid.arrange(plot1, plot2, plot4, nrow = 1);
ggsave(file = args[4], combined_plot, width = 16, height = 4, dpi = 400);
