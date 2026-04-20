## Unified correlation heatmap + scatterplot matrix.
##
## Supersedes the four legacy correlation_plot_multiway variants (v1/v2/v2e/v3).
## Each former variant maps to a --mode.
##
## Usage:
##   Rscript correlation_plot.R --mode basic <homo_file> <het_file> <out_prefix>
##     # legacy v1: merges homo+het by IID, correlates KB, KBAVG, HET_diff, F
##
##   Rscript correlation_plot.R --mode froh <diag_file> <het_file> <froh_file> <out_prefix>
##     # legacy v2: correlates F_SNP, F_ROH, D_ROH, D_STD
##
##   Rscript correlation_plot.R --mode froh-cons <diag_file> <het_file> <froh_file> <cons_file> <out_prefix>
##     # legacy v2e: correlates F_SNP, F_ROH, D_ROH, D_STD, ROH_sh (consensus share)
##
##   Rscript correlation_plot.R --mode pairwise <diag_file> <kingkin_file> <eucl_file> <out_prefix>
##     # legacy v3: pairwise-level; correlates G_STD, G_ROH, King_kin, IBS
##
## The output is always two PNGs:
##   <out_prefix>.correlation_heatmap.png
##   <out_prefix>.pairplot.png

require(ggplot2);
require(reshape2);
require(GGally);

args <- commandArgs(trailingOnly = TRUE)

usage <- function(msg = NULL) {
    if (!is.null(msg)) cat("ERROR: ", msg, "\n", sep = "")
    cat("Usage: Rscript correlation_plot.R --mode <basic|froh|froh-cons|pairwise> <inputs...> <out_prefix>\n")
    quit(status = 1)
}

if (length(args) < 4 || args[1] != "--mode") usage("expected --mode as first arg")
mode <- args[2]
rest <- args[-c(1, 2)]

if (length(rest) < 3) usage("not enough positional args for this mode")

out_prefix <- rest[length(rest)]
files <- rest[-length(rest)]
for (f in files) if (!file.exists(f)) stop(paste("Input file not found:", f), call. = FALSE)

## Shared plotting routine -- identical heatmap + pairplot in every legacy variant,
## parameterised by plot size (which varied per legacy script).
render_plots <- function(vars, out_prefix, heatmap_size = 6, pairplot_size = 10,
                         heatmap_dpi = 400, pairplot_dpi = 300) {
    cor_matrix <- cor(vars, use = "pairwise.complete.obs")
    melted_cor <- melt(cor_matrix)

    p_heatmap <- ggplot(melted_cor, aes(Var1, Var2, fill = value)) +
        geom_tile(color = "white") +
        scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0) +
        geom_text(aes(label = sprintf("%.2f", value)), size = 4) +
        theme_minimal(base_size = 14) +
        labs(title = "Correlation Matrix", x = "", y = "", fill = "Correlation")
    ggsave(paste0(out_prefix, ".correlation_heatmap.png"),
           p_heatmap, width = heatmap_size, height = heatmap_size, dpi = heatmap_dpi)

    p_pairs <- ggpairs(vars,
                       upper = list(continuous = "cor"),
                       diag  = list(continuous = "densityDiag"),
                       lower = list(continuous = "smooth"))
    ggsave(paste0(out_prefix, ".pairplot.png"),
           p_pairs, width = pairplot_size, height = pairplot_size, dpi = pairplot_dpi)
}

## ---- Mode dispatch -----------------------------------------------------

create_sorted_key <- function(col1, col2) {
    apply(cbind(as.character(col1), as.character(col2)), 1,
          function(x) paste(sort(x), collapse = "_"))
}

if (mode == "basic") {
    if (length(files) != 2) usage("mode=basic needs <homo_file> <het_file>")
    homo <- read.table(files[1], header = TRUE, comment.char = "", stringsAsFactors = FALSE)
    het  <- read.table(files[2], header = TRUE, comment.char = "", stringsAsFactors = FALSE)
    dat  <- merge(homo, het, by = "IID")
    dat$HET_diff <- dat$O.HET - dat$E.HET
    vars <- dat[, c("KB", "KBAVG", "HET_diff", "F")]
    render_plots(vars, out_prefix, heatmap_size = 6, pairplot_size = 10)

} else if (mode == "froh") {
    if (length(files) != 3) usage("mode=froh needs <diag_file> <het_file> <froh_file>")
    diag <- read.csv(files[1], header = TRUE, comment.char = "", stringsAsFactors = FALSE)
    het  <- read.table(files[2], header = TRUE, comment.char = "", stringsAsFactors = FALSE)
    names(het)[names(het) == "F"] <- "F_SNP"
    froh <- read.table(files[3], header = TRUE, comment.char = "", stringsAsFactors = FALSE)
    dat <- Reduce(function(x, y) merge(x, y, by = "IID", all = TRUE), list(het, froh, diag))
    vars <- dat[, c("F_SNP", "F_ROH", "D_ROH", "D_STD")]
    render_plots(vars, out_prefix, heatmap_size = 6, pairplot_size = 10)

} else if (mode == "froh-cons") {
    if (length(files) != 4) usage("mode=froh-cons needs <diag_file> <het_file> <froh_file> <cons_file>")
    diag <- read.csv(files[1], header = TRUE, comment.char = "", stringsAsFactors = FALSE)
    het  <- read.table(files[2], header = TRUE, comment.char = "", stringsAsFactors = FALSE)
    names(het)[names(het) == "F"] <- "F_SNP"
    froh <- read.table(files[3], header = TRUE, comment.char = "", stringsAsFactors = FALSE)
    cons <- read.table(files[4], header = TRUE, comment.char = "", stringsAsFactors = FALSE)
    names(cons)[names(cons) == "Percent_of_Consensus_ROH"] <- "ROH_sh"
    dat <- Reduce(function(x, y) merge(x, y, by = "IID", all = TRUE), list(het, froh, diag, cons))
    vars <- dat[, c("F_SNP", "F_ROH", "D_ROH", "D_STD", "ROH_sh")]
    render_plots(vars, out_prefix, heatmap_size = 13, pairplot_size = 13)

} else if (mode == "pairwise") {
    if (length(files) != 3) usage("mode=pairwise needs <diag_file> <kingkin_file> <eucl_file>")
    diag <- read.csv(files[1], header = TRUE, comment.char = "", stringsAsFactors = FALSE)
    diag$key <- create_sorted_key(diag$ID1, diag$ID2)
    names(diag)[names(diag) == "Kinship_Std"] <- "G_STD"
    names(diag)[names(diag) == "Kinship_ROH"] <- "G_ROH"

    kingkin_wIBS <- read.table(files[2], header = TRUE, comment.char = "", stringsAsFactors = FALSE)
    kingkin_wIBS$key <- create_sorted_key(kingkin_wIBS$IID1, kingkin_wIBS$IID2)
    names(kingkin_wIBS)[names(kingkin_wIBS) == "KINSHIP"] <- "King_kin"

    euclDist <- read.table(files[3], header = TRUE, comment.char = "", stringsAsFactors = FALSE)
    euclDist$key <- create_sorted_key(euclDist$IID1, euclDist$IID2)
    names(euclDist)[names(euclDist) == "PCA_EUCLIDEAN_DIST"] <- "Eucl_dist"

    dat <- Reduce(function(x, y) merge(x, y, by = "key", all = FALSE),
                  list(diag, kingkin_wIBS, euclDist))
    vars <- dat[, c("G_STD", "G_ROH", "King_kin", "IBS")]
    render_plots(vars, out_prefix, heatmap_size = 8, pairplot_size = 12)

} else {
    usage(paste("unknown mode:", mode))
}
