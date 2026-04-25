args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Usage: Rscript pca_plots.R <prefix> <eigenvec_suffix> <color_col> <output.png> [n_pcs=3] [color_type=factor]",
       call. = FALSE)
}

require(ggplot2)
require(gridExtra)

prefix         <- args[1]
eigenvec_suffix <- args[2]
color_col      <- args[3]
out_png        <- args[4]
n_pcs          <- if (length(args) >= 5) as.integer(args[5]) else 3L
color_type     <- if (length(args) >= 6) args[6] else "factor"

if (!n_pcs %in% c(3L, 6L)) stop("n_pcs must be 3 or 6", call. = FALSE)
if (!color_type %in% c("factor", "numeric")) stop("color_type must be 'factor' or 'numeric'", call. = FALSE)

eigenvec_file <- paste(prefix, "eigenvec", eigenvec_suffix, sep = ".")
eigenval_file <- paste(prefix, "eigenval", sep = ".")

if (!file.exists(eigenvec_file)) stop(paste("eigenvec file not found:", eigenvec_file), call. = FALSE)
if (!file.exists(eigenval_file)) stop(paste("eigenval file not found:", eigenval_file), call. = FALSE)

eigenvec <- read.table(eigenvec_file, header = TRUE, comment.char = "")
eigenval <- scan(eigenval_file, quiet = TRUE)

var_explained <- round(100 * eigenval / sum(eigenval), 2)
pc_labs <- paste0("PC", seq_along(var_explained), " (", var_explained, "%)")

if (color_type == "factor") {
  eigenvec[[color_col]] <- as.factor(eigenvec[[color_col]])
  color_scale <- NULL
} else {
  require(viridis)
  eigenvec[[color_col]] <- as.numeric(eigenvec[[color_col]])
  color_scale <- scale_color_viridis_c(option = "plasma")
}

make_plot <- function(xpc, ypc) {
  p <- ggplot(eigenvec, aes(x = .data[[paste0("PC", xpc)]],
                            y = .data[[paste0("PC", ypc)]],
                            col = .data[[color_col]])) +
    geom_point() +
    labs(title = "PCA Plot",
         x = pc_labs[xpc], y = pc_labs[ypc],
         color = if (color_type == "numeric") color_col else waiver())
  if (!is.null(color_scale)) p <- p + color_scale
  p
}

if (n_pcs == 3L) {
  plots <- list(make_plot(1, 2), make_plot(1, 3), make_plot(2, 3))
  combined <- arrangeGrob(grobs = plots, nrow = 1)
  ggsave(file = out_png, combined, width = 16, height = 4, dpi = 400)
} else {
  plots <- list(make_plot(1, 2), make_plot(1, 3), make_plot(1, 4),
                make_plot(2, 3), make_plot(2, 4), make_plot(3, 4))
  combined <- arrangeGrob(grobs = plots, nrow = 2)
  ggsave(file = out_png, combined, width = 16, height = 8, dpi = 400)
}
