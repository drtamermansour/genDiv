#!/usr/bin/env Rscript

# Load necessary library for plotting
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2", repos = "http://cran.us.r-project.org")
  library(ggplot2)
} else {
  library(ggplot2)
}

# 1. Receive command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# UPDATED: Check for 5 arguments: input, x, y, color factor, and shape factor
if (length(args) != 5) {
  stop("Usage: Rscript plot_correlation.R <input_file.tsv> <column_name_x> <column_name_y> <color_factor_col> <shape_factor_col>\n", call. = FALSE)
}

input_file <- args[1]
col_x_name <- args[2]
col_y_name <- args[3]
col_color_name <- args[4]
col_shape_name <- args[5]

# 2. Read the tab-separated file
if (!file.exists(input_file)) {
  stop("Input file not found: ", input_file, call. = FALSE)
}

data <- read.table(input_file, sep = "\t", comment.char = "", header = TRUE)

# 3. Validate column names
if (!(col_x_name %in% names(data) && col_y_name %in% names(data) && 
      col_color_name %in% names(data) && col_shape_name %in% names(data))) {
  stop("One or more column names not found in the data file.\n", call. = FALSE)
}

# Extract the columns and ensure they are numeric/factor
x_var <- as.numeric(data[[col_x_name]])
y_var <- as.numeric(data[[col_y_name]])
color_var <- as.factor(data[[col_color_name]])
shape_var <- as.factor(data[[col_shape_name]])

# Remove NA values for correlation calculation and plotting
# Use complete.cases to ensure rows are aligned across all variables
keep_idx <- complete.cases(x_var, y_var, color_var, shape_var)
cleaned_data <- data.frame(
  x = x_var[keep_idx], 
  y = y_var[keep_idx], 
  color = color_var[keep_idx],
  shape = shape_var[keep_idx]
)

# 4. Calculate the correlation (Global)
correlation_val <- cor(cleaned_data$x, cleaned_data$y, method = "pearson")
correlation_text <- paste("Correlation (Pearson):", round(correlation_val, 3))

# 5. Define Output Path (Saved to the directory of the input file)
input_dir <- dirname(input_file)
base_name <- paste0("correlation_plot_", col_x_name, "_vs_", col_y_name, "_doubleAnn" , ".png")
output_path <- file.path(input_dir, base_name)

# 6. Plot the correlation using ggplot2
plot_title <- paste("Correlation Plot of", col_x_name, "vs", col_y_name)

Sys.setenv(TZ = "America/Edmonton") 
ggplot(cleaned_data, aes(x = x, y = y)) +
  geom_point(aes(color = color, shape = shape), alpha = 0.7, size = 2) + # Added shape mapping
  geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed") + # Global trend line
  labs(title = plot_title,
       x = col_x_name,
       y = col_y_name,
       color = col_color_name,
       shape = col_shape_name,
       subtitle = correlation_text) +
  theme_minimal()

# 7. Save the plot to the same directory as the input file
ggsave(output_path)

cat("Plot saved as:", output_path, "\n")
cat(correlation_text, "\n")
