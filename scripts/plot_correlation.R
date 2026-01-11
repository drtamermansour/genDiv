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

# Check if the correct number of arguments is provided
if (length(args) != 3) {
  stop("Usage: Rscript plot_correlation.R <input_file.tsv> <column_name_x> <column_name_y>\n", call. = FALSE)
}

input_file <- args[1]
col_x_name <- args[2]
col_y_name <- args[3]

# 2. Read the tab-separated file
if (!file.exists(input_file)) {
  stop("Input file not found: ", input_file, call. = FALSE)
}

data <- read.table(input_file, sep = "\t", comment.char = "", header = TRUE)
# print the header of the data for debugging
print(head(data))

# 3. Validate column names
if (!(col_x_name %in% names(data) && col_y_name %in% names(data))) {
  stop("One or both column names not found in the data file.\n", call. = FALSE)
}

# Extract the columns and ensure they are numeric
x_var <- as.numeric(data[[col_x_name]])
y_var <- as.numeric(data[[col_y_name]])

if (any(is.na(x_var)) || any(is.na(y_var))) {
    warning("NA values introduced during conversion to numeric. Non-numeric data may be present.")
}

# Remove NA values for correlation calculation and plotting
cleaned_data <- data.frame(x = na.omit(x_var), y = na.omit(y_var))

# 4. Calculate the correlation
correlation_val <- cor(cleaned_data$x, cleaned_data$y, method = "pearson")
correlation_text <- paste("Correlation (Pearson):", round(correlation_val, 3))

# 5. Plot the correlation using ggplot2
plot_title <- paste("Correlation Plot of", col_x_name, "vs", col_y_name)
output_filename <- paste0("correlation_plot_", col_x_name, "_vs_", col_y_name, ".png")

ggplot(cleaned_data, aes(x = x, y = y)) +
  geom_point() + # Scatter plot points
  geom_smooth(method = "lm", se = FALSE, color = "blue") + # Add a linear regression line
  labs(title = plot_title,
       x = col_x_name,
       y = col_y_name,
       subtitle = correlation_text) +
  theme_minimal()

# 6. Save the plot to a file
ggsave(output_filename)

cat("Plot saved as:", output_filename, "\n")
cat(correlation_text, "\n")

