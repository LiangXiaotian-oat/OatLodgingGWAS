# ==============================================================================
# Project: Multi-Trait Correlation Analysis Pipeline
# Description: 
#   This script performs correlation matrix visualization for multiple phenotypic 
#   traits across different environments/targets using the GGally package.
#   It includes synthetic data generation to ensure privacy and reproducibility.
#
# Author: Liang Xiaotian
# Email: 494382219@qq.com
# Date: 2025-12-06
# License: MIT
# ==============================================================================

# ------------------------------------------------------------------------------
# Section 1: Environment Setup & Dependency Management
# ------------------------------------------------------------------------------

# List of required packages
required_packages <- c("ggplot2", "GGally")

# Function to check and install missing packages automatically
install_if_missing <- function(packages) {
  new_packages <- packages[!(packages %in% installed.packages()[, "Package"])]
  if (length(new_packages)) {
    message("Installing missing packages: ", paste(new_packages, collapse = ", "))
    install.packages(new_packages)
  }
  invisible(sapply(packages, library, character.only = TRUE))
}

# Load packages
install_if_missing(required_packages)

# Create output directory to keep workspace clean
output_dir <- "Correlation_Plots_Output"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# ------------------------------------------------------------------------------
# Section 2: Data Simulation (Privacy Protection)
# NOTE: In real analysis, replace this section with read.csv("your_data.csv")
# ------------------------------------------------------------------------------

message("Generating synthetic correlation data...")

# Set seed for reproducibility
set.seed(123)

# 2.1 Define Simulation Parameters
n_samples <- 300
groups <- c("Env_A", "Env_B", "Env_C", "Env_D", "Env_E")
# Define traits (Generalized names)
trait_names <- c("Trait_1", "Trait_2", "Trait_3", "Trait_4", "Trait_5", "Trait_6")

# 2.2 Generate Random Data with Correlation Structure
# Create base dataframe
data <- data.frame(
  ID = 1:n_samples,
  Target = sample(groups, n_samples, replace = TRUE)
)

# Simulate correlated traits (e.g., Trait_1 is correlated with Trait_2)
data$Trait_1 <- rnorm(n_samples, mean = 50, sd = 10)
data$Trait_2 <- data$Trait_1 * 0.8 + rnorm(n_samples, mean = 0, sd = 5) # Positive correlation
data$Trait_3 <- rnorm(n_samples, mean = 20, sd = 4)
data$Trait_4 <- -data$Trait_3 * 0.6 + rnorm(n_samples, mean = 100, sd = 15) # Negative correlation
data$Trait_5 <- runif(n_samples, 0, 10)
data$Trait_6 <- data$Trait_1 * 0.3 + data$Trait_5 * 0.5 + rnorm(n_samples) # Complex relationship

# Convert grouping variable to factor
data$Target <- factor(data$Target)

# ------------------------------------------------------------------------------
# Section 3: Correlation Matrix Visualization
# ------------------------------------------------------------------------------

message("Generating correlation matrix plot...")

# 3.1 Define Plot Parameters
# Define custom colors for different groups
custom_colors <- c("steelblue", "yellowgreen", "#e31a1c", "purple", "#ff7f00")

# Select columns for analysis (excluding ID and Target)
analysis_cols <- trait_names 

# 3.2 Create ggpairs Plot
# Using 'Target' as the grouping variable for color mapping
p <- ggpairs(
  data = data,
  columns = analysis_cols,       # Columns to include in the matrix
  aes(color = Target, alpha = 0.6), # Grouping by 'Target'
  
  # Upper triangle: Correlation coefficients
  upper = list(continuous = wrap("cor", size = 3, alignPercent = 0.8)),
  
  # Lower triangle: Scatter plots
  lower = list(continuous = wrap("points", alpha = 0.6, size = 1.2)),
  
  # Diagonal: Density plots
  diag = list(continuous = wrap("densityDiag", alpha = 0.5)),
  
  title = "Multi-Trait Correlation Analysis Matrix"
) +
  # Apply custom colors
  scale_color_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors) +
  
  # Customize theme for publication quality
  theme_bw() +
  theme(
    axis.text = element_text(colour = "black", size = 8),
    axis.title = element_text(size = 10, face = "bold"),
    strip.text = element_text(face = "bold", size = 10), # Facet label size
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold", margin = margin(b = 15)),
    legend.position = "bottom",
    legend.title = element_text(face = "bold")
  )

# ------------------------------------------------------------------------------
# Section 4: Save Results
# ------------------------------------------------------------------------------

# Define output filename
output_filename <- paste0(output_dir, "/Correlation_Matrix_Plot.png")

# Save as high-resolution PNG
ggsave(
  filename = output_filename,
  plot = p,
  width = 14,  # Adjusted for 6x6 matrix visibility
  height = 12,
  dpi = 300,
  bg = "white" # Ensure white background
)

message("Plot successfully saved to: ", output_filename)
message("Analysis Complete.")
