# ==============================================================================
# Script Name: 03_correlation_visualization_GGally.R
# Description: 
#    This script performs multi-environment correlation matrix visualization.
#    It uses 'GGally::ggpairs' to display pairwise correlations, scatter plots,
#    and density distributions, colored by Environment (Env).
#    This visualizes the stability of trait relationships across environments.
#
#    Includes synthetic data generation for privacy protection and reproducibility.
#
# Author:      Liang Xiaotian
# Email:       494382219@qq.com
# Date:        2025-12-06
# License:     MIT
# ==============================================================================

# ------------------------------------------------------------------------------
# Section 1: Environment Setup & Dependency Management
# ------------------------------------------------------------------------------

# List of required packages
required_packages <- c("ggplot2", "GGally")

# Function to check and install packages automatically (Consistent with previous files)
check_install_packages <- function(pkgs) {
  for (pkg in pkgs) {
    if (!require(pkg, character.only = TRUE)) {
      message(paste0("[Setup] Package '", pkg, "' not found. Installing..."))
      install.packages(pkg, dependencies = TRUE)
      library(pkg, character.only = TRUE)
    } else {
      message(paste0("[Setup] Package '", pkg, "' is loaded."))
    }
  }
}

# Load packages
check_install_packages(required_packages)

# Create output directory
output_dir <- "results"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# ------------------------------------------------------------------------------
# Section 2: Data Simulation (Privacy Protection)
# ------------------------------------------------------------------------------

message("\n[Data] Generating synthetic multi-environment data...")

# Set seed for reproducibility (Consistent with previous files)
set.seed(2025)

# 2.1 Define Simulation Parameters
n_samples <- 300
# Using 'Env' to represent Environments (Wenjiang, Chongzhou, etc.)
groups <- c("Env_A", "Env_B", "Env_C", "Env_D", "Env_E")

# 2.2 Generate Random Data
data <- data.frame(
  ID = 1:n_samples,
  Env = sample(groups, n_samples, replace = TRUE)
)

# Simulate correlated traits (Using generic names consistent with File 02)
# Trait_A serves as a hub trait (e.g., Plant Height)
data$Trait_A <- rnorm(n_samples, mean = 50, sd = 10)
# Trait_B: Positive correlation with A
data$Trait_B <- data$Trait_A * 0.8 + rnorm(n_samples, mean = 0, sd = 5) 
# Trait_C: Independent
data$Trait_C <- rnorm(n_samples, mean = 20, sd = 4)
# Trait_D: Negative correlation with C
data$Trait_D <- -data$Trait_C * 0.6 + rnorm(n_samples, mean = 100, sd = 15) 
# Trait_E: Weak interaction
data$Trait_E <- data$Trait_A * 0.3 + runif(n_samples, 0, 10)

# Convert grouping variable to factor
data$Env <- factor(data$Env)

# Display data structure
print(str(data))

# ------------------------------------------------------------------------------
# Section 3: Correlation Matrix Visualization
# ------------------------------------------------------------------------------

message("\n[Plotting] Generating GGally correlation matrix...")

# 3.1 Define Plot Parameters
# Define custom colors for the 5 environments (Publication-ready palette)
custom_colors <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd")

# Select columns for analysis (Exclude ID and Env columns)
# Dynamically select all numeric columns or specific trait names
plot_cols <- c("Trait_A", "Trait_B", "Trait_C", "Trait_D", "Trait_E")

# 3.2 Create ggpairs Plot
p <- ggpairs(
  data = data,
  columns = plot_cols,                # Columns to include
  aes(color = Env, alpha = 0.6),      # Group by Environment
  
  # Upper triangle: Correlation coefficients
  upper = list(continuous = wrap("cor", size = 3, alignPercent = 0.8)),
  
  # Lower triangle: Scatter plots
  lower = list(continuous = wrap("points", alpha = 0.6, size = 1.2)),
  
  # Diagonal: Density plots
  diag = list(continuous = wrap("densityDiag", alpha = 0.5)),
  
  title = "Multi-Environment Trait Correlations"
) +
  # Apply custom colors
  scale_color_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors) +
  
  # Customize theme
  theme_bw() +
  theme(
    axis.text = element_text(color = "black", size = 8),
    axis.title = element_text(size = 10, face = "bold"),
    strip.text = element_text(face = "bold", size = 10),
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold", margin = margin(b = 15)),
    legend.position = "bottom",
    legend.title = element_text(face = "bold")
  )

# ------------------------------------------------------------------------------
# Section 4: Save Results
# ------------------------------------------------------------------------------

output_filename <- file.path(output_dir, "03_GGally_Correlation_Matrix.png")

# Save as high-resolution PNG
ggsave(
  filename = output_filename,
  plot = p,
  width = 12, 
  height = 12,
  dpi = 300,
  bg = "white"
)

message(paste("[Success] Plot saved to:", output_filename))
message("File 03 Analysis Complete.")
