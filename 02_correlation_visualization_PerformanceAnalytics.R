# ==============================================================================
# Script Name: 02_correlation_visualization_PerformanceAnalytics.R
# Description: 
#    This script performs comprehensive correlation analysis.
#    1. Calculates Pearson correlation coefficients (r) and p-values using Hmisc.
#    2. Exports a flattened summary table of correlations.
#    3. Visualizes correlations using 'corrplot' (heatmap) and 'PerformanceAnalytics'.
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
required_packages <- c("Hmisc", "corrplot", "PerformanceAnalytics")

# Function to check and install packages automatically (Consistent with File 01)
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
output_dir <- "results" # Consistent with File 01
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# ------------------------------------------------------------------------------
# Section 2: Data Preparation (Synthetic Data)
# ------------------------------------------------------------------------------

message("\n[Data] Generating synthetic data for correlation analysis...")

# Set seed for reproducibility (Consistent with File 01)
set.seed(2025)

# Simulate a dataset with 100 samples and 6 variables (Traits)
n_samples <- 100
mydata <- data.frame(
  SampleID = paste0("Sample_", 1:n_samples), # Added ID column to test robustness
  Trait_A = rnorm(n_samples, mean = 50, sd = 10),
  Trait_B = rnorm(n_samples, mean = 20, sd = 5),
  Trait_C = runif(n_samples, 0, 100),
  Trait_D = rnorm(n_samples, mean = 10, sd = 2),
  Trait_E = rnorm(n_samples, mean = 0, sd = 1),
  Trait_F = runif(n_samples, 10, 20)
)

# Introduce some correlations artificially
# Trait_B is positively correlated with Trait_A
mydata$Trait_B <- mydata$Trait_A * 0.7 + rnorm(n_samples, 0, 5)
# Trait_E is negatively correlated with Trait_D
mydata$Trait_E <- -mydata$Trait_D * 0.5 + rnorm(n_samples, 0, 1)

# Display first few rows
print(head(mydata))

# ------------------------------------------------------------------------------
# Section 3: Correlation Calculation (Hmisc)
# ------------------------------------------------------------------------------

message("\n[Analysis] Calculating correlation matrix and p-values...")

# IMPORTANT: Select ONLY numeric columns for correlation analysis
# This prevents errors if your data has 'SampleID' or 'Env' columns
numeric_data <- Filter(is.numeric, mydata)

# Use rcorr from Hmisc to get both r and p values
res <- rcorr(as.matrix(numeric_data), type = "pearson") 

# Function to flatten the correlation matrix into a table format
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor = (cormat)[ut],
    p = pmat[ut]
  )
}

# Apply the function
flat_corr_data <- flattenCorrMatrix(res$r, res$P)

# View top 10 significant correlations
message("Top 10 significant correlations:")
print(head(flat_corr_data[order(flat_corr_data$p), ], 10))

# Save flattened matrix to CSV
write.csv(flat_corr_data, file = file.path(output_dir, "02_Correlation_Summary_Table.csv"), row.names = FALSE)

# ------------------------------------------------------------------------------
# Section 4: Visualization - Method 1: corrplot
# ------------------------------------------------------------------------------

message("\n[Plotting] Generating visualization: corrplot...")

plot1_path <- file.path(output_dir, "02_Corrplot_Visualization.png")

# Save plot to PNG
png(filename = plot1_path, width = 2000, height = 2000, res = 300)

# Plot 1: Correlation matrix with significance levels
corrplot(res$r, 
         type = "upper", 
         order = "hclust", 
         p.mat = res$P, 
         sig.level = 0.01, 
         insig = "blank", 
         tl.col = "black", 
         tl.srt = 45,
         title = "Correlation Matrix (p < 0.01)",
         mar = c(0,0,2,0))

dev.off()
message(paste("Saved plot to:", plot1_path))

# ------------------------------------------------------------------------------
# Section 5: Visualization - Method 2: PerformanceAnalytics
# ------------------------------------------------------------------------------

message("\n[Plotting] Generating visualization: PerformanceAnalytics chart...")

plot2_path <- file.path(output_dir, "02_PerformanceAnalytics_Chart.png")

# Save plot to PNG
png(filename = plot2_path, width = 2400, height = 2400, res = 300)

# Plot 2: Comprehensive correlation chart
chart.Correlation(numeric_data, histogram = TRUE, pch = 19)

dev.off()
message(paste("Saved plot to:", plot2_path))

message("\n[Success] File 02 Analysis Complete.")
