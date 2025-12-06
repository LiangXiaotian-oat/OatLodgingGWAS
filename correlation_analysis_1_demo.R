# ==============================================================================
# Project: Comprehensive Correlation Analysis & Visualization Pipeline
# Description: 
#   This script demonstrates how to perform correlation analysis using the Hmisc package.
#   It calculates Pearson/Spearman correlation coefficients and p-values, flattens the 
#   correlation matrix into a readable table, and visualizes the results using 
#   corrplot and PerformanceAnalytics.
#
#   Includes synthetic data generation for privacy protection and reproducibility.
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
required_packages <- c("Hmisc", "corrplot", "PerformanceAnalytics")

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

# Create output directory
output_dir <- "Correlation_Analysis_Output"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# ------------------------------------------------------------------------------
# Section 2: Data Preparation (Synthetic Data)
# NOTE: In real analysis, replace this section with: mydata <- read.csv("your_data.csv")
# ------------------------------------------------------------------------------

message("Generating synthetic data for correlation analysis...")

# Set seed for reproducibility
set.seed(123)

# Simulate a dataset with 100 samples and 6 variables (Traits)
n_samples <- 100
mydata <- data.frame(
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

message("Calculating correlation matrix and p-values...")

# Use rcorr from Hmisc to get both r and p values
# Input must be a matrix
res <- rcorr(as.matrix(mydata), type = "pearson") # Options: "pearson", "spearman"

# Extract correlation coefficients (r) and p-values (P)
cor_matrix <- res$r
p_matrix <- res$P

# Function to flatten the correlation matrix into a table format
# Arguments:
#   cormat: matrix of the correlation coefficients
#   pmat: matrix of the correlation p-values
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
print("Top significant correlations:")
print(head(flat_corr_data[order(flat_corr_data$p), ], 10))

# Save flattened matrix to CSV
write.csv(flat_corr_data, file = paste0(output_dir, "/Correlation_Summary_Table.csv"), row.names = FALSE)

# ------------------------------------------------------------------------------
# Section 4: Visualization - Method 1: corrplot
# ------------------------------------------------------------------------------

message("Generating visualization: corrplot...")

# Define output filename
plot1_path <- paste0(output_dir, "/Corrplot_Visualization.png")

# Save plot to PNG
png(filename = plot1_path, width = 2000, height = 2000, res = 300)

# Plot 1: Correlation matrix with significance levels
# Insignificant correlations (p > 0.01) are left blank
corrplot(res$r, 
         type = "upper", 
         order = "hclust", 
         p.mat = res$P, 
         sig.level = 0.01, 
         insig = "blank", 
         tl.col = "black", 
         tl.srt = 45,
         title = "Correlation Matrix (Insignificant Hidden)",
         mar = c(0,0,2,0)) # Adjust margin for title

dev.off()
message("Saved plot to: ", plot1_path)

# ------------------------------------------------------------------------------
# Section 5: Visualization - Method 2: PerformanceAnalytics
# ------------------------------------------------------------------------------

message("Generating visualization: PerformanceAnalytics chart...")

# Define output filename
plot2_path <- paste0(output_dir, "/PerformanceAnalytics_Chart.png")

# Save plot to PNG
png(filename = plot2_path, width = 2400, height = 2400, res = 300)

# Plot 2: Comprehensive correlation chart
# Includes histograms (diagonal), scatter plots (lower), and correlation values (upper)
chart.Correlation(mydata, histogram = TRUE, pch = 19)

dev.off()
message("Saved plot to: ", plot2_path)

message("Analysis Complete. Check the output folder.")
