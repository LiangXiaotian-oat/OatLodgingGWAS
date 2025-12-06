# ==============================================================================
# Project: GWAS Visualization Pipeline (Manhattan & QQ Plots)
# Description: 
#   This script demonstrates how to generate publication-quality Manhattan and 
#   QQ plots from GWAS summary statistics using the CMplot package. It includes 
#   data simulation to ensure privacy protection and reproducibility.
#
# Author: Liang Xiaotian
# Email: 494382219@qq.com
# Date: 2025/12/06
# License: MIT
# ==============================================================================

# ------------------------------------------------------------------------------
# Section 1: Environment Setup & Dependency Management
# ------------------------------------------------------------------------------

# Function to check and install missing packages automatically
install_if_missing <- function(packages) {
  new_packages <- packages[!(packages %in% installed.packages()[, "Package"])]
  if (length(new_packages)) {
    message("Installing missing packages: ", paste(new_packages, collapse = ", "))
    install.packages(new_packages)
  }
  invisible(sapply(packages, library, character.only = TRUE))
}

# Load required packages
# CMplot is the core package for GWAS visualization
install_if_missing(c("CMplot"))

# Create output directory to keep workspace clean
output_dir <- "GWAS_Plots_Output"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# Set working directory to output folder for CMplot (it saves files to current wd)
original_wd <- getwd()
setwd(output_dir)

# ------------------------------------------------------------------------------
# Section 2: Data Simulation (Privacy Protection)
# NOTE: Replace this section with real data loading in actual analysis.
# Example: file_list <- list.files(pattern = "*.csv")
# ------------------------------------------------------------------------------

message("Generating synthetic GWAS summary statistics...")

# Define parameters for simulation
n_snps <- 50000
n_chromosomes <- 7
trait_names <- c("Trait_A_Env1_BLINK", "Trait_A_Env1_FarmCPU", "Trait_B_Env2_BLINK")

# Function to generate dummy GWAS results
generate_gwas_data <- function(n, n_chr) {
  data.frame(
    SNP = paste0("SNP_", 1:n),
    Chromosome = sort(rep(1:n_chr, length.out = n)),
    Position = unlist(lapply(table(sort(rep(1:n_chr, length.out = n))), function(x) seq(1, x * 1000, length.out = x))),
    # Simulate P-values: mostly uniform distribution (null), with few significant peaks
    P.value = c(runif(n * 0.99, 0, 1), runif(n * 0.01, 0, 1e-5)) 
  )
}

# Create a list of dataframes to mimic reading multiple CSV files
gwas_results_list <- lapply(trait_names, function(name) {
  df <- generate_gwas_data(n_snps, n_chromosomes)
  return(list(name = name, data = df))
})

# ------------------------------------------------------------------------------
# Section 3: Batch Visualization (Manhattan & QQ Plots)
# ------------------------------------------------------------------------------

# Define significance thresholds
# Stricter threshold (P=0.0001, corresponding to -log10(P)=4.0)
stricter_threshold <- 0.0001 
# Relaxed threshold: arbitrary cutoff (e.g., -log10(P) = 3)
relaxed_threshold <- 0.001

message("Starting batch processing...")
message("Bonferroni Threshold: ", bonferroni_threshold)
message("Relaxed Threshold: ", relaxed_threshold)

for (item in gwas_results_list) {
  tryCatch({
    # 3.1 Data Preparation
    prefix <- item$name
    data <- item$data
    
    # Rename columns to standard format required by CMplot (SNP, Chromosome, Position, P-value)
    # Note: Adjust column indices or names based on your actual CSV format
    plot_data <- data[, c("SNP", "Chromosome", "Position", "P.value")]
    colnames(plot_data) <- c("SNP", "Chromosome", "Position", "P")
    
    cat("\nProcessing: ", prefix, "\n")
    
    # 3.2 Generate QQ Plot (PNG & PDF)
    # QQ plot assesses if the statistical model controls for population structure
    for (fmt in c("jpg", "pdf")) {
      CMplot(
        plot_data,
        plot.type = "q",
        box = FALSE,
        file = fmt,
        file.name = paste0(prefix, "_QQ"),
        dpi = 300,
        conf.int = TRUE,
        conf.int.col = NULL,
        threshold.col = "red",
        threshold.lty = 2,
        file.output = TRUE,
        verbose = FALSE,
        width = 5,
        height = 5
      )
    }
    
    # 3.3 Generate Manhattan Plot (PNG & PDF)
    # Manhattan plot visualizes significant SNPs across chromosomes
    for (fmt in c("jpg", "pdf")) {
      CMplot(
        plot_data,
        plot.type = "m",
        LOG10 = TRUE, # Automatically convert P-values to -log10(P)
        ylim = NULL,
        # Set multiple threshold lines
        threshold = c(bonferroni_threshold, relaxed_threshold),
        threshold.lty = c(1, 2),      # Solid line for Bonferroni, Dashed for Relaxed
        threshold.lwd = c(1, 1),
        threshold.col = c("black", "grey"),
        amplify = FALSE,              # Disable highlighting specific points to keep it clean
        bin.size = 1e6,               # Density bin size
        chr.den.col = c("darkgreen", "yellow", "red"), # Chromosome density color
        signal.col = c("red", "green"), # Color for significant points above thresholds
        signal.cex = c(1, 1),
        signal.pch = c(19, 19),
        file = fmt,
        dpi = 300,
        file.output = TRUE,
        verbose = FALSE,
        file.name = paste0(prefix, "_Manhattan"),
        width = 12,  # Wider for Manhattan plots
        height = 6
      )
    }
    
    cat("Successfully generated plots for: ", prefix, "\n")
    
  }, error = function(e) {
    cat("Error processing ", prefix, ": ", e$message, "\n")
  })
}

# Restore original working directory
setwd(original_wd)

message("\nAll tasks completed. Check the 'GWAS_Plots_Output' folder for results.")
