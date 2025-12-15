# ==============================================================================
# Script Name: 08_Manhattan_QQ_Plots.R
# Description: 
#    This script generates publication-quality Manhattan and QQ plots using CMplot.
#    It visualizes GWAS summary statistics with dual thresholds as defined in 
#    the manuscript (Section 2.4.2).
#
#    Thresholds:
#    - Suggestive:  P < 1e-3 (Blue dashed line)
#    - Significant: P < 1e-4 (Red solid line)
#
# Author:      Liang Xiaotian
# Email:       494382219@qq.com
# Date:        2025/12/06
# License:     MIT
# ==============================================================================

# ------------------------------------------------------------------------------
# Section 1: Environment Setup & Dependency Management
# ------------------------------------------------------------------------------

# Function to check and install missing packages automatically
check_install_packages <- function(pkgs) {
  for (pkg in pkgs) {
    if (!require(pkg, character.only = TRUE)) {
      message(paste0("[Setup] Package '", pkg, "' not found. Installing..."))
      install.packages(pkg, dependencies = TRUE)
      library(pkg, character.only = TRUE)
    }
  }
}

# CMplot is the core package for GWAS visualization
check_install_packages("CMplot")

# Create output directory
output_dir <- file.path("results", "08_Manhattan_QQ_Plots")
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Save original working directory
original_wd <- getwd()
# Safety mechanism: ensure we go back even if script fails
on.exit(setwd(original_wd))

# CMplot saves files to the *current* working directory, so we must switch
setwd(output_dir)

# ------------------------------------------------------------------------------
# Section 2: Data Simulation (Privacy Protection)
# ------------------------------------------------------------------------------

message("\n[Data] Generating synthetic GWAS summary statistics...")

set.seed(2025)

# Define parameters for simulation
n_snps <- 20000
n_chromosomes <- 21 # Oat is hexaploid (6x), 21 chromosomes
trait_names <- c("Trait_A_Env1_BLINK", "Trait_A_Env1_FarmCPU", "Trait_B_Env2_BLINK")

# Function to generate dummy GWAS results
generate_gwas_data <- function(n, n_chr) {
  data.frame(
    SNP = paste0("SNP_", 1:n),
    Chromosome = sort(rep(1:n_chr, length.out = n)),
    Position = unlist(lapply(table(sort(rep(1:n_chr, length.out = n))), function(x) seq(1, x * 10000, length.out = x))),
    # Simulate P-values: mostly null (uniform), with few significant peaks (beta distribution-ish)
    P.value = c(runif(n * 0.99, 0, 1), runif(n * 0.01, 0, 1e-5)) 
  )
}

# Create a list of dataframes to mimic reading multiple GWAS result files
gwas_results_list <- lapply(trait_names, function(name) {
  df <- generate_gwas_data(n_snps, n_chromosomes)
  return(list(name = name, data = df))
})

# ------------------------------------------------------------------------------
# Section 3: Batch Visualization
# ------------------------------------------------------------------------------

# Define significance thresholds (Matched to Manuscript Section 2.4.2)
# "suggestive threshold of P < 1×10-3 and significance threshold of P < 1×10-4"
sig_threshold <- 1e-4
sug_threshold <- 1e-3

message(paste0("\n[Config] Significance Threshold: ", sig_threshold))
message(paste0("[Config] Suggestive Threshold:  ", sug_threshold))
message("\n[Plotting] Starting batch processing...")

for (item in gwas_results_list) {
  
  prefix <- item$name
  data <- item$data
  
  # Standardize columns for CMplot: SNP, Chromosome, Position, P
  plot_data <- data[, c("SNP", "Chromosome", "Position", "P.value")]
  colnames(plot_data) <- c("SNP", "Chromosome", "Position", "P")
  
  message(paste0("  -> Processing: ", prefix))
  
  tryCatch({
    
    # 3.1 Generate QQ Plot
    CMplot(
      plot_data,
      plot.type = "q",
      conf.int = TRUE,
      box = FALSE,
      file = "jpg",
      file.name = prefix,
      dpi = 300,
      width = 5, height = 5,
      verbose = FALSE
    )
    
    # 3.2 Generate Manhattan Plot (Rectangular)
    CMplot(
      plot_data,
      plot.type = "m",
      LOG10 = TRUE,
      ylim = NULL,
      # Dual Thresholds
      threshold = c(sig_threshold, sug_threshold),
      threshold.lty = c(1, 2),        # Solid for Significant, Dashed for Suggestive
      threshold.lwd = c(1, 1),
      threshold.col = c("red", "blue"),
      amplify = FALSE,                # Disable point amplification for clean look
      bin.size = 1e6,
      # Custom Chromosome Colors (Alternating)
      chr.den.col = NULL,
      col = c("grey30", "grey60"),    # Classic publication style
      signal.col = c("red", "blue"),  # Color for significant/suggestive points
      signal.cex = c(1.2, 1.2),
      file = "jpg",
      file.name = prefix,
      dpi = 300,
      width = 14, height = 6,
      verbose = FALSE
    )
    
    # 3.3 Circular Manhattan Plot (Optional but cool for cover images)
    # Only generated for the first trait to save time
    if (item$name == trait_names[1]) {
       CMplot(
         plot_data,
         plot.type = "c",
         r = 0.4,
         cir.legend = TRUE,
         outward = FALSE, 
         cir.legend.col = "black",
         cir.chr.h = 1.3,
         chr.den.col = "black",
         file = "jpg",
         file.name = paste0(prefix, "_Circular"),
         dpi = 300,
         width = 10, height = 10,
         verbose = FALSE
       )
    }
    
    message(paste0("     [Success] Plots saved."))
    
  }, error = function(e) {
    message(paste0("     [Error] Failed: ", e$message))
  })
}

message("\n[Success] File 08 Analysis Complete.")
message(paste0("Plots are available in: ", output_dir))
