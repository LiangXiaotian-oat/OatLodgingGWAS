# ==============================================================================
# Script Name: 07_gapit_gwas.R
# Description: 
#    This script demonstrates a standardized workflow for Genome-Wide Association 
#    Studies (GWAS) using the GAPIT package.
#    
#    Models: FarmCPU & BLINK (Multi-locus GWAS models).
#    Input:  Synthetic HapMap Genotype (HMP) and Phenotype Data.
#    Output: Manhattan plots, QQ plots, and GWAS result tables organized by trait.
#
# Author:      Liang Xiaotian
# Email:       494382219@qq.com
# Date:        2025-12-06
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

# Install GAPIT (Version 3) from GitHub if not present
if (!require("GAPIT")) {
  message("[Setup] Installing GAPIT3 from GitHub...")
  check_install_packages("devtools")
  devtools::install_github("jiabowang/GAPIT3", force = TRUE)
  library(GAPIT)
}

# Load utility packages
check_install_packages(c("data.table", "MASS"))

# Create main output directory
# Standardized output path
main_output_dir <- file.path("results", "07_GWAS_FarmCPU_BLINK")
if (!dir.exists(main_output_dir)) {
  dir.create(main_output_dir, recursive = TRUE)
}

# Save original working directory to restore later
original_wd <- getwd()
# Safety mechanism: ensure we go back even if script fails
on.exit(setwd(original_wd))

# Switch to output dir
setwd(main_output_dir)

# ------------------------------------------------------------------------------
# Section 2: Data Simulation (Privacy Protection)
# ------------------------------------------------------------------------------

message("\n[Data] Generating synthetic HapMap Genotype and Phenotype data...")

# Set seed for reproducibility (Consistent with project)
set.seed(2025)

# 2.1 Simulate Genotype Data (HapMap Format)
n_samples <- 100
n_snps <- 1000

# Create header info (first 11 columns required by HapMap format)
# Format: rs# alleles chrom pos strand assembly# center protLSID assayLSID panelLSID QCcode
snp_ids <- paste0("SNP_", 1:n_snps)
chroms <- sort(rep(1:5, length.out = n_snps))
pos <- unlist(lapply(table(chroms), function(x) seq(1, x * 1000, length.out = x)))

geno_header <- data.frame(
  rs = snp_ids,
  alleles = "A/T",
  chrom = chroms,
  pos = pos,
  strand = "+",
  assembly = "NA",
  center = "NA",
  protLSID = "NA",
  assayLSID = "NA",
  panelLSID = "NA",
  QCcode = "NA"
)

# Create Genotype matrix
# GAPIT HapMap accepts characters like "A", "T", "G", "C" or "A/T" etc.
# We simulate random calls
geno_matrix <- matrix(
  sample(c("A", "T", "C", "G"), n_samples * n_snps, replace = TRUE),
  nrow = n_snps,
  ncol = n_samples
)
colnames(geno_matrix) <- paste0("Sample_", 1:n_samples)

# Combine to create myG (The G parameter for GAPIT)
myG <- cbind(geno_header, as.data.frame(geno_matrix))

# 2.2 Simulate Phenotype Data
# Format: Taxa Trait1 Trait2 ...
pheno_data <- data.frame(Taxa = colnames(geno_matrix))

# Simulate 3 traits (e.g., Yield, Height, Flowering Time)
traits_list <- c("Trait_A", "Trait_B", "Trait_C")
for (t in traits_list) {
  # Simulate random normal distribution
  pheno_data[[t]] <- rnorm(n_samples, mean = 50, sd = 10)
}

myY <- pheno_data

message(paste0("[Data] Generation complete. Samples: ", n_samples, " | SNPs: ", n_snps))

# ------------------------------------------------------------------------------
# Section 3: Batch GWAS Analysis
# ------------------------------------------------------------------------------

message("\n[Analysis] Starting GAPIT Batch Analysis (FarmCPU & BLINK)...")

# Define analysis parameters
pca_components <- 3 
models_to_run <- c("FarmCPU", "BLINK")

# Loop through each trait
# Note: Trait columns start from index 2
trait_cols_indices <- 2:ncol(myY) 

for (i in trait_cols_indices) {
  
  trait_name <- colnames(myY)[i]
  message(paste0("\n  -> Processing Trait: ", trait_name))
  
  # Create a dedicated sub-directory for this trait to avoid file clutter
  if (!dir.exists(trait_name)) {
    dir.create(trait_name)
  }
  
  # Change working directory to the trait folder
  setwd(trait_name)
  
  # Extract current phenotype (Taxa + Trait column)
  current_Y <- myY[, c(1, i)]
  
  # Run GAPIT
  tryCatch({
    myGAPIT <- GAPIT(
      Y = current_Y,
      G = myG,                # Using HapMap format directly
      PCA.total = pca_components,
      model = models_to_run,
      Multiple_analysis = TRUE,
      file.output = TRUE
    )
    message(paste0("  [Success] ", trait_name, " finished."))
    
  }, error = function(e) {
    message(paste0("  [Error] Failed for ", trait_name, ": ", e$message))
  })
  
  # Return to the main output directory for the next iteration
  setwd("..")
}

message("\n[Success] File 07 Analysis Complete.")
message(paste0("Results are saved in: ", main_output_dir))
