# ==============================================================================
# Project: GWAS Pipeline using GAPIT (FarmCPU & BLINK Models)
# Description: 
#   This script demonstrates a standardized workflow for Genome-Wide Association 
#   Studies (GWAS) using the GAPIT package. It performs batch analysis for multiple 
#   traits using FarmCPU and BLINK models.
#
#   Key Features:
#   - Synthetic Genotype (HMP) and Phenotype Data Generation (Privacy Protection)
#   - Automatic Package Installation (Robustness)
#   - Batch Processing with Organized Output Directories (Reproducibility)
#
# Author: Liang Xiaotian
# Email: 494382219@qq.com
# Date: 2025-12-06
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

# GAPIT requires specific installation from GitHub if not present
if (!require("GAPIT")) {
  message("Installing GAPIT from GitHub (jiabowang/GAPIT)...")
  if (!require("devtools")) install.packages("devtools")
  devtools::install_github("jiabowang/GAPIT")
  library(GAPIT)
}

# Load other utilities
install_if_missing(c("data.table", "MASS")) # MASS for multivariate normal distribution

# Create main output directory
main_output_dir <- "GAPIT_GWAS_Output"
if (!dir.exists(main_output_dir)) {
  dir.create(main_output_dir)
}

# Set working directory to output folder to keep project clean
original_wd <- getwd()
setwd(main_output_dir)

# ------------------------------------------------------------------------------
# Section 2: Data Simulation (Privacy Protection)
# NOTE: Replace this section with real data loading:
# myG <- read.table("your_genotype.hmp.txt", head = FALSE)
# myY <- read.table("your_phenotype.txt", head = TRUE, sep="\t")
# ------------------------------------------------------------------------------

message("Generating synthetic Genotype and Phenotype data...")

# 2.1 Simulate Genotype Data (HapMap Format)
# Format: rs# alleles chrom pos strand assembly# center protLSID assayLSID panelLSID QCcode Sample1 ... SampleN
n_samples <- 100
n_snps <- 1000

# Create header info (first 11 columns)
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

# Create Genotype matrix (0, 1, 2 encoding converted to characters for HapMap)
# Simple simulation: random genotypes
geno_matrix <- matrix(
  sample(c("A", "T", "H", "N"), n_samples * n_snps, replace = TRUE, prob = c(0.4, 0.4, 0.15, 0.05)),
  nrow = n_snps,
  ncol = n_samples
)
colnames(geno_matrix) <- paste0("Sample_", 1:n_samples)

# Combine to create myG
myG <- cbind(geno_header, as.data.frame(geno_matrix))

# 2.2 Simulate Phenotype Data
# Format: Taxa Trait1 Trait2 ...
pheno_data <- data.frame(Taxa = colnames(geno_matrix))

# Simulate 3 traits (e.g., Yield, Height, Flowering Time)
traits_list <- c("Trait_A", "Trait_B", "Trait_C")
for (t in traits_list) {
  # Simulate random normal distribution for traits
  pheno_data[[t]] <- rnorm(n_samples, mean = 50, sd = 10)
}

myY <- pheno_data

message("Data generation complete.")
message("Samples: ", n_samples, " | SNPs: ", n_snps, " | Traits: ", length(traits_list))

# ------------------------------------------------------------------------------
# Section 3: Batch GWAS Analysis
# ------------------------------------------------------------------------------

message("Starting GAPIT Batch Analysis...")

# Define analysis parameters
pca_components <- 3   # Number of PCA components to control population structure
models_to_run <- c("FarmCPU", "BLINK") # Efficient models for large datasets

# Loop through each trait
# Note: Trait columns start from index 2 (Index 1 is Taxa)
trait_cols_indices <- 2:ncol(myY) 

for (i in trait_cols_indices) {
  
  # 3.1 Prepare Data for Current Trait
  trait_name <- colnames(myY)[i]
  message("\nProcessing Trait: ", trait_name, " ...")
  
  # Create a dedicated sub-directory for this trait
  if (!dir.exists(trait_name)) {
    dir.create(trait_name)
  }
  
  # Change working directory to the trait folder
  # This ensures all GAPIT output files (Manhattan plots, CSVs) are stored here
  setwd(trait_name)
  
  # Extract current phenotype (Taxa + Trait column)
  current_Y <- myY[, c(1, i)]
  
  # 3.2 Run GAPIT
  # Using tryCatch to prevent the entire loop from stopping if one trait fails
  tryCatch({
    myGAPIT <- GAPIT(
      Y = current_Y,
      G = myG,
      PCA.total = pca_components,
      model = models_to_run,
      Multiple_analysis = TRUE, # Enable multi-model comparison
      file.output = TRUE        # Ensure plots and results are saved
    )
    message("Successfully completed: ", trait_name)
    
  }, error = function(e) {
    message("Error running GAPIT for ", trait_name, ": ", e$message)
  })
  
  # 3.3 Return to the main output directory for the next iteration
  setwd("..")
}

# Restore original working directory
setwd(original_wd)

message("\nGWAS Pipeline Completed Successfully.")
message("Results are saved in: ", file.path(original_wd, main_output_dir))
