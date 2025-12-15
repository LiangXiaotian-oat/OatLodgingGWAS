# -------------------------------------------------------------------------
# Script Name: 01_phenotype_stats.R
# Description: Generates synthetic phenotype data (multi-environment) and 
#              calculates descriptive statistics (Mean, SD, CV, etc.).
#              Demonstrates data preprocessing (averaging replicates) and
#              statistical summary without requiring raw datasets.
#
# Author:      Liang Xiaotian
# Email:       494382219@qq.com
# Date:        2025/12/6
# License:     MIT
# -------------------------------------------------------------------------

# =========================================================================
# 1. Environment Setup & Dependency Management
# =========================================================================

# List of required packages
required_packages <- c("tidyverse", "moments")

# Function to check and install packages automatically
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

# Run setup
check_install_packages(required_packages)

# =========================================================================
# 2. Configuration & Parameters
# =========================================================================

# Simulation Parameters
set.seed(2025)       # Ensure reproducibility
n_samples <- 386     # Number of accessions
n_envs <- 5          # Number of environments

# Define Environment Names (Generic)
env_names <- paste0("Environment_", 1:n_envs)

# Define Traits (Generic)
# Traits with replicates (e.g., Trait_A_1, Trait_A_2, Trait_A_3)
traits_with_reps <- c("Trait_A", "Trait_B", "Trait_C", "Trait_D", "Trait_E")
# Single value traits (e.g., Trait_F)
traits_single <- c("Trait_F")

# Output Directory
output_dir <- "results"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  message(paste0("[Setup] Created output directory: ", output_dir))
}

# =========================================================================
# 3. Module: Synthetic Data Generation
# =========================================================================

generate_synthetic_data <- function(env_name, n_samples) {
  # Create a data frame with Sample IDs
  df <- data.frame(Sample = paste0("Sample_", 1:n_samples))
  
  # Simulate Single Trait (e.g., Lodging Score 0-3)
  # Using rounded normal distribution clamped between 0 and 3
  val <- round(rnorm(n_samples, mean = 1.5, sd = 0.8))
  val[val < 0] <- 0
  val[val > 3] <- 3
  df[[traits_single[1]]] <- val
  
  # Simulate Traits with Replicates (e.g., Plant Height, Stem Strength)
  for (trait in traits_with_reps) {
    # Base value for the sample (genetic effect)
    base_value <- runif(n_samples, min = 50, max = 150)
    
    # Generate 3 replicates with random noise (environmental error)
    for (i in 1:3) {
      col_name <- paste0(trait, "_", i)
      noise <- rnorm(n_samples, mean = 0, sd = 5)
      df[[col_name]] <- base_value + noise
    }
  }
  
  return(df)
}

# =========================================================================
# 4. Module: Data Processing & Statistics
# =========================================================================

process_data <- function(df, env_name) {
  
  message(paste0("[Processing] Calculating stats for: ", env_name, "..."))
  
  # --- Step A: Calculate Row Means for Replicates ---
  # Consolidate columns like Trait_A_1, Trait_A_2, Trait_A_3 into Trait_A
  
  df_means <- df %>% select(Sample)
  
  # 1. Handle Single Traits
  for (trait in traits_single) {
    if(trait %in% colnames(df)) {
      df_means[[trait]] <- df[[trait]]
    }
  }
  
  # 2. Handle Replicated Traits (Calculate Mean)
  for (trait in traits_with_reps) {
    # Find columns matching the pattern (Trait_A_1, etc.)
    cols <- grep(paste0("^", trait, "_[0-9]+"), names(df), value = TRUE)
    
    if (length(cols) > 0) {
      df_means[[trait]] <- rowMeans(df[, cols], na.rm = TRUE)
    }
  }
  
  # --- Step B: Calculate Descriptive Statistics ---
  stats_result <- df_means %>%
    pivot_longer(cols = -Sample, names_to = "Trait", values_to = "Value") %>%
    group_by(Trait) %>%
    summarise(
      Env = env_name,
      Count = sum(!is.na(Value)),
      Min = min(Value, na.rm = TRUE),
      Max = max(Value, na.rm = TRUE),
      Mean = mean(Value, na.rm = TRUE),
      SD = sd(Value, na.rm = TRUE),
      CV_percent = (sd(Value, na.rm = TRUE) / mean(Value, na.rm = TRUE)) * 100,
      Skewness = skewness(Value, na.rm = TRUE),
      Kurtosis = kurtosis(Value, na.rm = TRUE),
      .groups = 'drop'
    ) %>%
    select(Env, Trait, everything()) # Reorder columns
  
  return(stats_result)
}

# =========================================================================
# 5. Main Execution Loop
# =========================================================================

message("\n[Start] Starting Multi-Environment Analysis Demo...\n")

all_stats_list <- list()

for (env in env_names) {
  # 1. Generate Data (Simulating reading a CSV)
  dummy_df <- generate_synthetic_data(env, n_samples)
  
  # Optional: You could save this dummy data to show users what input looks like
  # write.csv(dummy_df, file.path(output_dir, paste0(env, "_raw_demo.csv")), row.names=F)
  
  # 2. Process Data
  env_stats <- process_data(dummy_df, env)
  
  # 3. Store Result
  all_stats_list[[env]] <- env_stats
}

# Merge all results
final_stats_df <- bind_rows(all_stats_list)

# =========================================================================
# 6. Save Final Output
# =========================================================================

output_file <- file.path(output_dir, "01_descriptive_statistics.csv")
write.csv(final_stats_df, output_file, row.names = FALSE)

message(paste0("\n[Success] Analysis complete!"))
message(paste0("[Output] Statistics saved to: ", output_file))

# Preview results
print(head(final_stats_df))
