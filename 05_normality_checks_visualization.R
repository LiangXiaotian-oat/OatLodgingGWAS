# ==============================================================================
# Script Name: 05_normality_checks_visualization.R
# Description: 
#    This script performs normality checks and visualization for phenotypic traits.
#    It generates histograms with fitted normal distribution curves for each 
#    Trait x Environment combination to assess data distribution.
#
#    Outputs: 
#    1. Individual PNG/PDF plots for each combination.
#    2. A single combined PDF containing all plots for easy reviewing.
#
# Author:      Liang Xiaotian
# Email:       494382219@qq.com
# Date:        2025/12/06
# License:     MIT
# ==============================================================================

# ------------------------------------------------------------------------------
# Section 1: Environment Setup & Dependency Management
# ------------------------------------------------------------------------------

# List of required packages
required_packages <- c("ggplot2", "dplyr", "purrr")

# Function to check and install packages automatically (Consistent with File 01-04)
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
# Using a specific subdirectory in 'results' to keep the project clean
output_dir <- file.path("results", "05_distribution_plots")
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# ------------------------------------------------------------------------------
# Section 2: Data Preparation (Synthetic Data)
# ------------------------------------------------------------------------------

message("\n[Data] Generating synthetic data for normality check...")

# Set seed for reproducibility (Consistent with previous files)
set.seed(2025)

# 2.1 Generate Synthetic Data
n_samples <- 300
demo_targets <- c("Env_A", "Env_B", "Env_C")
# Using generic trait names consistent with File 02/03
demo_traits  <- c("Trait_1", "Trait_2", "Trait_3", "Trait_4")

# Create a dummy dataframe
data <- data.frame(
  ID = 1:n_samples,
  Target = sample(demo_targets, n_samples, replace = TRUE)
)

# Add trait columns with random normal distribution data
for (t in demo_traits) {
  # Simulate normal distribution with random Mean (10-50) and SD (2-10)
  data[[t]] <- rnorm(n_samples, mean = runif(1, 10, 50), sd = runif(1, 2, 10))
}

# 2.2 Define Variables for Analysis
analysis_traits  <- demo_traits   
analysis_targets <- demo_targets  

# ------------------------------------------------------------------------------
# Section 3: Visualization Logic
# ------------------------------------------------------------------------------

# 3.1 Generate all combinations of Trait x Target
combinations <- expand.grid(
  trait = analysis_traits, 
  target = analysis_targets, 
  stringsAsFactors = FALSE
)

# 3.2 Define the Visualization Function
plot_distribution <- function(trait_name, target_name, source_data) {
  
  # Filter data: select specific trait in specific target
  plot_data <- source_data %>% 
    filter(Target == target_name) %>% 
    select(value = all_of(trait_name), Target) %>% 
    na.omit()
  
  # Skip if insufficient data points
  if(nrow(plot_data) < 10) return(NULL)
  
  # Calculate statistics for annotation
  mean_val <- mean(plot_data$value, na.rm = TRUE)
  sd_val <- sd(plot_data$value, na.rm = TRUE)
  
  # Generate Plot
  p <- ggplot(plot_data, aes(x = value)) +
    # Histogram
    geom_histogram(aes(y = after_stat(density)), 
                   fill = "#66c2a5", alpha = 0.7, color = "black", bins = 20) +
    # Normal Distribution Curve
    stat_function(fun = dnorm, 
                  args = list(mean = mean_val, sd = sd_val), 
                  color = "#e41a1c", linewidth = 1.2) +
    # Labels and Themes
    labs(title = paste0(trait_name, " in ", target_name),
         subtitle = paste0("N=", nrow(plot_data), " | Mean=", round(mean_val, 2), " | SD=", round(sd_val, 2)),
         x = trait_name, 
         y = "Density") +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      plot.subtitle = element_text(hjust = 0.5, color = "dimgray", size = 10),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black")
    )
  
  return(p)
}

# ------------------------------------------------------------------------------
# Section 4: Execution and Saving
# ------------------------------------------------------------------------------

message("\n[Plotting] Starting batch plotting...")

# 4.1 Batch Generate Plots
# Using map2 to iterate through combinations
all_plots_list <- map2(
  .x = combinations$trait, 
  .y = combinations$target, 
  .f = ~ plot_distribution(.x, .y, data)
)

# Filter out NULLs (failed plots)
valid_indices <- which(!sapply(all_plots_list, is.null))
final_plots <- all_plots_list[valid_indices]
final_combos <- combinations[valid_indices, ]

if (length(final_plots) == 0) {
  stop("No valid plots generated! Check data input.")
}

# 4.2 Save Individual Files (PNG)
pwalk(list(final_combos$trait, final_combos$target, final_plots), 
      function(tr, ta, p) {
        
        # Define clean filename
        base_name <- file.path(output_dir, paste0(tr, "_", ta, "_Dist"))
        
        # Save as PNG
        ggsave(filename = paste0(base_name, ".png"), plot = p, 
               width = 6, height = 4, dpi = 300, bg = "white")
      })

message(paste("[Output] Individual plots saved to:", output_dir))

# 4.3 Save Combined PDF (Optional but useful)
combined_pdf_path <- file.path(output_dir, "05_Combined_Distribution_Plots.pdf")
pdf(combined_pdf_path, width = 6, height = 4)
  invisible(lapply(final_plots, print))
dev.off()

message(paste("[Output] Combined PDF saved to:", combined_pdf_path))
message("\n[Success] File 05 Analysis Complete.")
