# ==============================================================================
# Project: Phenotypic Trait Normal Distribution Analysis Pipeline
# Description: 
#   This script performs normality checks and visualization for phenotypic traits 
#   across different environments/targets. It generates histograms with fitted 
#   normal distribution curves and saves them as high-resolution images.
#
# Author: [Liang Xiaotian/494382219@qq.com]
# Date: 2025/12/6
# License: MIT
# ==============================================================================

# ------------------------------------------------------------------------------
# Section 1: Environment Setup & Dependency Management
# ------------------------------------------------------------------------------

# List of required packages
required_packages <- c("ggplot2", "dplyr", "purrr")

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

# Create output directory to keep the workspace clean
output_dir <- "Distribution_Plots_Output"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# ------------------------------------------------------------------------------
# Section 2: Data Preparation (Demo Data Generation)
# NOTE: In real analysis, replace this section with read.csv("your_data.csv")
# ------------------------------------------------------------------------------

# Setting seed for reproducibility
set.seed(123)

# 2.1 Generate Synthetic Data (To protect privacy of original data)
# Simulating 300 samples across 3 environments
n_samples <- 300
demo_targets <- c("Env_2024_A", "Env_2024_B", "Env_2025_A")
demo_traits  <- c("Trait_Length", "Trait_Width", "Trait_Density", "Trait_Yield")

# Create a dummy dataframe
data <- data.frame(
  ID = 1:n_samples,
  Target = sample(demo_targets, n_samples, replace = TRUE)
)

# Add trait columns with random normal distribution data
for (t in demo_traits) {
  data[[t]] <- rnorm(n_samples, mean = runif(1, 10, 50), sd = runif(1, 2, 10))
}

# 2.2 Define Variables for Analysis
# Replace these with your actual column names when using real data
analysis_traits <- demo_traits    # e.g., c("LS", "TIL", ...)
analysis_targets <- demo_targets  # e.g., c("24WJ", "25BC", ...)

# ------------------------------------------------------------------------------
# Section 3: Analysis Logic
# ------------------------------------------------------------------------------

# 3.1 Generate all combinations of Trait x Target
combinations <- expand.grid(
  trait = analysis_traits, 
  target = analysis_targets, 
  stringsAsFactors = FALSE
)

# 3.2 Define the Visualization Function
# Arguments:
#   trait_name: Character, column name of the trait
#   target_name: Character, name of the environment/target
plot_distribution <- function(trait_name, target_name, source_data) {
  
  # Filter data: select specific trait in specific target
  # Use .data[[col]] or all_of() for tidy evaluation
  plot_data <- source_data %>% 
    filter(Target == target_name) %>% 
    select(value = all_of(trait_name), Target) %>% 
    na.omit()
  
  # Skip if no data
  if(nrow(plot_data) < 3) return(NULL)
  
  # Calculate statistics
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
    labs(title = paste(trait_name, "@", target_name),
         subtitle = paste0("N=", nrow(plot_data), " | Mean=", round(mean_val, 2), " | SD=", round(sd_val, 2)),
         x = paste(trait_name, "Value"), 
         y = "Density") +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      plot.subtitle = element_text(hjust = 0.5, color = "dimgray", size = 10),
      panel.grid.minor = element_blank()
    )
  
  return(p)
}

# ------------------------------------------------------------------------------
# Section 4: Execution and Saving
# ------------------------------------------------------------------------------

message("Starting batch plotting...")

# 4.1 Batch Generate Plots
# Using map2 to iterate through combinations
all_plots_list <- map2(
  .x = combinations$trait, 
  .y = combinations$target, 
  .f = ~ plot_distribution(.x, .y, data)
)

# Remove NULLs (in case of empty data)
valid_indices <- which(!sapply(all_plots_list, is.null))
final_plots <- all_plots_list[valid_indices]
final_combos <- combinations[valid_indices, ]

# 4.2 Save Individual Files (PDF & PNG)
# Using pwalk (parallel walk) for side-effect operations (saving)
pwalk(list(final_combos$trait, final_combos$target, final_plots), 
      function(tr, ta, p) {
        
        # Define filenames (sanitize file names to avoid illegal characters)
        base_name <- paste0(output_dir, "/", tr, "_", ta, "_Distribution")
        
        # Save as PNG
        ggsave(filename = paste0(base_name, ".png"), plot = p, 
               width = 6, height = 4, dpi = 300, bg = "white")
        
        # Save as PDF
        ggsave(filename = paste0(base_name, ".pdf"), plot = p, 
               width = 6, height = 4, device = "pdf")
      })

message("Individual plots saved in: ", output_dir)

# 4.3 Save Combined PDF
combined_pdf_path <- paste0(output_dir, "/Combined_All_Distributions.pdf")
pdf(combined_pdf_path, width = 6, height = 4)
  invisible(lapply(final_plots, print))
dev.off()

message("Combined PDF saved as: ", combined_pdf_path)
message("Analysis Complete.")
