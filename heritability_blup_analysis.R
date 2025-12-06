# ==============================================================================
# Project: Heritability (H²) and BLUP Calculation Pipeline
# Description: 
#   This script fits Linear Mixed Models (LMM) to phenotypic data to calculate 
#   Best Linear Unbiased Predictions (BLUPs) and Broad-Sense Heritability (H²).
#   It performs model comparison (BIC) and Likelihood Ratio Tests (LRT) for 
#   variance components. Includes synthetic data generation for reproducibility.
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
required_packages <- c("lme4", "emmeans", "data.table", "tidyverse", "writexl")

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
output_dir <- "Heritability_BLUP_Output"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# ------------------------------------------------------------------------------
# Section 2: Data Simulation (Privacy Protection)
# NOTE: In real analysis, replace this section with read.csv("your_data.csv")
# ------------------------------------------------------------------------------

message("Generating synthetic phenotypic data...")

# Set seed for reproducibility
set.seed(123)

# 2.1 Simulation Parameters
n_lines <- 100  # Number of genotypes
n_years <- 2    # Number of years
n_locs  <- 3    # Number of locations
n_reps  <- 2    # Replications per environment

# Generate experimental design structure
dat <- expand.grid(
  Line = paste0("Line_", 1:n_lines),
  Year = paste0("Year_", 1:n_years),
  Loc  = paste0("Loc_", 1:n_locs),
  Rep  = 1:n_reps
)

# Define traits (Generalized names)
traits <- c("Trait_1", "Trait_2", "Trait_3", "Trait_4")

# Simulate trait values with variance components
# Model: y = mu + Line + Year + Loc + Year:Loc + e
for (t in traits) {
  # Random effects
  eff_line <- rnorm(n_lines, mean = 0, sd = runif(1, 5, 10))
  names(eff_line) <- paste0("Line_", 1:n_lines)
  
  eff_year <- rnorm(n_years, mean = 0, sd = 2)
  names(eff_year) <- paste0("Year_", 1:n_years)
  
  eff_loc <- rnorm(n_locs, mean = 0, sd = 4)
  names(eff_loc) <- paste0("Loc_", 1:n_locs)
  
  # Construct phenotype
  dat[[t]] <- 50 + # Intercept
    eff_line[dat$Line] +
    eff_year[dat$Year] +
    eff_loc[dat$Loc] +
    rnorm(nrow(dat), mean = 0, sd = 3) # Residual error
}

# Convert categorical variables to factors
cols_to_factor <- c("Line", "Year", "Loc")
dat[cols_to_factor] <- lapply(dat[cols_to_factor], as.factor)

# ------------------------------------------------------------------------------
# Section 3: Model Fitting and Analysis
# ------------------------------------------------------------------------------

message("Starting BLUP and Heritability analysis...")

# Initialize storage lists
best_blup_list <- list()
model_compare_list <- list()
heritability_list <- list()

for (trait in traits) {
  message("Processing trait: ", trait)
  
  # 3.1 Define Linear Mixed Models (LMM)
  # Model 1: Line + Year random
  model1 <- lmer(as.formula(paste(trait, "~ 1 + (1 | Line) + (1 | Year)")), 
                 data = dat, control = lmerControl(optimizer = "bobyqa"))
  
  # Model 2: Line + Loc random
  model2 <- lmer(as.formula(paste(trait, "~ 1 + (1 | Line) + (1 | Loc)")), 
                 data = dat, control = lmerControl(optimizer = "bobyqa"))
  
  # Model 3: Line + Year + Loc random
  model3 <- lmer(as.formula(paste(trait, "~ 1 + (1 | Line) + (1 | Year) + (1 | Loc)")), 
                 data = dat, control = lmerControl(optimizer = "bobyqa"))
  
  # Model 4: Line + Year + Loc + Year:Loc random (Full model)
  model4 <- lmer(as.formula(paste(trait, "~ 1 + (1 | Line) + (1 | Year) + (1 | Loc) + (1 | Year:Loc)")), 
                 data = dat, control = lmerControl(optimizer = "bobyqa"))
  
  # 3.2 Model Comparison (BIC)
  bic_values <- c(BIC(model1), BIC(model2), BIC(model3), BIC(model4))
  model_names <- c("Model1", "Model2", "Model3", "Model4")
  models_list <- list(model1, model2, model3, model4)
  
  model_compare <- data.frame(Model = model_names, BIC = bic_values)
  model_compare_list[[trait]] <- model_compare
  
  # Select Best Model
  best_idx <- which.min(bic_values)
  best_model <- models_list[[best_idx]]
  best_model_name <- model_names[best_idx]
  
  # 3.3 Likelihood Ratio Test (LRT) for Genetic Variance
  # Null model: Remove (1 | Line) from best model
  best_formula <- formula(best_model)
  null_formula <- update(best_formula, . ~ . - (1 | Line))
  null_model <- lmer(null_formula, data = dat, control = lmerControl(optimizer = "bobyqa"))
  
  lrt_res <- anova(best_model, null_model, refit = FALSE)
  lrt_chisq <- lrt_res$Chisq[2]
  lrt_pvalue <- 1 - pchisq(lrt_chisq, df = 1) # 1 DF difference
  
  # 3.4 Variance Components Extraction
  var_comp <- as.data.frame(VarCorr(best_model)) %>%
    select(grp, vcov) %>%
    rename(Component = grp, Variance = vcov)
  
  # Helper to extract variance safely
  get_var <- function(comp_name) {
    val <- var_comp$Variance[var_comp$Component == comp_name]
    if(length(val) == 0) return(0) else return(val)
  }
  
  sigma2_G <- get_var("Line")
  sigma2_Year <- get_var("Year")
  sigma2_Loc <- get_var("Loc")
  sigma2_YearLoc <- get_var("Year:Loc")
  sigma2_e <- get_var("Residual")
  
  # 3.5 Calculate Broad-Sense Heritability (H2)
  # Experimental design parameters
  n_y <- length(unique(dat$Year))
  n_l <- length(unique(dat$Loc))
  # Calculate harmonic mean of reps (for unbalanced data)
  n_r <- dat %>% 
    group_by(Line, Year, Loc) %>% summarise(n = n(), .groups = "drop") %>% 
    pull(n) %>% mean()
  
  # Phenotypic Variance (Vp)
  # Formula: Vg + V_y/y + V_l/l + V_yl/yl + Ve/ylr
  Vp <- sigma2_G + (sigma2_Year/n_y) + (sigma2_Loc/n_l) + 
        (sigma2_YearLoc/(n_y*n_l)) + (sigma2_e/(n_y*n_l*n_r))
  
  H2 <- sigma2_G / Vp
  
  # Store Results
  heritability_df <- data.frame(
    Trait = trait,
    Best_Model = best_model_name,
    Sigma2_G = round(sigma2_G, 4),
    Sigma2_Year = round(sigma2_Year, 4),
    Sigma2_Loc = round(sigma2_Loc, 4),
    Sigma2_YearLoc = round(sigma2_YearLoc, 4),
    Sigma2_Residual = round(sigma2_e, 4),
    Total_Vp = round(Vp, 4),
    H2 = round(H2, 3),
    LRT_ChiSq = round(lrt_chisq, 3),
    LRT_Pvalue = format(lrt_pvalue, scientific = TRUE, digits = 6)
  )
  heritability_list[[trait]] <- heritability_df
  
  # 3.6 Extract BLUP Values
  blups <- ranef(best_model)$Line %>%
    rownames_to_column("Line") %>%
    rename(BLUP = `(Intercept)`) %>%
    mutate(Trait = trait)
  
  best_blup_list[[trait]] <- blups
}

# ------------------------------------------------------------------------------
# Section 4: Result Compilation and Export
# ------------------------------------------------------------------------------

# 4.1 Combine BLUPs (Wide Format)
# Start with the first trait's BLUPs
blup_combined <- best_blup_list[[1]] %>% select(Line, BLUP)
colnames(blup_combined)[2] <- paste0("BLUP_", traits[1])

if (length(traits) > 1) {
  for (i in 2:length(traits)) {
    temp <- best_blup_list[[i]] %>% select(Line, BLUP)
    colnames(temp)[2] <- paste0("BLUP_", traits[i])
    blup_combined <- merge(blup_combined, temp, by = "Line", all = TRUE)
  }
}

# 4.2 Combine Other Results
model_comparison_all <- bind_rows(model_compare_list, .id = "Trait")
heritability_all <- bind_rows(heritability_list)

# 4.3 Export to Excel
output_file <- paste0(output_dir, "/BLUP_Heritability_Results.xlsx")

write_xlsx(
  list(
    BLUP_Values = blup_combined,
    Heritability_LRT = heritability_all,
    Model_Comparison = model_comparison_all
  ),
  path = output_file
)

message("Analysis Complete.")
message("Results saved to: ", output_file)
