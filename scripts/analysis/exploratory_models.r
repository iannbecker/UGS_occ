##############################
#
# Multi-Season Occupancy Modeling
# Ian Becker
#
##############################

# Fit multi-season occupancy models (tPGOcc) to test species
# Model: occupancy ~ habitat covariates
#        detection ~ year

library(spOccupancy)
library(dplyr)
library(coda)

####################
#   Load Data
####################

cat("=== LOADING DATA ===\n\n")

# Load site covariates
site_covs <- readRDS("site_covariates.rds")

cat("Loaded site covariates for", nrow(site_covs), "sites\n")

# Check site_id order
cat("Site IDs range:", min(site_covs$site_id), "to", max(site_covs$site_id), "\n\n")

####################
#   Species to Model
####################

# Choose which species to run
test_species <- "Black-bellied Whistling-Duck"  # Change to "Black-bellied Whistling-Duck" for BBWD

cat("=== MODELING:", test_species, "===\n\n")

# Load detection matrix
species_filename <- gsub(" ", "_", test_species)
detection_matrix <- readRDS(paste0("detection_matrix_", species_filename, ".rds"))
metadata <- readRDS(paste0("detection_metadata_", species_filename, ".rds"))

cat("Loaded detection matrix:", dim(detection_matrix)[1], "sites ×", 
    dim(detection_matrix)[2], "years ×", 
    dim(detection_matrix)[3], "months\n")
cat("Total observations:", metadata$n_observations, "\n")
cat("Sites with detections:", metadata$n_sites_detected, "\n")
cat("Naive occupancy:", round(metadata$naive_occupancy, 3), "\n\n")

####################
#   Prepare Occupancy Covariates
####################

cat("=== PREPARING COVARIATES ===\n\n")

# Select covariates (excluding urban_pct due to collinearity)
occ_covs <- site_covs %>%
  select(trees_pct, grass_pct, shrub_pct, flooded_veg_pct, 
         crops_pct, log_area, habitat_diversity)

# Scale covariates (important for model convergence!)
occ_covs_scaled <- as.data.frame(scale(occ_covs))

cat("Occupancy covariates:\n")
print(names(occ_covs_scaled))
cat("\n")

# Check for NAs
if (any(is.na(occ_covs_scaled))) {
  cat("WARNING: NAs detected in covariates\n")
  cat("Sites with NAs:", which(rowSums(is.na(occ_covs_scaled)) > 0), "\n\n")
}

####################
#   Prepare Detection Covariates
####################

cat("=== PREPARING DETECTION COVARIATES ===\n\n")

# Year effect on detection (iNat usage increases over time)
years <- metadata$years
n_years <- metadata$n_years
n_months <- metadata$n_months
n_sites <- metadata$n_sites

# Create year covariate matrix [sites × years × months]
year_array <- array(NA, dim = c(n_sites, n_years, n_months))

for (i in 1:n_years) {
  year_array[, i, ] <- years[i]
}

# Scale year
year_array_scaled <- (year_array - mean(years)) / sd(years)

det_covs <- list(
  year = year_array_scaled
)

cat("Detection covariate: year (scaled)\n")
cat("Year range:", min(years), "to", max(years), "\n\n")

####################
#   Prepare Data for spOccupancy
####################

cat("=== PREPARING DATA FOR MODEL ===\n\n")

# spOccupancy data format
data_list <- list(
  y = detection_matrix,
  occ.covs = occ_covs_scaled,
  det.covs = det_covs
)

cat("Data structure:\n")
cat("  Detection matrix:", paste(dim(data_list$y), collapse = " × "), "\n")
cat("  Occupancy covariates:", nrow(data_list$occ.covs), "sites ×", 
    ncol(data_list$occ.covs), "variables\n")
cat("  Detection covariates:", length(data_list$det.covs), "variable(s)\n\n")

####################
#   Set Priors
####################

cat("=== SETTING PRIORS ===\n\n")

# Number of occupancy parameters (intercept + covariates)
n_occ_params <- ncol(occ_covs_scaled) + 1

# Priors: Normal(0, 2.72) for fixed effects (weakly informative)
beta_prior_mean <- rep(0, n_occ_params)
beta_prior_var <- rep(2.72, n_occ_params)

alpha_prior_mean <- 0
alpha_prior_var <- 2.72

cat("Using weakly informative priors:\n")
cat("  Occupancy coefficients: Normal(0, 2.72)\n")
cat("  Detection coefficients: Normal(0, 2.72)\n\n")

####################
#   Fit Multi-Season Occupancy Model
####################

cat("=== FITTING MULTI-SEASON OCCUPANCY MODEL ===\n\n")
cat("This may take several minutes...\n\n")

start_time <- Sys.time()

model <- tPGOcc(
  # Occupancy formula
  occ.formula = ~ trees_pct + grass_pct + shrub_pct + 
    flooded_veg_pct + crops_pct + log_area + habitat_diversity,
  
  # Detection formula
  det.formula = ~ year,
  
  # Data
  data = data_list,
  
  # Priors (optional - will use defaults if not specified)
  # beta.prior = list(mu = beta_prior_mean, var = beta_prior_var),
  # alpha.prior = list(mu = alpha_prior_mean, var = alpha_prior_var),
  
  # MCMC settings
  n.batch = 1000,     # Number of batches
  batch.length = 10,  # Length of each batch (total samples = 10,000)
  n.burn = 5000,
  n.thin = 5,
  n.chains = 3,
  
  # Computation
  n.omp.threads = 1,
  verbose = TRUE,
  n.report = 100
)

end_time <- Sys.time()
run_time <- difftime(end_time, start_time, units = "mins")

cat("\nModel fitting complete!\n")
cat("Run time:", round(run_time, 2), "minutes\n\n")

####################
#   Model Summary
####################

cat("=== MODEL SUMMARY ===\n\n")

model_summary <- summary(model)

####################
#   Convergence Diagnostics
####################

cat("\n=== CONVERGENCE DIAGNOSTICS ===\n\n")

# Check Rhat values
beta_rhat <- model_summary$beta[, "Rhat"]
alpha_rhat <- model_summary$alpha[, "Rhat"]

cat("Occupancy parameters (beta):\n")
print(round(beta_rhat, 3))
cat("\nDetection parameters (alpha):\n")
print(round(alpha_rhat, 3))

# Flag convergence issues
if (any(beta_rhat > 1.1) || any(alpha_rhat > 1.1)) {
  cat("\n⚠ WARNING: Some Rhat values > 1.1 (convergence issues)\n")
  cat("   Consider running longer chains or increasing burn-in\n")
} else {
  cat("\n✓ All Rhat values < 1.1 (good convergence)\n")
}

# Check ESS
beta_ess <- model_summary$beta[, "ESS"]
alpha_ess <- model_summary$alpha[, "ESS"]

cat("\nEffective sample sizes:\n")
cat("  Occupancy parameters: min =", min(beta_ess), ", median =", median(beta_ess), "\n")
cat("  Detection parameters: min =", min(alpha_ess), ", median =", median(alpha_ess), "\n")

if (any(beta_ess < 100) || any(alpha_ess < 100)) {
  cat("\n⚠ WARNING: Some ESS < 100 (low effective sample size)\n")
  cat("   Consider running longer chains\n")
} else {
  cat("\n✓ All ESS > 100 (adequate)\n")
}

####################
#   Extract Results
####################

cat("\n=== EXTRACTING RESULTS ===\n\n")

# Mean occupancy across sites and years
psi_samples <- model$psi.samples
mean_occupancy <- mean(psi_samples)

cat("Mean occupancy probability:", round(mean_occupancy, 3), "\n")

# Occupancy by year
psi_by_year <- apply(psi_samples, 2, function(x) {
  # Each column is a site-year combination
  mean(x)
})

cat("Occupancy trend over time: [would need to reshape to show by year]\n\n")

####################
#   Save Results
####################

cat("=== SAVING RESULTS ===\n\n")

# Save model object
saveRDS(model, paste0("model_", species_filename, ".rds"))
cat("Saved: model_", species_filename, ".rds\n")

# Save summary
model_results <- list(
  species = test_species,
  model_summary = model_summary,
  mean_occupancy = mean_occupancy,
  run_time_mins = as.numeric(run_time),
  convergence = list(
    beta_rhat = beta_rhat,
    alpha_rhat = alpha_rhat,
    beta_ess = beta_ess,
    alpha_ess = alpha_ess
  ),
  date_run = Sys.Date()
)

saveRDS(model_results, paste0("results_", species_filename, ".rds"))
cat("Saved: results_", species_filename, ".rds\n\n")

####################
#   Create Coefficient Table
####################

cat("=== COEFFICIENT SUMMARY TABLE ===\n\n")

# Extract coefficients
coef_table <- data.frame(
  parameter = rownames(model_summary$beta),
  mean = model_summary$beta[, "Mean"],
  sd = model_summary$beta[, "SD"],
  ci_lower = model_summary$beta[, "2.5%"],
  ci_upper = model_summary$beta[, "97.5%"],
  rhat = model_summary$beta[, "Rhat"],
  ess = model_summary$beta[, "ESS"],
  significant = !(model_summary$beta[, "2.5%"] < 0 & model_summary$beta[, "97.5%"] > 0)
)

print(coef_table)

write.csv(coef_table, paste0("coefficients_", species_filename, ".csv"), row.names = FALSE)
cat("\nSaved: coefficients_", species_filename, ".csv\n\n")

####################
#   Interpretation
####################

cat("=== QUICK INTERPRETATION ===\n\n")

cat("Significant habitat associations (95% CI excludes 0):\n")
sig_covs <- coef_table %>%
  filter(significant == TRUE, parameter != "(Intercept)")

if (nrow(sig_covs) > 0) {
  for (i in 1:nrow(sig_covs)) {
    direction <- ifelse(sig_covs$mean[i] > 0, "positive", "negative")
    cat("  ", sig_covs$parameter[i], ":", direction, 
        "(β =", round(sig_covs$mean[i], 3), ")\n")
  }
} else {
  cat("  None detected\n")
}

cat("\n")

####################
#   Done
####################

cat("=== MODEL FITTING COMPLETE ===\n")
cat("Results saved. Next steps:\n")
cat("  1. Check convergence diagnostics\n")
cat("  2. Examine coefficient estimates\n")
cat("  3. Run model for additional species\n")
cat("  4. Compare results across species/guilds\n")

ppc_result <- ppcOcc(model, fit.stat = 'freeman-tukey', group = 2)
summary(ppc_result)
