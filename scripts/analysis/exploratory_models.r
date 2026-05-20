##############################
#
# Batch Occupancy Modeling - All Species
# Ian Becker
# December 2025
#
##############################

# This script runs occupancy models for all species with detection matrices
# Automatically processes all species and saves results

library(spOccupancy)
library(dplyr)
library(coda)

####################
#   Setup
####################

cat("=== BATCH OCCUPANCY MODELING ===\n\n")
cat("Starting:", as.character(Sys.time()), "\n\n")

# Create output directory for model results
if (!dir.exists("model_results")) {
  dir.create("model_results")
  cat("Created model_results/ directory\n\n")
}

# Load site covariates once (used by all models)
cat("Loading site covariates...\n")
site_covs <- readRDS("site_covariates.rds")
cat("Loaded covariates for", nrow(site_covs), "sites\n\n")

####################
#   Species List
####################

# Get all species with detection matrices
detection_files <- list.files("detection_matrices/", 
                              pattern = "^detection_matrix_.*\\.rds$")

species_list <- gsub("detection_matrix_|.rds", "", detection_files)
species_list <- gsub("_", " ", species_list)

cat("Found detection matrices for", length(species_list), "species:\n")
for (sp in species_list) {
  cat("  -", sp, "\n")
}
cat("\n")

####################
#   Model Settings
####################

# Common settings for all models
model_settings <- list(
  n_batch = 1000,
  batch_length = 10,
  n_burn = 5000,
  n_thin = 5,
  n_chains = 3,
  n_omp_threads = 1
)

cat("Model settings:\n")
cat("  Total iterations per chain:", model_settings$n_batch * model_settings$batch_length, "\n")
cat("  Burn-in:", model_settings$n_burn, "\n")
cat("  Thinning:", model_settings$n_thin, "\n")
cat("  Chains:", model_settings$n_chains, "\n")
cat("  Posterior samples:", (model_settings$n_batch * model_settings$batch_length - 
                               model_settings$n_burn) / model_settings$n_thin * 
      model_settings$n_chains, "\n\n")

####################
#   Prepare Covariates (Once)
####################

cat("Preparing occupancy covariates...\n")

# Select covariates (7 total - excluding urban_pct and water_pct)
occ_covs <- site_covs %>%
  select(trees_pct, grass_pct, shrub_pct, flooded_veg_pct, 
         crops_pct, water_pct, habitat_diversity, log_area)

# Scale covariates
occ_covs_scaled <- as.data.frame(scale(occ_covs))

cat("Occupancy covariates (n=7):\n")
cat("  -", paste(names(occ_covs_scaled), collapse = "\n  - "), "\n\n")

# Check for NAs
if (any(is.na(occ_covs_scaled))) {
  cat("WARNING: NAs detected in covariates!\n")
  stop("Fix covariate NAs before running models")
}

####################
#   Batch Model Fitting (Simplified)
####################

all_results <- list()
batch_start <- Sys.time()

cat("========================================\n")
cat("STARTING BATCH MODEL FITTING\n")
cat("========================================\n\n")

for (i in seq_along(species_list)) {
  
  species_name <- species_list[i]
  species_filename <- gsub(" ", "_", species_name)
  
  cat("\n[", i, "/", length(species_list), "]", species_name, "... ")
  
  tryCatch({
    
    # Load detection matrix
    detection_matrix <- readRDS(paste0("detection_matrices/detection_matrix_", 
                                       species_filename, ".rds"))
    metadata <- readRDS(paste0("detection_matrices/detection_metadata_", 
                               species_filename, ".rds"))
    
    # Prepare detection covariates
    years <- metadata$years
    n_sites <- metadata$n_sites
    n_years <- metadata$n_years
    n_months <- metadata$n_months
    
    year_array <- array(NA, dim = c(n_sites, n_years, n_months))
    for (j in 1:n_years) {
      year_array[, j, ] <- years[j]
    }
    year_array_scaled <- (year_array - mean(years)) / sd(years)
    
    #Peak season array (Nov-Apr = 1, May-Oct = 0)
    peak_months <- c(1, 2, 3, 4, 11, 12)  # Jan-Apr and Nov-Dec
    peak_season_array <- array(NA, dim = c(n_sites, n_years, n_months))
    for (m in 1:n_months) {
      peak_season_array[, , m] <- ifelse(m %in% peak_months, 1, 0)
    }
    
    # Create data list - ADD peak_season here
    data_list <- list(
      y = detection_matrix,
      occ.covs = occ_covs_scaled,
      det.covs = list(
        year = year_array_scaled,
        peak_season = peak_season_array
      )
    )
    
    # Fit model
    model_start <- Sys.time()
    
    model <- tPGOcc(
      occ.formula = ~ trees_pct + grass_pct + shrub_pct + flooded_veg_pct + 
        crops_pct + water_pct + habitat_diversity + log_area,
      det.formula = ~ year + peak_season,
      data = data_list,
      n.batch = model_settings$n_batch,
      batch.length = model_settings$batch_length,
      n.burn = model_settings$n_burn,
      n.thin = model_settings$n_thin,
      n.chains = model_settings$n_chains,
      n.omp.threads = model_settings$n_omp_threads,
      verbose = FALSE
    )
    
    run_time <- as.numeric(difftime(Sys.time(), model_start, units = "secs"))
    
    # Store results
    all_results[[species_name]] <- list(
      species = species_name,
      model = model,
      metadata = metadata,
      run_time = run_time,
      success = TRUE
    )
    
    # Save individual model
    saveRDS(all_results[[species_name]], 
            paste0("model_results/model_", species_filename, ".rds"))
    
    cat("done (", round(run_time, 1), "s)\n")
    
    # Clean up
    rm(model, detection_matrix, data_list)
    gc(verbose = FALSE)
    
  }, error = function(e) {
    cat("FAILED:", e$message, "\n")
    all_results[[species_name]] <<- list(
      species = species_name,
      success = FALSE,
      error = e$message
    )
  })
}

# Save all results
saveRDS(all_results, paste0("model_results/all_results_", Sys.Date(), ".rds"))

cat("\n\nDone!", sum(sapply(all_results, function(x) x$success)), "/", 
    length(species_list), "models fitted\n")
cat("Total time:", round(difftime(Sys.time(), batch_start, units = "mins"), 1), "min\n")
batch_end <- Sys.time()
total_time <- difftime(batch_end, batch_start, units = "mins")

########## Checking fit and results

for (sp in names(all_results)) {
  
  if (!all_results[[sp]]$success) {
    cat("\n", sp, "- FAILED\n")
    next
  }
  
  cat("\n========================================\n")
  cat(sp, "\n")
  cat("========================================\n")
  
  # Model summary
  print(summary(all_results[[sp]]$model))
  
  # Bayesian p-value
  ppc <- ppcOcc(all_results[[sp]]$model, fit.stat = 'freeman-tukey', group = 1)
  ppc_pval <- mean(ppc$fit.y.rep > ppc$fit.y)
  
  cat("\nBayesian p-value:", round(ppc_pval, 3))
  if (ppc_pval > 0.1 & ppc_pval < 0.9) cat(" ✓\n") else cat(" ⚠\n")
}



