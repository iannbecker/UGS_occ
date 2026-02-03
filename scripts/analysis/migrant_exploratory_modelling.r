##############################
#
# Spring Migrant Occupancy Models
# Ian Becker
# January 2025
#
##############################

library(spOccupancy)
library(dplyr)
library(coda)

cat("=== SPRING MIGRANT OCCUPANCY MODELING ===\n\n")
cat("Starting:", as.character(Sys.time()), "\n\n")

####################
#   Setup
####################

# Load site covariates
cat("Loading site covariates...\n")
site_covs <- readRDS("site_covariates.rds")
cat("Loaded covariates for", nrow(site_covs), "sites\n\n")

# Prepare occupancy covariates
occ_covs <- site_covs %>%
  select(trees_pct, grass_pct, shrub_pct, flooded_veg_pct, 
         crops_pct, habitat_diversity, log_area)

occ_covs_scaled <- as.data.frame(scale(occ_covs))

cat("Occupancy covariates (n=7):\n")
cat("  -", paste(names(occ_covs_scaled), collapse = "\n  - "), "\n\n")

####################
#   Model Settings
####################

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
cat("  Chains:", model_settings$n_chains, "\n\n")

####################
#   Species List
####################

spring_migrants <- c(
  "Summer_Tanager",          # 5 sites, Mar-May
  "Baltimore_Oriole",        # 6 sites, Apr-May
  "Black-and-white_Warbler"  # 4 sites, Mar-May
)

cat("Spring migrants to model:", length(spring_migrants), "\n")
for (sp in spring_migrants) {
  cat("  -", gsub("_", " ", sp), "\n")
}
cat("\n")

####################
#   Batch Model Fitting
####################

all_results <- list()
batch_start <- Sys.time()

cat("========================================\n")
cat("STARTING BATCH MODEL FITTING\n")
cat("========================================\n\n")

for (species_filename in spring_migrants) {
  
  species_name <- gsub("_", " ", species_filename)
  
  cat("\n========================================\n")
  cat(species_name, "\n")
  cat("========================================\n\n")
  
  tryCatch({
    
    # Load detection matrix and metadata
    detection_matrix <- readRDS(paste0("detection_matrices/detection_matrix_spring_", 
                                       species_filename, ".rds"))
    metadata <- readRDS(paste0("detection_matrices/detection_metadata_spring_", 
                               species_filename, ".rds"))
    
    cat("Loaded detection matrix\n")
    cat("  Sites:", metadata$n_sites, "\n")
    cat("  Years:", metadata$n_years, "\n")
    cat("  Months:", metadata$n_months, "(", paste(metadata$month_names, collapse = ", "), ")\n")
    cat("  Window:", metadata$phenology_window, "\n")
    cat("  Observations:", metadata$n_observations, "\n")
    cat("  Sites detected:", metadata$n_sites_detected, "\n")
    cat("  Naive occupancy:", round(metadata$naive_occupancy, 3), "\n\n")
    
    # Get dimensions
    years <- metadata$years
    n_sites <- metadata$n_sites
    n_years <- metadata$n_years
    n_months <- metadata$n_months
    
    # Warning for sparse data
    if (metadata$n_sites_detected < 5) {
      cat("⚠️  WARNING: Only", metadata$n_sites_detected, "sites with detections!\n")
      cat("   Model may not converge or have very high uncertainty.\n\n")
    }
    
    ####################
    #   Detection Covariates
    ####################
    
    cat("Preparing detection covariates...\n")
    
    # Year covariate
    year_array <- array(NA, dim = c(n_sites, n_years, n_months))
    for (j in 1:n_years) {
      year_array[, j, ] <- years[j]
    }
    year_array_scaled <- (year_array - mean(years)) / sd(years)
    
    # Month covariate (species-specific)
    month_array <- array(NA, dim = c(n_sites, n_years, n_months))
    
    if (n_months == 2) {
      # Baltimore Oriole: April (0), May (1)
      month_array[, , 1] <- 0
      month_array[, , 2] <- 1
      cat("  Month coding: April=0, May=1\n")
      
    } else if (n_months == 3) {
      # SUTA & BWWA: March (0), April (1), May (2)
      month_array[, , 1] <- 0
      month_array[, , 2] <- 1
      month_array[, , 3] <- 2
      cat("  Month coding: March=0, April=1, May=2\n")
    }
    
    # Scale month (centers on middle of window)
    month_vals <- 0:(n_months - 1)
    month_array_scaled <- (month_array - mean(month_vals)) / sd(month_vals)
    
    cat("  Month covariate: scaled\n")
    cat("  Year covariate: scaled\n\n")
    
    ####################
    #   Create Data List
    ####################
    
    data_list <- list(
      y = detection_matrix,
      occ.covs = occ_covs_scaled,
      det.covs = list(
        year = year_array_scaled,
        month = month_array_scaled
      )
    )
    
    cat("Data list created\n\n")
    
    ####################
    #   Fit Model
    ####################
    
    cat("Fitting model...\n")
    cat("(This may take 5-15 minutes depending on data)\n\n")
    
    model_start <- Sys.time()
    
    model <- tPGOcc(
      occ.formula = ~ trees_pct + grass_pct + shrub_pct + flooded_veg_pct + 
        crops_pct + habitat_diversity + log_area,
      det.formula = ~ year + month,
      data = data_list,
      n.batch = model_settings$n_batch,
      batch.length = model_settings$batch_length,
      n.burn = model_settings$n_burn,
      n.thin = model_settings$n_thin,
      n.chains = model_settings$n_chains,
      n.omp.threads = model_settings$n_omp_threads,
      verbose = TRUE
    )
    
    run_time <- as.numeric(difftime(Sys.time(), model_start, units = "secs"))
    
    cat("\n✓ Model fitted (", round(run_time, 1), "s)\n\n")
    
    ####################
    #   Convergence Diagnostics
    ####################
    
    cat("=== CONVERGENCE DIAGNOSTICS ===\n\n")
    
    # Extract Rhat values
    rhat_vals <- model$rhat
    
    # Occupancy parameters
    beta_rhat <- rhat_vals$beta.samples
    max_beta_rhat <- max(beta_rhat)
    
    cat("Occupancy parameters (beta):\n")
    cat("  Max Rhat:", round(max_beta_rhat, 3))
    if (max_beta_rhat < 1.1) cat(" ✓\n") else cat(" ⚠️  (>1.1!)\n")
    
    # Detection parameters  
    alpha_rhat <- rhat_vals$alpha.samples
    max_alpha_rhat <- max(alpha_rhat)
    
    cat("Detection parameters (alpha):\n")
    cat("  Max Rhat:", round(max_alpha_rhat, 3))
    if (max_alpha_rhat < 1.1) cat(" ✓\n") else cat(" ⚠️  (>1.1!)\n")
    
    # Occupancy probabilities
    psi_rhat <- rhat_vals$psi.samples
    max_psi_rhat <- max(psi_rhat)
    
    cat("Occupancy probabilities (psi):\n")
    cat("  Max Rhat:", round(max_psi_rhat, 3))
    if (max_psi_rhat < 1.1) cat(" ✓\n") else cat(" ⚠️  (>1.1!)\n")
    
    # Overall assessment
    all_converged <- max_beta_rhat < 1.1 & max_alpha_rhat < 1.1 & max_psi_rhat < 1.1
    
    cat("\nOverall convergence: ")
    if (all_converged) {
      cat("✓ All parameters converged (Rhat < 1.1)\n\n")
    } else {
      cat("⚠️  Some parameters did not converge!\n")
      cat("   Consider longer chains or checking priors.\n\n")
    }
    
    ####################
    #   Model Summary
    ####################
    
    cat("=== MODEL SUMMARY ===\n\n")
    print(summary(model))
    
    ####################
    #   Bayesian P-Value
    ####################
    
    cat("\n=== MODEL FIT (BAYESIAN P-VALUE) ===\n\n")
    
    ppc <- ppcOcc(model, fit.stat = 'freeman-tukey', group = 1)
    ppc_pval <- mean(ppc$fit.y.rep > ppc$fit.y)
    
    cat("Bayesian p-value:", round(ppc_pval, 3))
    
    if (ppc_pval > 0.1 & ppc_pval < 0.9) {
      cat(" ✓ (Good fit)\n\n")
    } else {
      cat(" ⚠️  (Poor fit - p-value should be 0.1-0.9)\n\n")
    }
    
    ####################
    #   Store Results
    ####################
    
    all_results[[species_name]] <- list(
      species = species_name,
      model = model,
      metadata = metadata,
      run_time = run_time,
      converged = all_converged,
      max_rhat = max(c(max_beta_rhat, max_alpha_rhat, max_psi_rhat)),
      bayesian_pval = ppc_pval,
      success = TRUE
    )
    
    # Save individual model
    saveRDS(all_results[[species_name]], 
            paste0("model_results/model_spring_", species_filename, ".rds"))
    
    cat("Saved: model_spring_", species_filename, ".rds\n")
    
    # Clean up
    rm(model, detection_matrix, data_list, ppc)
    gc(verbose = FALSE)
    
  }, error = function(e) {
    cat("✗ ERROR:", e$message, "\n\n")
    all_results[[species_name]] <<- list(
      species = species_name,
      success = FALSE,
      error = e$message
    )
  })
}

####################
#   Final Summary
####################

batch_end <- Sys.time()
total_time <- difftime(batch_end, batch_start, units = "mins")

cat("\n========================================\n")
cat("BATCH MODELING COMPLETE\n")
cat("========================================\n\n")

cat("Total time:", round(total_time, 1), "minutes\n\n")

# Summary table
cat("Results summary:\n")
for (sp in names(all_results)) {
  res <- all_results[[sp]]
  if (res$success) {
    status <- if (res$converged) "✓ Converged" else "⚠️  Convergence issues"
    cat(sprintf("  %-25s: %s (Rhat=%.3f, p=%.3f)\n",
                sp, status, res$max_rhat, res$bayesian_pval))
  } else {
    cat(sprintf("  %-25s: ✗ Failed\n", sp))
  }
}

# Save all results
saveRDS(all_results, paste0("model_results/all_spring_migrants_", Sys.Date(), ".rds"))
cat("\nSaved: all_spring_migrants_", Sys.Date(), ".rds\n")

cat("\n=== COMPLETE ===\n")
