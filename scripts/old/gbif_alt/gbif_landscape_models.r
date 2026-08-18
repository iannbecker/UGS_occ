##############################
#
# Batch Occupancy Modeling - GBIF Data
# Ian Becker
# May 2026
#
##############################

# Runs landscape-level occupancy models using GBIF detection matrices
# Same model structure as original iNat analysis
# All outputs saved to model_results_gbif/ to avoid overwriting iNat results

library(spOccupancy)
library(dplyr)
library(coda)

setwd("~/Desktop/project_code/UGS_occ/data")

####################
#   Setup
####################

cat("=== BATCH OCCUPANCY MODELING — GBIF DATA ===\n\n")
cat("Starting:", as.character(Sys.time()), "\n\n")

if (!dir.exists("model_results_gbif")) {
  dir.create("model_results_gbif")
  cat("Created model_results_gbif/ directory\n\n")
}

####################
#   Load Site Covariates
####################

cat("Loading site covariates...\n")
site_covs_all <- readRDS("site_covariates_gbif.rds")
cat("Loaded covariates for", nrow(site_covs_all), "sites\n\n")

####################
#   Filter to Sites with GBIF Observations
####################

cat("Filtering to sites with GBIF observations...\n")

gbif_with_sites <- read.csv("gbif_in_sites.csv")

sites_with_gbif <- gbif_with_sites %>%
  filter(!is.na(site_id)) %>%
  pull(site_id) %>%
  unique() %>%
  sort()

cat("Sites with GBIF observations:", length(sites_with_gbif), "\n")
cat("Sites without GBIF observations:",
    nrow(site_covs_all) - length(sites_with_gbif), "\n\n")

# Filter covariates to GBIF sites only
site_covs <- site_covs_all %>%
  filter(site_id %in% sites_with_gbif)

saveRDS(site_covs, "site_covariates_gbif.rds")
cat("Saved: site_covariates_gbif.rds\n\n")

####################
#   Species List
####################

detection_files <- list.files("detection_matrices_gbif/",
                              pattern = "^detection_matrix_.*\\.rds$")

species_list <- gsub("detection_matrix_|.rds", "", detection_files)
species_list <- gsub("_", " ", species_list)

cat("Found detection matrices for", length(species_list), "species:\n")
for (sp in species_list) cat("  -", sp, "\n")
cat("\n")

####################
#   Model Settings
####################

model_settings <- list(
  n_batch       = 1000,
  batch_length  = 10,
  n_burn        = 5000,
  n_thin        = 5,
  n_chains      = 3,
  n_omp_threads = 1
)

cat("Model settings:\n")
cat("  Total iterations per chain:",
    model_settings$n_batch * model_settings$batch_length, "\n")
cat("  Burn-in:", model_settings$n_burn, "\n")
cat("  Thinning:", model_settings$n_thin, "\n")
cat("  Chains:", model_settings$n_chains, "\n")
cat("  Posterior samples:",
    (model_settings$n_batch * model_settings$batch_length - model_settings$n_burn) /
      model_settings$n_thin * model_settings$n_chains, "\n\n")

####################
#   Prepare Covariates
####################

cat("Preparing occupancy covariates...\n")

occ_covs <- site_covs %>%
  select(trees_pct, grass_pct, shrub_pct, flooded_veg_pct,
         crops_pct, water_pct, habitat_diversity, log_area)

occ_covs_scaled <- as.data.frame(scale(occ_covs))

cat("Occupancy covariates (n=8):\n")
cat("  -", paste(names(occ_covs_scaled), collapse = "\n  - "), "\n\n")

if (any(is.na(occ_covs_scaled))) {
  cat("WARNING: NAs detected in covariates!\n")
  stop("Fix covariate NAs before running models")
}

####################
#   Batch Model Fitting
####################

all_results <- list()
batch_start <- Sys.time()

cat("========================================\n")
cat("STARTING BATCH MODEL FITTING\n")
cat("========================================\n\n")

for (i in seq_along(species_list)) {
  
  species_name     <- species_list[i]
  species_filename <- gsub(" ", "_", species_name)
  
  cat("\n[", i, "/", length(species_list), "]", species_name, "... ")
  
  tryCatch({
    
    # Load detection matrix from GBIF directory
    detection_matrix <- readRDS(paste0("detection_matrices_gbif/detection_matrix_",
                                       species_filename, ".rds"))
    metadata         <- readRDS(paste0("detection_matrices_gbif/detection_metadata_",
                                       species_filename, ".rds"))
    
    years    <- metadata$years
    n_sites  <- metadata$n_sites
    n_years  <- metadata$n_years
    n_months <- metadata$n_months
    
    # Year detection covariate
    year_array <- array(NA, dim = c(n_sites, n_years, n_months))
    for (j in 1:n_years) year_array[, j, ] <- years[j]
    year_array_scaled <- (year_array - mean(years)) / sd(years)
    
    # Peak season detection covariate (Nov-Apr = 1, May-Oct = 0)
    peak_months       <- c(1, 2, 3, 4, 11, 12)
    peak_season_array <- array(NA, dim = c(n_sites, n_years, n_months))
    for (m in 1:n_months) {
      peak_season_array[, , m] <- ifelse(m %in% peak_months, 1, 0)
    }
    
    data_list <- list(
      y        = detection_matrix,
      occ.covs = occ_covs_scaled,
      det.covs = list(
        year        = year_array_scaled,
        peak_season = peak_season_array
      )
    )
    
    model_start <- Sys.time()
    
    model <- tPGOcc(
      occ.formula = ~ trees_pct + grass_pct + shrub_pct + flooded_veg_pct +
        crops_pct + water_pct + habitat_diversity + log_area,
      det.formula = ~ year + peak_season,
      data          = data_list,
      n.batch       = model_settings$n_batch,
      batch.length  = model_settings$batch_length,
      n.burn        = model_settings$n_burn,
      n.thin        = model_settings$n_thin,
      n.chains      = model_settings$n_chains,
      n.omp.threads = model_settings$n_omp_threads,
      verbose       = FALSE
    )
    
    run_time <- as.numeric(difftime(Sys.time(), model_start, units = "secs"))
    
    all_results[[species_name]] <- list(
      species     = species_name,
      model       = model,
      metadata    = metadata,
      run_time    = run_time,
      success     = TRUE,
      data_source = "GBIF"
    )
    
    # Save to gbif directory
    saveRDS(all_results[[species_name]],
            paste0("model_results_gbif/model_", species_filename, ".rds"))
    
    cat("done (", round(run_time, 1), "s)\n")
    
    rm(model, detection_matrix, data_list)
    gc(verbose = FALSE)
    
  }, error = function(e) {
    cat("FAILED:", e$message, "\n")
    all_results[[species_name]] <<- list(
      species = species_name,
      success = FALSE,
      error   = e$message
    )
  })
}

# Save combined results
saveRDS(all_results,
        paste0("model_results_gbif/all_results_gbif_", Sys.Date(), ".rds"))

cat("\n\nDone!", sum(sapply(all_results, function(x) x$success)), "/",
    length(species_list), "models fitted\n")
cat("Total time:",
    round(difftime(Sys.time(), batch_start, units = "mins"), 1), "min\n\n")

####################
#   Check Fit and Results
####################

for (sp in names(all_results)) {
  
  if (!all_results[[sp]]$success) {
    cat("\n", sp, "- FAILED\n")
    next
  }
  
  cat("\n========================================\n")
  cat(sp, "\n")
  cat("========================================\n")
  
  print(summary(all_results[[sp]]$model))
  
  ppc      <- ppcOcc(all_results[[sp]]$model, fit.stat = "freeman-tukey", group = 1)
  ppc_pval <- mean(ppc$fit.y.rep > ppc$fit.y)
  
  cat("\nBayesian p-value:", round(ppc_pval, 3))
  if (ppc_pval > 0.1 & ppc_pval < 0.9) cat(" ✓\n") else cat(" ⚠\n")
}