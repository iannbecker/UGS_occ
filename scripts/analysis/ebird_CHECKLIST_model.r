##############################
#
# Batch Occupancy Modeling — eBird Checklist Level
# Ian Becker
# August 2026
#
##############################

# Runs landscape-level occupancy models using ebird data
# with checklist-level detection matrices

library(spOccupancy)
library(dplyr)
library(coda)

setwd("~/Desktop/project_code/UGS_occ/data")

# ============================================================================
# 1. SETUP AND LOAD DATA
# ============================================================================

dir.create("model_results_ebird_checklist")

# Load in site covariates

site_covs <- readRDS("site_covariates_ebird.rds")

# Detection matrices file path

detection_files <- list.files("detection_matrices_ebird_checklist/",
                              pattern = "^detection_matrix_.*\\.rds$")

# List of species based on detection matrices

species_list <- gsub("detection_matrix_|.rds", "", detection_files)
species_list <- gsub("_", " ", species_list)

# Check detection matrix species 

length(species_list)

# Load in effort data

effort_lookup <- readRDS("detection_matrices_ebird_checklist/effort_lookup.rds")

# ============================================================================
# 2. PREP MODEL PARAMETERS AND VARIABLES
# ============================================================================

# Model settings

model_settings <- list(
  n_batch       = 1000,
  batch_length  = 10,
  n_burn        = 5000,
  n_thin        = 5,
  n_chains      = 3,
  n_omp_threads = 1
)


# Setting up occupancy covariates 

occ_covs <- site_covs %>%
  select(trees_pct, grass_pct, shrub_pct, flooded_veg_pct,
         crops_pct, water_pct, habitat_diversity, log_area)

# Scale occupancy covariates

occ_covs_scaled <- as.data.frame(scale(occ_covs))

# ============================================================================
# 3. LOOP THROUGH AND FIT MODELS FOR EACH SPECIES 
# ============================================================================

# Setup results structure

all_results <- list()
batch_start <- Sys.time()

# Main loop through species 

for (i in seq_along(species_list)) {
  
  species_name     <- species_list[i]
  species_filename <- gsub(" ", "_", species_name)
  
  cat("\n[", i, "/", length(species_list), "]", species_name, "... ")
  
  tryCatch({
    
    detection_matrix <- readRDS(
      paste0("detection_matrices_ebird_checklist/detection_matrix_",
             species_filename, ".rds"))
    metadata <- readRDS(
      paste0("detection_matrices_ebird_checklist/detection_metadata_",
             species_filename, ".rds"))
    
    years    <- metadata$years
    n_sites  <- metadata$n_sites
    n_years  <- metadata$n_years
    n_occ    <- metadata$n_occ
    site_ids <- metadata$site_ids
    
    ### Year detection covariate 
    
    year_array <- array(NA, dim = c(n_sites, n_years, n_occ))
    for (j in 1:n_years) year_array[, j, ] <- years[j]
    year_array_scaled <- (year_array - mean(years)) / sd(years)
    
    ### Peak season covariate at checklist level

    peak_months       <- c(1, 2, 3, 4, 11, 12)
    peak_season_array <- array(NA, dim = c(n_sites, n_years, n_occ))
    
    ### Effort covariates
    
    duration_array  <- array(NA, dim = c(n_sites, n_years, n_occ))
    distance_array  <- array(NA, dim = c(n_sites, n_years, n_occ))
    observers_array <- array(NA, dim = c(n_sites, n_years, n_occ))
    
    ### Loop through and add effort covariates
    
    for (k in 1:nrow(effort_lookup)) {
      row      <- effort_lookup[k, ]
      site_idx <- which(site_ids == row$site_id)
      year_idx <- which(years == row$year)
      occ_idx  <- row$occ_idx
      
      if (length(site_idx) == 0 || length(year_idx) == 0) next
      if (is.na(occ_idx) || occ_idx < 1 || occ_idx > n_occ) next
      
      peak_season_array[site_idx, year_idx, occ_idx] <-
        ifelse(row$month %in% peak_months, 1, 0)
      
      duration_array [site_idx, year_idx, occ_idx] <- row$duration_minutes
      distance_array [site_idx, year_idx, occ_idx] <- row$effort_distance_km
      observers_array[site_idx, year_idx, occ_idx] <- row$number_observers
    }
    
    ### Scale effort arrays and impute missing values 
    
    scale_array <- function(arr) {
      vals <- arr[!is.na(arr)]
      (arr - mean(vals)) / sd(vals)
    }
    
    impute_array <- function(arr) {
      arr[is.na(arr)] <- 0  
      arr
    }
    
    
    duration_ready  <- impute_array(scale_array(duration_array))
    distance_ready  <- impute_array(scale_array(distance_array))
    observers_ready <- impute_array(scale_array(observers_array))
    peak_season_array[is.na(peak_season_array)] <- 0
    
    ### Prepare data list for occupancy model
    
    data_list <- list(
      y        = detection_matrix,
      occ.covs = occ_covs_scaled,
      det.covs = list(
        year        = year_array_scaled,
        peak_season = peak_season_array,
        duration    = duration_ready,
        distance    = distance_ready,
        observers   = observers_ready
      )
    )
    
    # Start model time
    
    model_start <- Sys.time()
    
    # Model formula (single-species, multi-season)

    model <- tPGOcc(
      occ.formula = ~ trees_pct + grass_pct + shrub_pct + flooded_veg_pct +
        crops_pct + water_pct + habitat_diversity + log_area,
      det.formula = ~ year + peak_season + duration + distance + observers,
      data          = data_list,
      n.batch       = model_settings$n_batch,
      batch.length  = model_settings$batch_length,
      n.burn        = model_settings$n_burn,
      n.thin        = model_settings$n_thin,
      n.chains      = model_settings$n_chains,
      n.omp.threads = model_settings$n_omp_threads,
      verbose       = FALSE
    )
    
    # Calculate run time for each model
    
    run_time <- as.numeric(difftime(Sys.time(), model_start, units = "secs"))
    
    ### Evaluate model fit for each model
    
    ppc_occ  <- ppcOcc(model, fit.stat = "freeman-tukey", group = 1)
    ppc_site <- ppcOcc(model, fit.stat = "freeman-tukey", group = 2)
    
    ppc_pval_occ  <- mean(ppc_occ$fit.y.rep  > ppc_occ$fit.y)
    ppc_pval_site <- mean(ppc_site$fit.y.rep > ppc_site$fit.y)
    
    # Store results for each species
    
    all_results[[species_name]] <- list(
      species        = species_name,
      model          = model,
      metadata       = metadata,
      run_time       = run_time,
      ppc_pval_occ   = ppc_pval_occ,
      ppc_pval_site  = ppc_pval_site,
      success        = TRUE,
      data_source    = "eBird_checklist"
    )
    
    # Save results for each species
    
    saveRDS(all_results[[species_name]],
            paste0("model_results_ebird_checklist/model_",
                   species_filename, ".rds"))
    
    cat("done (", round(run_time, 1), "s) | p-val occ:",
        round(ppc_pval_occ, 3), "| p-val site:", round(ppc_pval_site, 3))
    if (ppc_pval_occ > 0.1 & ppc_pval_occ < 0.9) cat(" ✓\n") else cat(" ⚠\n")
    
    # Cleanup to save space
    
    rm(model, detection_matrix, data_list,
       year_array, peak_season_array,
       duration_array, distance_array, observers_array,
       duration_ready, distance_ready, observers_ready)
    gc(verbose = FALSE)
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

# ============================================================================
# 4. SAVE RESULTS
# ============================================================================

# Save all results

saveRDS(all_results,
        paste0("model_results_ebird_checklist/all_results_ebird_checklist_",
               Sys.Date(), ".rds"))

# Summarize model fit across species

for (sp in names(all_results)) {
  if (!all_results[[sp]]$success) {
    cat(sp, "- FAILED\n")
    next
  }
  p_occ  <- round(all_results[[sp]]$ppc_pval_occ,  3)
  p_site <- round(all_results[[sp]]$ppc_pval_site, 3)
  pass   <- p_occ > 0.1 & p_occ < 0.9
  cat(sprintf("%-35s occ: %.3f  site: %.3f  %s\n",
              sp, p_occ, p_site, ifelse(pass, "✓", "⚠")))
}
