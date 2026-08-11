##############################
#
# Batch Occupancy Modeling — eBird Checklist Level
# Ian Becker
# August 2026
#
##############################

# Runs landscape-level occupancy models using checklist-level detection matrices
# Secondary occasions = individual checklists (subsampled, not months)
# Effort covariates at checklist level — no aggregation needed
# Both occasion-level (group=1) and site-level (group=2) PPC reported
# All outputs saved to model_results_ebird_checklist/

library(spOccupancy)
library(dplyr)
library(coda)

setwd("~/Desktop/project_code/UGS_occ/data")

# ============================================================================
# 1. SETUP
# ============================================================================

cat("=== BATCH OCCUPANCY MODELING — eBird CHECKLIST LEVEL ===\n\n")
cat("Starting:", as.character(Sys.time()), "\n\n")

if (!dir.exists("model_results_ebird_checklist")) {
  dir.create("model_results_ebird_checklist")
  cat("Created model_results_ebird_checklist/ directory\n\n")
}

# ============================================================================
# 2. LOAD SITE COVARIATES
# ============================================================================

cat("Loading site covariates...\n")
site_covs <- readRDS("site_covariates_ebird.rds")
cat("Loaded covariates for", nrow(site_covs), "sites\n\n")

# ============================================================================
# 3. SPECIES LIST
# ============================================================================

detection_files <- list.files("detection_matrices_ebird_checklist/",
                              pattern = "^detection_matrix_.*\\.rds$")

species_list <- gsub("detection_matrix_|.rds", "", detection_files)
species_list <- gsub("_", " ", species_list)

cat("Found detection matrices for", length(species_list), "species:\n")
for (sp in species_list) cat("  -", sp, "\n")
cat("\n")

# ============================================================================
# 4. LOAD EFFORT LOOKUP
# ============================================================================

cat("Loading effort lookup...\n")
effort_lookup <- readRDS("detection_matrices_ebird_checklist/effort_lookup.rds")
cat("Effort lookup rows:", nrow(effort_lookup), "\n\n")

# ============================================================================
# 5. MODEL SETTINGS
# ============================================================================

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
cat("  Chains:", model_settings$n_chains, "\n\n")

# ============================================================================
# 6. PREPARE OCCUPANCY COVARIATES
# ============================================================================

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

# ============================================================================
# 7. BATCH MODEL FITTING
# ============================================================================

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
    
    # ── Year detection covariate ──────────────────────────────────────────────
    # One value per site-year-occasion — same year value across all occasions
    year_array <- array(NA, dim = c(n_sites, n_years, n_occ))
    for (j in 1:n_years) year_array[, j, ] <- years[j]
    year_array_scaled <- (year_array - mean(years)) / sd(years)
    
    # ── Month (peak season) covariate at checklist level ─────────────────────
    # Built from effort lookup — each occasion has a month value
    peak_months       <- c(1, 2, 3, 4, 11, 12)
    peak_season_array <- array(NA, dim = c(n_sites, n_years, n_occ))
    
    duration_array  <- array(NA, dim = c(n_sites, n_years, n_occ))
    distance_array  <- array(NA, dim = c(n_sites, n_years, n_occ))
    observers_array <- array(NA, dim = c(n_sites, n_years, n_occ))
    
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
    
    # ── Scale effort arrays ───────────────────────────────────────────────────
    scale_array <- function(arr) {
      vals <- arr[!is.na(arr)]
      (arr - mean(vals)) / sd(vals)
    }
    
    impute_array <- function(arr) {
      arr[is.na(arr)] <- 0  # 0 = mean after scaling
      arr
    }
    
    duration_ready  <- impute_array(scale_array(duration_array))
    distance_ready  <- impute_array(scale_array(distance_array))
    observers_ready <- impute_array(scale_array(observers_array))
    peak_season_array[is.na(peak_season_array)] <- 0
    
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
    
    model_start <- Sys.time()
    
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
    
    run_time <- as.numeric(difftime(Sys.time(), model_start, units = "secs"))
    
    # ── Model fit — both occasion and site level ───────────────────────────────
    ppc_occ  <- ppcOcc(model, fit.stat = "freeman-tukey", group = 1)
    ppc_site <- ppcOcc(model, fit.stat = "freeman-tukey", group = 2)
    
    ppc_pval_occ  <- mean(ppc_occ$fit.y.rep  > ppc_occ$fit.y)
    ppc_pval_site <- mean(ppc_site$fit.y.rep > ppc_site$fit.y)
    
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
    
    saveRDS(all_results[[species_name]],
            paste0("model_results_ebird_checklist/model_",
                   species_filename, ".rds"))
    
    cat("done (", round(run_time, 1), "s) | p-val occ:",
        round(ppc_pval_occ, 3), "| p-val site:", round(ppc_pval_site, 3))
    if (ppc_pval_occ > 0.1 & ppc_pval_occ < 0.9) cat(" ✓\n") else cat(" ⚠\n")
    
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
# 8. SAVE AND SUMMARISE
# ============================================================================

saveRDS(all_results,
        paste0("model_results_ebird_checklist/all_results_ebird_checklist_",
               Sys.Date(), ".rds"))

cat("\n\nDone!", sum(sapply(all_results, function(x) x$success)), "/",
    length(species_list), "models fitted\n")
cat("Total time:",
    round(difftime(Sys.time(), batch_start, units = "mins"), 1), "min\n\n")

cat("=== FIT SUMMARY ===\n\n")
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
