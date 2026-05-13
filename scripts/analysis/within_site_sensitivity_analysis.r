##############################
#
# Within-Site Sensitivity Analysis
# Random Coordinate Jittering (0-100m)
# Ian Becker
# May 2026
#
##############################

# Sensitivity analysis for within-site occupancy models
# Tests robustness of results to positional uncertainty in iNat data
# Randomly jitters all observations 0-100m in a random direction
# Runs full within-site batch model at 100m cell size
# Results stored separately for comparison with primary models

library(sf)
library(terra)
library(dplyr)
library(lubridate)
library(spOccupancy)

####################
#   Settings
####################

cell_size        <- 100    # meters — finer grid for sensitivity analysis
min_viable_covs  <- 3
min_cells_effort <- 10
jitter_seed      <- 42     # for reproducibility
max_jitter_m     <- 100    # maximum jitter distance in meters

cat("=== WITHIN-SITE SENSITIVITY ANALYSIS ===\n")
cat("Cell size:", cell_size, "m\n")
cat("Max jitter distance:", max_jitter_m, "m\n")
cat("Jitter seed:", jitter_seed, "\n\n")

####################
#   Setup Output Directory
####################

sensitivity_dir <- "within_site_models_sensitivity"

if (!dir.exists(sensitivity_dir)) {
  dir.create(sensitivity_dir)
  cat("Created", sensitivity_dir, "directory\n\n")
}

####################
#   Load Data
####################

cat("Loading data...\n")

viability <- read.csv("within_site_viability_assessment.csv")
cat("Loaded viability assessment:", nrow(viability), "combinations\n")

viable <- viability %>%
  filter(n_viable_covariates >= min_viable_covs,
         n_cells_with_effort >= min_cells_effort)
cat("Viable combinations after filtering:", nrow(viable), "\n\n")

sites     <- st_read("lrgv_green_spaces_detection_filtered", quiet = TRUE)
sites_utm <- st_transform(sites, crs = 32614)

inat <- read.csv("inat_observations_with_sites.csv")
inat$date  <- as.Date(inat$observed_on)
inat$year  <- year(inat$date)
inat$month <- month(inat$date)
inat <- inat %>% filter(year >= 2015, year <= 2025)

landcover <- rast("lrgv_dynamic_world_cover.tif")

cat("Data loaded:", nrow(inat), "observations\n\n")

####################
#   Apply Random Jitter
####################

cat("Applying random coordinate jitter (0-", max_jitter_m, "m)...\n")

set.seed(jitter_seed)

n_obs    <- nrow(inat)
distance <- runif(n_obs, 0, max_jitter_m)  # random distance 0-100m
angle    <- runif(n_obs, 0, 2 * pi)         # random direction in radians

# Convert meter offsets to decimal degrees
# 1 degree latitude  ≈ 111,320m
# 1 degree longitude ≈ 111,320m * cos(latitude)
lat_rad <- inat$latitude * pi / 180

inat$longitude_j <- inat$longitude + (distance * sin(angle)) / (111320 * cos(lat_rad))
inat$latitude_j  <- inat$latitude  + (distance * cos(angle)) / 111320

cat("Jitter applied\n")
cat("Mean jitter distance:", round(mean(distance), 1), "m\n")
cat("Max jitter distance:", round(max(distance), 1), "m\n\n")

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

years    <- 2015:2025
n_years  <- length(years)
months   <- 1:12
n_months <- length(months)

cat("Model settings:\n")
cat("  Iterations:", model_settings$n_batch * model_settings$batch_length, "\n")
cat("  Burn-in:", model_settings$n_burn, "\n")
cat("  Chains:", model_settings$n_chains, "\n\n")

####################
#   Results Storage
####################

all_results <- list()
summary_results <- data.frame(
  species          = character(),
  site_id          = integer(),
  n_obs            = integer(),
  n_cells          = integer(),
  n_cells_effort   = integer(),
  n_cells_detected = integer(),
  naive_occ        = numeric(),
  n_covariates     = integer(),
  covariates_used  = character(),
  model_converged  = logical(),
  ppc_pval         = numeric(),
  ppc_pass         = logical(),
  intercept_mean   = numeric(),
  intercept_lower  = numeric(),
  intercept_upper  = numeric(),
  run_time_min     = numeric(),
  stringsAsFactors = FALSE
)

####################
#   Main Modeling Loop
####################

cat("========================================\n")
cat("STARTING SENSITIVITY ANALYSIS MODELS\n")
cat("========================================\n\n")

batch_start <- Sys.time()

for (i in 1:nrow(viable)) {
  
  sp          <- viable$species[i]
  sid         <- viable$site_id[i]
  viable_covs <- strsplit(viable$viable_covariates[i], ", ")[[1]]
  
  cat("\n[", i, "/", nrow(viable), "]", sp, "at site", sid, "\n")
  
  tryCatch({
    
    model_start <- Sys.time()
    
    ####################
    #   Get Site and Create Grid
    ####################
    
    site <- sites_utm %>% filter(site_id == sid)
    
    if (nrow(site) == 0) {
      cat("  Site not found - skipping\n")
      next
    }
    
    # Create hex grid over full site
    hex_grid     <- st_make_grid(site, cellsize = cell_size, square = FALSE)
    hex_sf       <- st_sf(geometry = hex_grid)
    hex_filtered <- hex_sf[lengths(st_intersects(hex_sf, site)) > 0, ]
    hex_filtered$cell_id <- 1:nrow(hex_filtered)
    
    n_cells <- nrow(hex_filtered)
    cat("  Created", n_cells, "cells\n")
    
    ####################
    #   Assign JITTERED Observations to Cells
    ####################
    
    # Species observations at this site — using jittered coordinates
    sp_obs <- inat %>% filter(common_name == sp, site_id == sid)
    
    sp_sf <- st_as_sf(sp_obs,
                      coords = c("longitude_j", "latitude_j"),
                      crs = 4326)
    sp_sf <- st_transform(sp_sf, crs = 32614)
    
    obs_in_cells  <- st_intersects(sp_sf, hex_filtered)
    sp_sf$cell_id <- sapply(obs_in_cells, function(x) {
      if (length(x) > 0) return(x[1]) else return(NA)
    })
    sp_sf <- sp_sf %>% filter(!is.na(cell_id))
    
    # All observations at site (effort) — also jittered
    all_site_obs <- inat %>% filter(site_id == sid)
    all_sf <- st_as_sf(all_site_obs,
                       coords = c("longitude_j", "latitude_j"),
                       crs = 4326)
    all_sf <- st_transform(all_sf, crs = 32614)
    
    all_in_cells   <- st_intersects(all_sf, hex_filtered)
    all_sf$cell_id <- sapply(all_in_cells, function(x) {
      if (length(x) > 0) return(x[1]) else return(NA)
    })
    
    cells_with_effort <- unique(all_sf$cell_id[!is.na(all_sf$cell_id)])
    n_cells_effort    <- length(cells_with_effort)
    
    cat("  Cells with effort:", n_cells_effort, "of", n_cells, "\n")
    
    ####################
    #   Extract Covariates (ALL cells)
    ####################
    
    cov_data <- data.frame(cell_id = hex_filtered$cell_id)
    
    for (j in 1:nrow(hex_filtered)) {
      cell    <- hex_filtered[j, ]
      lc_vals <- terra::extract(landcover, vect(cell), df = TRUE)
      
      if (nrow(lc_vals) == 0 || all(is.na(lc_vals[, 2]))) {
        cov_data$trees_pct[j]       <- NA
        cov_data$grass_pct[j]       <- NA
        cov_data$shrub_pct[j]       <- NA
        cov_data$flooded_veg_pct[j] <- NA
        cov_data$crops_pct[j]       <- NA
        cov_data$water_pct[j]       <- NA
        cov_data$habitat_div[j]     <- NA
      } else {
        lc_col   <- names(lc_vals)[2]
        total_px <- sum(!is.na(lc_vals[[lc_col]]))
        
        lc_sum <- lc_vals %>%
          filter(!is.na(.data[[lc_col]])) %>%
          count(.data[[lc_col]]) %>%
          mutate(pct = n / total_px * 100)
        names(lc_sum)[1] <- "class"
        
        get_pct <- function(cls) {
          ifelse(cls %in% lc_sum$class, lc_sum$pct[lc_sum$class == cls], 0)
        }
        
        cov_data$water_pct[j]       <- get_pct(0)
        cov_data$trees_pct[j]       <- get_pct(1)
        cov_data$grass_pct[j]       <- get_pct(2)
        cov_data$flooded_veg_pct[j] <- get_pct(3)
        cov_data$crops_pct[j]       <- get_pct(4)
        cov_data$shrub_pct[j]       <- get_pct(5)
        
        props               <- lc_sum$pct / 100
        props               <- props[props > 0]
        cov_data$habitat_div[j] <- -sum(props * log(props), na.rm = TRUE)
      }
    }
    
    names(cov_data) <- gsub("habitat_div$", "habitat_diversity", names(cov_data))
    names(cov_data) <- gsub("flooded_veg_pct", "flooded_veg", names(cov_data))
    names(cov_data) <- gsub("_pct", "", names(cov_data))
    
    cov_name_map <- c(
      "trees"       = "trees",
      "grass"       = "grass",
      "shrub"       = "shrub",
      "flooded_veg" = "flooded_veg",
      "crops"       = "crops",
      "water"       = "water",
      "habitat_div" = "habitat_diversity"
    )
    
    covs_to_use <- cov_name_map[viable_covs]
    covs_to_use <- covs_to_use[covs_to_use %in% names(cov_data)]
    
    occ_covs <- cov_data[, covs_to_use, drop = FALSE]
    
    var_check <- sapply(occ_covs, var, na.rm = TRUE)
    keep_covs <- names(var_check)[!is.na(var_check) & var_check > 0.01]
    
    if (length(keep_covs) < 2) {
      cat("  Insufficient covariate variance - skipping\n")
      next
    }
    
    occ_covs        <- occ_covs[, keep_covs, drop = FALSE]
    occ_covs_scaled <- as.data.frame(scale(occ_covs))
    
    cat("  Using covariates:", paste(keep_covs, collapse = ", "), "\n")
    
    ####################
    #   Build Detection Matrix (ALL cells)
    ####################
    
    n_model_cells <- n_cells
    
    detection_matrix <- array(0,
                              dim = c(n_model_cells, n_years, n_months),
                              dimnames = list(
                                cell  = 1:n_model_cells,
                                year  = years,
                                month = months
                              ))
    
    sp_df <- st_drop_geometry(sp_sf)
    
    for (k in 1:nrow(sp_df)) {
      obs       <- sp_df[k, ]
      cell_idx  <- obs$cell_id
      year_idx  <- which(years == obs$year)
      month_idx <- obs$month
      
      if (is.na(cell_idx) || length(year_idx) == 0 || is.na(month_idx)) next
      detection_matrix[cell_idx, year_idx, month_idx] <- 1
    }
    
    cells_detected   <- apply(detection_matrix, 1, function(x) any(x == 1))
    n_cells_detected <- sum(cells_detected)
    naive_occ        <- n_cells_detected / n_model_cells
    
    cat("  Cells detected:", n_cells_detected, "| Naive occ:", round(naive_occ, 3), "\n")
    
    ####################
    #   Prepare Detection Covariates
    ####################
    
    year_array <- array(NA, dim = c(n_model_cells, n_years, n_months))
    for (y in 1:n_years) year_array[, y, ] <- years[y]
    year_array_scaled <- (year_array - mean(years)) / sd(years)
    
    peak_months       <- c(1, 2, 3, 4, 11, 12)
    peak_season_array <- array(NA, dim = c(n_model_cells, n_years, n_months))
    for (m in 1:n_months) {
      peak_season_array[, , m] <- ifelse(m %in% peak_months, 1, 0)
    }
    
    ####################
    #   Fit Model
    ####################
    
    data_list <- list(
      y        = detection_matrix,
      occ.covs = occ_covs_scaled,
      det.covs = list(
        year        = year_array_scaled,
        peak_season = peak_season_array
      )
    )
    
    occ_formula <- as.formula(paste("~", paste(keep_covs, collapse = " + ")))
    
    cat("  Fitting model... ")
    
    model <- tPGOcc(
      occ.formula   = occ_formula,
      det.formula   = ~ year + peak_season,
      data          = data_list,
      n.batch       = model_settings$n_batch,
      batch.length  = model_settings$batch_length,
      n.burn        = model_settings$n_burn,
      n.thin        = model_settings$n_thin,
      n.chains      = model_settings$n_chains,
      n.omp.threads = model_settings$n_omp_threads,
      verbose       = FALSE
    )
    
    ####################
    #   Assess Convergence & Fit
    ####################
    
    model_summary <- summary(model)
    rhat_vals     <- model_summary$beta[, "Rhat"]
    converged     <- all(rhat_vals < 1.1, na.rm = TRUE)
    
    ppc      <- ppcOcc(model, fit.stat = 'freeman-tukey', group = 1)
    ppc_pval <- mean(ppc$fit.y.rep > ppc$fit.y)
    ppc_pass <- ppc_pval > 0.1 & ppc_pval < 0.9
    
    run_time <- as.numeric(difftime(Sys.time(), model_start, units = "mins"))
    
    cat("done (", round(run_time, 2), "min) | p-value:", round(ppc_pval, 3))
    if (ppc_pass) cat(" ✓\n") else cat(" ⚠\n")
    
    ####################
    #   Store Results
    ####################
    
    result_name <- paste0(gsub(" ", "_", sp), "_site", sid)
    
    all_results[[result_name]] <- list(
      species             = sp,
      site_id             = sid,
      model               = model,
      model_summary       = model_summary,
      ppc_pval            = ppc_pval,
      detection_matrix    = detection_matrix,
      covariates          = occ_covs_scaled,
      covariates_unscaled = occ_covs,
      covariate_names     = keep_covs,
      n_cells             = n_model_cells,
      n_cells_effort      = n_cells_effort,
      n_cells_detected    = n_cells_detected,
      naive_occ           = naive_occ,
      converged           = converged,
      run_time            = run_time,
      # Sensitivity analysis metadata
      sensitivity = list(
        jitter_seed    = jitter_seed,
        max_jitter_m   = max_jitter_m,
        cell_size_m    = cell_size
      )
    )
    
    saveRDS(all_results[[result_name]],
            file.path(sensitivity_dir, paste0(result_name, ".rds")))
    
    intercept_stats <- model_summary$beta["(Intercept)", ]
    
    summary_row <- data.frame(
      species          = sp,
      site_id          = sid,
      n_obs            = viable$n_obs[i],
      n_cells          = n_model_cells,
      n_cells_effort   = n_cells_effort,
      n_cells_detected = n_cells_detected,
      naive_occ        = round(naive_occ, 3),
      n_covariates     = length(keep_covs),
      covariates_used  = paste(keep_covs, collapse = ", "),
      model_converged  = converged,
      ppc_pval         = round(ppc_pval, 3),
      ppc_pass         = ppc_pass,
      intercept_mean   = round(intercept_stats["Mean"], 3),
      intercept_lower  = round(intercept_stats["2.5%"], 3),
      intercept_upper  = round(intercept_stats["97.5%"], 3),
      run_time_min     = round(run_time, 2)
    )
    
    summary_results <- rbind(summary_results, summary_row)
    
    rm(model, detection_matrix, data_list, sp_sf, all_sf)
    gc(verbose = FALSE)
    
  }, error = function(e) {
    cat("  ERROR:", e$message, "\n")
  })
}

####################
#   Final Summary
####################

batch_end  <- Sys.time()
total_time <- difftime(batch_end, batch_start, units = "mins")

cat("\n========================================\n")
cat("SENSITIVITY ANALYSIS COMPLETE\n")
cat("========================================\n\n")

cat("Total models attempted:", nrow(viable), "\n")
cat("Models completed:", nrow(summary_results), "\n")
cat("Models converged:", sum(summary_results$model_converged), "\n")
cat("Models passing PPC:", sum(summary_results$ppc_pass), "\n")
cat("Total run time:", round(total_time, 1), "minutes\n\n")

write.csv(summary_results,
          file.path(sensitivity_dir, "sensitivity_model_summary.csv"),
          row.names = FALSE)
cat("Saved: sensitivity_model_summary.csv\n\n")

saveRDS(all_results,
        file.path(sensitivity_dir, paste0("all_results_sensitivity_", Sys.Date(), ".rds")))
cat("Saved: all_results_sensitivity_", Sys.Date(), ".rds\n")

####################
#   Extract Coefficients
####################

cat("\n=== EXTRACTING COEFFICIENTS ===\n\n")

model_files <- list.files(sensitivity_dir, pattern = "\\.rds$", full.names = TRUE)
model_files <- model_files[!grepl("all_results", model_files)]

cat("Found", length(model_files), "model files\n\n")

all_coefs <- data.frame()

for (f in model_files) {
  tryCatch({
    result <- readRDS(f)
    sp     <- result$species
    sid    <- result$site_id
    model  <- result$model
    
    beta_samples <- as.matrix(model$beta.samples)
    
    for (p in colnames(beta_samples)) {
      vals <- beta_samples[, p]
      coef_row <- data.frame(
        species   = sp,
        site_id   = sid,
        parameter = p,
        Mean      = mean(vals),
        SD        = sd(vals),
        lower     = quantile(vals, 0.025),
        upper     = quantile(vals, 0.975),
        Rhat      = ifelse(!is.null(model$rhat$beta[p]), model$rhat$beta[p], NA)
      )
      all_coefs <- rbind(all_coefs, coef_row)
    }
    cat("Extracted:", sp, "site", sid, "\n")
  }, error = function(e) {
    cat("Error with", f, ":", e$message, "\n")
  })
}

cat("\n=== SUMMARY ===\n")
cat("Total coefficients:", nrow(all_coefs), "\n")
cat("Species:", length(unique(all_coefs$species)), "\n")
cat("Site-species combinations:", nrow(distinct(all_coefs, species, site_id)), "\n")

write.csv(all_coefs,
          file.path(sensitivity_dir, "sensitivity_within_site_results.csv"),
          row.names = FALSE)
cat("\nSaved: sensitivity_within_site_results.csv\n")

cat("\n=== COMPLETE ===\n")
cat("Compare sensitivity_within_site_results.csv with\n")
cat("within_site_models/within_site_results.csv\n")
cat("to assess robustness to positional uncertainty.\n")