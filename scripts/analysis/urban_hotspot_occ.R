##############################
#
# LRGV Flycatcher Occupancy Modeling
# 12/20/2025
# Ian Becker
#
##############################

# Two-species occupancy analysis:
# - Pitangus sulphuratus (Great Kiskadee): 231 obs
# - Tyrannus melancholicus (Tropical Kingbird): 29 obs

library(spOccupancy)
library(dplyr)
library(tidyr)
library(sf)
library(terra)

####################
#   Load Data
####################

# Load filtered iNaturalist data
species_data <- readRDS("lrgv_flycatcher_inat_data.rds")

cat("=== SPECIES DATA LOADED ===\n")
for (sp in names(species_data)) {
  cat(sp, ":", nrow(species_data[[sp]]), "observations\n")
}
cat("\n")

# Load urban green space grid (25% vegetation threshold)
study_grid <- st_read("lrgv_osm_green_spaces_filtered")

cat("Study grid:", nrow(study_grid), "cells (25% vegetation threshold)\n\n")

# Transform to UTM if needed
if (st_crs(study_grid)$input != "EPSG:32614") {
  study_grid <- st_transform(study_grid, crs = 32614)
}

# Load landcover raster
landcover <- rast("lrgv_dynamic_world_cover.tif")
landcover <- project(landcover, "EPSG:32614", method = "near")

####################
#   Extract Habitat Covariates
####################

cat("Extracting habitat covariates...\n")

habitat_covariates <- data.frame(
  cell_id = integer(),
  trees_pct = numeric(),
  grass_pct = numeric(),
  shrub_pct = numeric(),
  flooded_veg_pct = numeric(),
  urban_pct = numeric(),
  crops_pct = numeric(),
  habitat_diversity = numeric(),
  stringsAsFactors = FALSE
)

for (i in 1:nrow(study_grid)) {
  
  if (i %% 100 == 0) cat("  Cell", i, "of", nrow(study_grid), "\n")
  
  grid_cell <- study_grid[i, ]
  lc_values <- terra::extract(landcover, vect(grid_cell), df = TRUE)
  
  if (nrow(lc_values) == 0) {
    cell_result <- data.frame(
      cell_id = i,
      trees_pct = 0, grass_pct = 0, shrub_pct = 0,
      flooded_veg_pct = 0, urban_pct = 0, crops_pct = 0,
      habitat_diversity = 0
    )
  } else {
    
    lc_column <- names(lc_values)[2]
    total_pixels <- nrow(lc_values)
    
    lc_summary <- lc_values %>%
      count(.data[[lc_column]]) %>%
      mutate(percentage = n / total_pixels * 100)
    
    names(lc_summary)[1] <- "class"
    
    trees_pct <- ifelse(1 %in% lc_summary$class, 
                        lc_summary$percentage[lc_summary$class == 1], 0)
    grass_pct <- ifelse(2 %in% lc_summary$class, 
                        lc_summary$percentage[lc_summary$class == 2], 0)
    flooded_veg_pct <- ifelse(3 %in% lc_summary$class,
                              lc_summary$percentage[lc_summary$class == 3], 0)
    crops_pct <- ifelse(4 %in% lc_summary$class, 
                        lc_summary$percentage[lc_summary$class == 4], 0)
    shrub_pct <- ifelse(5 %in% lc_summary$class,
                        lc_summary$percentage[lc_summary$class == 5], 0)
    urban_pct <- ifelse(6 %in% lc_summary$class, 
                        lc_summary$percentage[lc_summary$class == 6], 0)
    
    proportions <- lc_summary$percentage / 100
    proportions <- proportions[proportions > 0]
    habitat_diversity <- -sum(proportions * log(proportions), na.rm = TRUE)
    
    cell_result <- data.frame(
      cell_id = i,
      trees_pct = trees_pct,
      grass_pct = grass_pct,
      shrub_pct = shrub_pct,
      flooded_veg_pct = flooded_veg_pct,
      urban_pct = urban_pct,
      crops_pct = crops_pct,
      habitat_diversity = habitat_diversity
    )
  }
  
  habitat_covariates <- rbind(habitat_covariates, cell_result)
}

cat("\nHabitat covariates extracted for", nrow(habitat_covariates), "cells\n\n")

####################
#   Create Detection Histories
####################

cat("Creating detection histories for each species...\n\n")

detection_histories <- list()

for (species_name in names(species_data)) {
  
  cat("Processing", species_name, "...\n")
  
  current_species <- species_data[[species_name]]
  
  # Convert to sf and transform to UTM
  species_sf <- st_as_sf(current_species,
                         coords = c("longitude", "latitude"),
                         crs = 4326)
  species_sf <- st_transform(species_sf, crs = 32614)
  
  # Assign observations to grid cells
  obs_intersect <- st_intersects(species_sf, study_grid)
  species_sf$cell_id <- sapply(obs_intersect, function(x) ifelse(length(x) > 0, x[1], NA))
  
  # Remove observations not in grid
  species_sf <- species_sf[!is.na(species_sf$cell_id), ]
  
  cat("  ", nrow(species_sf), "observations assigned to grid cells\n")
  
  # Extract year from date
  date_col <- "observed_on"  # iNat standard column name
  species_sf$year <- as.numeric(format(as.Date(species_sf[[date_col]]), "%Y"))
  species_sf <- species_sf[!is.na(species_sf$year), ]
  
  # Define sampling occasions (years)
  sampling_occasions <- sort(unique(species_sf$year))
  n_occasions <- length(sampling_occasions)
  
  cat("  ", n_occasions, "sampling occasions (years):", 
      paste(sampling_occasions, collapse = ", "), "\n")
  
  # Create detection matrix (sites x occasions)
  n_sites <- nrow(study_grid)
  detection_matrix <- matrix(0, nrow = n_sites, ncol = n_occasions)
  
  # Fill in detections
  species_df <- as.data.frame(species_sf)
  for (i in 1:nrow(species_df)) {
    site_idx <- species_df$cell_id[i]
    occasion_idx <- which(sampling_occasions == species_df$year[i])
    if (length(occasion_idx) > 0) {
      detection_matrix[site_idx, occasion_idx] <- 1
    }
  }
  
  # Store results
  detection_histories[[species_name]] <- list(
    y = detection_matrix,
    species = species_name,
    n_sites = n_sites,
    n_occasions = n_occasions,
    occasions = sampling_occasions,
    n_detections = sum(detection_matrix),
    n_sites_occupied = sum(rowSums(detection_matrix) > 0)
  )
  
  cat("  Detection matrix:", n_sites, "sites x", n_occasions, "occasions\n")
  cat("  Sites with detections:", sum(rowSums(detection_matrix) > 0), "\n\n")
}

####################
#   Fit Occupancy Models
####################

cat("=== FITTING OCCUPANCY MODELS ===\n\n")

# Prepare covariates (standardized)
occ_covs <- habitat_covariates[, c("trees_pct", "grass_pct", "urban_pct",
                                   "crops_pct", "shrub_pct", "flooded_veg_pct",
                                   "habitat_diversity")]
occ_covs_scaled <- scale(occ_covs)

occupancy_results <- list()

for (species_name in names(detection_histories)) {
  
  cat("Fitting model for", species_name, "...\n")
  
  det_hist <- detection_histories[[species_name]]
  
  # Detection covariates (intercept only for now)
  det_covs <- list(
    intercept = matrix(1, nrow = det_hist$n_sites, ncol = det_hist$n_occasions)
  )
  
  # Fit model
  tryCatch({
    
    model <- PGOcc(
      occ.formula = ~ trees_pct + grass_pct + urban_pct + crops_pct +
        shrub_pct + flooded_veg_pct + habitat_diversity,
      det.formula = ~ 1,
      data = list(
        y = det_hist$y,
        occ.covs = data.frame(occ_covs_scaled),
        det.covs = det_covs
      ),
      n.samples = 10000,
      n.burn = 5000,
      n.thin = 5,
      n.chains = 3,
      n.omp.threads = 1,
      verbose = TRUE,
      n.report = 1000
    )
    
    # Store results
    occupancy_results[[species_name]] <- list(
      model = model,
      detection_history = det_hist,
      summary = summary(model),
      mean_occupancy = round(mean(model$psi.samples), 3),
      success = TRUE
    )
    
    cat("✓ Model converged. Mean occupancy:", 
        round(mean(model$psi.samples), 3), "\n\n")
    
  }, error = function(e) {
    cat("✗ Model failed:", e$message, "\n\n")
    occupancy_results[[species_name]] <- list(
      model = NULL,
      detection_history = det_hist,
      success = FALSE,
      error = e$message
    )
  })
}

####################
#   Model Summaries
####################

cat("\n=== MODEL RESULTS SUMMARY ===\n\n")

for (species_name in names(occupancy_results)) {
  
  result <- occupancy_results[[species_name]]
  
  cat(species_name, ":\n")
  
  if (result$success) {
    cat("  Status: ✓ Success\n")
    cat("  Mean occupancy:", result$mean_occupancy, "\n")
    cat("  Sites occupied:", result$detection_history$n_sites_occupied, 
        "of", result$detection_history$n_sites, "\n")
    
    # Print coefficient summary
    cat("\n  Occupancy coefficients:\n")
    beta_summary <- result$summary$beta
    print(round(beta_summary, 3))
    
  } else {
    cat("  Status: ✗ Failed\n")
    cat("  Error:", result$error, "\n")
  }
  
  cat("\n")
}

####################
#   Save Results
####################

save(occupancy_results, detection_histories, habitat_covariates, study_grid,
     file = paste0("lrgv_flycatcher_occupancy_", Sys.Date(), ".RData"))

cat("Results saved to: lrgv_flycatcher_occupancy_", Sys.Date(), ".RData\n")

# Save coefficient summaries
coef_summaries <- list()

for (species_name in names(occupancy_results)) {
  if (occupancy_results[[species_name]]$success) {
    
    beta_summary <- occupancy_results[[species_name]]$summary$beta
    
    coef_df <- data.frame(
      species = species_name,
      covariate = rownames(beta_summary),
      mean = beta_summary[, "Mean"],
      sd = beta_summary[, "SD"],
      ci_lower = beta_summary[, "2.5%"],
      ci_upper = beta_summary[, "97.5%"],
      significant = !(beta_summary[, "2.5%"] <= 0 & beta_summary[, "97.5%"] >= 0)
    )
    
    coef_summaries[[species_name]] <- coef_df
  }
}

if (length(coef_summaries) > 0) {
  all_coefs <- bind_rows(coef_summaries)
  write.csv(all_coefs, "lrgv_flycatcher_coefficients.csv", row.names = FALSE)
  cat("Coefficients saved to: lrgv_flycatcher_coefficients.csv\n")
}

cat("\n=== ANALYSIS COMPLETE ===\n")