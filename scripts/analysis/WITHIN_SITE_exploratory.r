##############################
#
# Within-Site Fine-Scale Occupancy Analysis
# Black-bellied Whistling-Duck at Site 24
# Ian Becker
# January 2026
#
##############################

# THIS IS EXPLORATORY WITHIN SITE TEST

library(sf)
library(terra)
library(dplyr)
library(lubridate)
library(spOccupancy)

####################
#   Settings
####################

target_site_id <- 191
target_species <- "Black-bellied Whistling-Duck"
cell_size <- 250  # meters

cat("=== WITHIN-SITE OCCUPANCY ANALYSIS ===\n")
cat("Species:", target_species, "\n")
cat("Site:", target_site_id, "\n")
cat("Grid cell size:", cell_size, "m hexagons\n\n")

####################
#   Part 1: Load Site and Create Grid
####################

cat("=== PART 1: CREATING HEXAGONAL GRID ===\n\n")

# Load sites shapefile
sites <- st_read("lrgv_green_spaces_detection_filtered")

# Extract target site
target_site <- sites %>% filter(site_id == target_site_id)

cat("Site", target_site_id, "area:", round(as.numeric(st_area(target_site)) / 10000, 2), "hectares\n")

# Transform to UTM for accurate grid creation
target_site_utm <- st_transform(target_site, crs = 32614)

# Create hexagonal grid
hex_grid <- st_make_grid(target_site_utm, 
                         cellsize = cell_size, 
                         square = FALSE)

# Convert to sf and clip to site boundary
hex_grid_sf <- st_sf(geometry = hex_grid)
hex_grid_filtered <- hex_grid_sf[lengths(st_intersects(hex_grid_sf, target_site_utm)) > 0, ]

# Add cell IDs
hex_grid_filtered$cell_id <- 1:nrow(hex_grid_filtered)

# Calculate cell areas (some edge cells will be smaller)
hex_grid_filtered$area_m2 <- as.numeric(st_area(hex_grid_filtered))
hex_grid_filtered$area_ha <- hex_grid_filtered$area_m2 / 10000

n_cells <- nrow(hex_grid_filtered)
cat("Created", n_cells, "hexagonal grid cells\n")

# Save grid
st_write(hex_grid_filtered, 
         paste0("site_", target_site_id, "_hex_grid_", cell_size, "m.shp"),
         delete_dsn = TRUE)

####################
#   Part 2: Extract Land Cover for Each Cell
####################

cat("=== PART 2: EXTRACTING LAND COVER COVARIATES ===\n\n")

# Load land cover raster
landcover <- rast("lrgv_dynamic_world_cover.tif")
landcover_utm <- project(landcover, "EPSG:32614", method = "near")

# Initialize covariate dataframe
cell_covariates <- data.frame(
  cell_id = integer(),
  area_ha = numeric(),
  log_area = numeric(),
  trees_pct = numeric(),
  grass_pct = numeric(),
  shrub_pct = numeric(),
  flooded_veg_pct = numeric(),
  crops_pct = numeric(),
  water_pct = numeric(),
  habitat_diversity = numeric(),
  stringsAsFactors = FALSE
)

# Extract for each cell
cat("Extracting land cover for", n_cells, "cells...\n")

for (i in 1:n_cells) {
  
  if (i %% 10 == 0) cat("  Cell", i, "of", n_cells, "\n")
  
  cell <- hex_grid_filtered[i, ]
  
  # Extract land cover values
  lc_values <- terra::extract(landcover_utm, vect(cell), df = TRUE)
  
  if (nrow(lc_values) == 0 || all(is.na(lc_values[,2]))) {
    cell_result <- data.frame(
      cell_id = cell$cell_id,
      area_ha = cell$area_ha,
      log_area = log(cell$area_ha + 0.01),
      trees_pct = 0, grass_pct = 0, shrub_pct = 0,
      flooded_veg_pct = 0, crops_pct = 0, habitat_diversity = 0, water_pct = 0
    )
  } else {
    lc_column <- names(lc_values)[2]
    total_pixels <- sum(!is.na(lc_values[[lc_column]]))
    
    lc_summary <- lc_values %>%
      filter(!is.na(.data[[lc_column]])) %>%
      count(.data[[lc_column]]) %>%
      mutate(percentage = n / total_pixels * 100)
    
    names(lc_summary)[1] <- "class"
    
    # Dynamic World classes: 1=Trees, 2=Grass, 3=Flooded, 4=Crops, 5=Shrub
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
    # In the land cover extraction loop, add:
    water_pct <- ifelse(0 %in% lc_summary$class, 
                        lc_summary$percentage[lc_summary$class == 0], 0)
    
    # Shannon diversity
    proportions <- lc_summary$percentage / 100
    proportions <- proportions[proportions > 0]
    habitat_diversity <- -sum(proportions * log(proportions), na.rm = TRUE)
    
    cell_result <- data.frame(
      cell_id = cell$cell_id,
      area_ha = cell$area_ha,
      log_area = log(cell$area_ha + 0.01),
      trees_pct = trees_pct,
      grass_pct = grass_pct,
      shrub_pct = shrub_pct,
      flooded_veg_pct = flooded_veg_pct,
      crops_pct = crops_pct,
      water_pct = water_pct,
      habitat_diversity = habitat_diversity
    )
  }
  
  cell_covariates <- rbind(cell_covariates, cell_result)
}

cat("\nCovariate summary:\n")
print(summary(cell_covariates[, -1]))

# Save covariates
saveRDS(cell_covariates, paste0("site_", target_site_id, "_cell_covariates.rds"))

####################
#   Part 3: Assign Observations to Grid Cells
####################

cat("\n=== PART 3: ASSIGNING OBSERVATIONS TO CELLS ===\n\n")

# Load iNat data with site assignments
inat_with_sites <- read.csv("inat_observations_with_sites.csv")

# Filter to target species and site
species_data <- inat_with_sites %>%
  filter(common_name == target_species,
                  site_id == 67)

cat("Total observations:", nrow(species_data), "\n")

# Convert to spatial
species_sf <- st_as_sf(species_data,
                       coords = c("longitude", "latitude"),
                       crs = 4326)
species_sf <- st_transform(species_sf, crs = 32614)

# Assign to grid cells
obs_in_cells <- st_intersects(species_sf, hex_grid_filtered)
species_sf$cell_id <- sapply(obs_in_cells, function(x) {
  if (length(x) > 0) return(hex_grid_filtered$cell_id[x[1]])
  else return(NA)
})

# Remove observations that fall outside cells (edge cases)
species_sf <- species_sf %>% filter(!is.na(cell_id))

cat("Observations assigned to cells:", nrow(species_sf), "\n")
cat("Cells with observations:", n_distinct(species_sf$cell_id), "of", n_cells, "\n\n")

# Add temporal info
species_sf$date <- as.Date(species_sf$observed_on)
species_sf$year <- year(species_sf$date)
species_sf$month <- month(species_sf$date)

# Filter to study period
species_sf <- species_sf %>%
  filter(year >= 2015, year <= 2025)

cat("After filtering to 2015-2025:", nrow(species_sf), "observations\n")

####################
#   Part 4: Build Detection Matrix
####################

cat("\n=== PART 4: BUILDING DETECTION MATRIX ===\n\n")

years <- 2015:2025
n_years <- length(years)
months <- 1:12
n_months <- length(months)

# Initialize 3D array: cells × years × months
detection_matrix <- array(0, 
                          dim = c(n_cells, n_years, n_months),
                          dimnames = list(
                            cell = 1:n_cells,
                            year = years,
                            month = months
                          ))

# Fill in detections
species_df <- st_drop_geometry(species_sf)

for (i in 1:nrow(species_df)) {
  obs <- species_df[i, ]
  
  cell_idx <- obs$cell_id
  year_idx <- which(years == obs$year)
  month_idx <- obs$month
  
  if (is.na(cell_idx) || cell_idx < 1 || cell_idx > n_cells) next
  if (length(year_idx) == 0) next
  if (is.na(month_idx) || month_idx < 1 || month_idx > 12) next
  
  detection_matrix[cell_idx, year_idx, month_idx] <- 1
}

# Summary
cells_with_detections <- apply(detection_matrix, 1, function(x) any(x == 1))
n_cells_detected <- sum(cells_with_detections)

cat("Detection matrix dimensions:", dim(detection_matrix), "\n")
cat("Cells with detections:", n_cells_detected, "of", n_cells, 
    "(", round(n_cells_detected/n_cells*100, 1), "%)\n")
cat("Naive occupancy:", round(n_cells_detected/n_cells, 3), "\n\n")

####################
#   Part 6: Fit Occupancy Model
####################

cat("=== PART 6: FITTING OCCUPANCY MODEL ===\n\n")

# Remove log_area since cells are ~same size
occ_covs <- cell_covariates %>%
  select(trees_pct, grass_pct, shrub_pct, flooded_veg_pct, 
         crops_pct, habitat_diversity, water_pct)

# Check variance BEFORE scaling
var_check <- sapply(occ_covs, var, na.rm = TRUE)
cat("Covariate variance:\n")
print(round(var_check, 4))

# Remove covariates with zero or near-zero variance
keep_covs <- names(var_check)[!is.na(var_check) & var_check > 0.01]
occ_covs <- occ_covs[, keep_covs, drop = FALSE]

# NOW scale
occ_covs_scaled <- as.data.frame(scale(occ_covs))

cat("\nUsing covariates:", paste(keep_covs, collapse = ", "), "\n\n")

# Remove covariates with zero variance
keep_covs <- names(var_check)[var_check > 0.01]
occ_covs_scaled <- occ_covs_scaled[, keep_covs, drop = FALSE]
cat("\nUsing covariates:", paste(keep_covs, collapse = ", "), "\n\n")

# Prepare detection covariates
year_array <- array(NA, dim = c(n_cells, n_years, n_months))
for (j in 1:n_years) {
  year_array[, j, ] <- years[j]
}
year_array_scaled <- (year_array - mean(years)) / sd(years)

peak_months <- c(1, 2, 3, 4, 11, 12)
peak_season_array <- array(NA, dim = c(n_cells, n_years, n_months))
for (m in 1:n_months) {
  peak_season_array[, , m] <- ifelse(m %in% peak_months, 1, 0)
}

# Create data list
data_list <- list(
  y = detection_matrix,
  occ.covs = occ_covs_scaled,
  det.covs = list(
    year = year_array_scaled,
    peak_season = peak_season_array
  )
)

# Model settings
model_settings <- list(
  n_batch = 1000,
  batch_length = 10,
  n_burn = 5000,
  n_thin = 5,
  n_chains = 3,
  n_omp_threads = 1
)

# Build formula dynamically based on available covariates
occ_formula <- as.formula(paste("~", paste(keep_covs, collapse = " + ")))
cat("Occupancy formula:", deparse(occ_formula), "\n")
cat("Detection formula: ~ year + peak_season\n\n")

# Fit model
cat("Fitting model...\n")
model_start <- Sys.time()

model <- tPGOcc(
  occ.formula = occ_formula,
  det.formula = ~ year + peak_season,
  data = data_list,
  n.batch = model_settings$n_batch,
  batch.length = model_settings$batch_length,
  n.burn = model_settings$n_burn,
  n.thin = model_settings$n_thin,
  n.chains = model_settings$n_chains,
  n.omp.threads = model_settings$n_omp_threads,
  verbose = TRUE
)

run_time <- difftime(Sys.time(), model_start, units = "mins")
cat("\nModel fitting complete in", round(run_time, 2), "minutes\n\n")

####################
#   Part 7: Model Results
####################

cat("=== PART 7: MODEL RESULTS ===\n\n")

print(summary(model))

# Posterior predictive check
ppc <- ppcOcc(model, fit.stat = 'freeman-tukey', group = 1)
ppc_pval <- mean(ppc$fit.y.rep > ppc$fit.y)
cat("\nBayesian p-value:", round(ppc_pval, 3))
if (ppc_pval > 0.1 & ppc_pval < 0.9) cat(" ✓\n") else cat(" ⚠\n")

####################
#   Save Results
####################

cat("\n=== SAVING RESULTS ===\n\n")

results <- list(
  species = target_species,
  site_id = target_site_id,
  cell_size = cell_size,
  n_cells_total = n_cells,
  n_cells_with_effort = n_cells_filtered,
  model = model,
  detection_matrix = detection_matrix_filtered,
  cell_covariates = cell_covariates_filtered,
  hex_grid = hex_grid_filtered,
  ppc_pval = ppc_pval,
  run_time = run_time
)

saveRDS(results, paste0("within_site_model_site", target_site_id, "_", 
                        gsub(" ", "_", target_species), "_", cell_size, "m.rds"))

cat("Results saved!\n")
cat("=== ANALYSIS COMPLETE ===\n")