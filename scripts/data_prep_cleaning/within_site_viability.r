##############################
#
# Within-Site Model Viability Assessment
# Ian Becker
# January 2026
#
##############################

# This script checks available within site species-site combos

library(sf)
library(terra)
library(dplyr)
library(tidyr)

####################
#   Settings
####################

cell_size <- 200  # meters
min_obs <- 10     # minimum observations per species-site

cat("=== WITHIN-SITE VIABILITY ASSESSMENT ===\n")
cat("Cell size:", cell_size, "m\n")
cat("Minimum observations:", min_obs, "\n\n")

####################
#   Load Data
####################

cat("Loading data...\n")

# Sites
sites <- st_read("lrgv_green_spaces_detection_filtered", quiet = TRUE)
sites_utm <- st_transform(sites, crs = 32614)
cat("Loaded", nrow(sites), "sites\n")

# iNat observations
inat <- read.csv("inat_observations_with_sites.csv")
cat("Loaded", nrow(inat), "observations\n")

# Land cover
landcover <- rast("lrgv_dynamic_world_cover.tif")
cat("Loaded land cover raster\n\n")

# Study species
waterbirds <- c("Black-bellied Whistling-Duck", "Green Heron", 
                "Neotropic Cormorant", "Great Blue Heron", "Great Egret")
generalists <- c("Great-tailed Grackle", "Northern Mockingbird", 
                 "White-winged Dove", "Inca Dove", "Mourning Dove")
specialists <- c("Great Kiskadee", "Tropical Kingbird", "Green Jay", 
                 "Long-billed Thrasher", "Golden-fronted Woodpecker")
study_species <- c(waterbirds, generalists, specialists)

####################
#   Identify Candidate Site-Species Combinations
####################

cat("=== IDENTIFYING CANDIDATES ===\n\n")

# Count observations per species-site
obs_counts <- inat %>%
  filter(common_name %in% study_species) %>%
  group_by(common_name, site_id) %>%
  summarise(n_obs = n(), .groups = "drop") %>%
  filter(n_obs >= min_obs)

cat("Site-species combinations with ≥", min_obs, "obs:", nrow(obs_counts), "\n\n")

####################
#   Assess Each Combination
####################

cat("=== ASSESSING VIABILITY ===\n\n")

# Initialize results
viability_results <- data.frame(
  species = character(),
  guild = character(),
  site_id = integer(),
  site_area_ha = numeric(),
  n_obs = integer(),
  n_cells_total = integer(),
  n_cells_with_effort = integer(),
  n_cells_detected = integer(),
  naive_occupancy = numeric(),
  trees_var = numeric(),
  grass_var = numeric(),
  shrub_var = numeric(),
  flooded_veg_var = numeric(),
  crops_var = numeric(),
  water_var = numeric(),
  habitat_div_var = numeric(),
  n_viable_covariates = integer(),
  viable_covariates = character(),
  stringsAsFactors = FALSE
)

# Loop through each candidate
for (i in 1:nrow(obs_counts)) {
  
  sp <- obs_counts$common_name[i]
  sid <- obs_counts$site_id[i]
  n_obs <- obs_counts$n_obs[i]
  
  cat("Processing:", sp, "at site", sid, "(", i, "/", nrow(obs_counts), ")\n")
  
  # Determine guild
  guild <- case_when(
    sp %in% waterbirds ~ "waterbird",
    sp %in% generalists ~ "generalist",
    sp %in% specialists ~ "specialist"
  )
  
  # Get site
  site <- sites_utm %>% filter(site_id == sid)
  
  if (nrow(site) == 0) {
    cat("  Site not found - skipping\n")
    next
  }
  
  site_area_ha <- as.numeric(st_area(site)) / 10000
  
  # Create hex grid (full hexagons only)
  hex_grid <- st_make_grid(site, cellsize = cell_size, square = FALSE)
  hex_sf <- st_sf(geometry = hex_grid)
  hex_filtered <- hex_sf[lengths(st_intersects(hex_sf, site)) > 0, ]
  hex_filtered$cell_id <- 1:nrow(hex_filtered)
  
  n_cells <- nrow(hex_filtered)
  
  if (n_cells < 3) {
    cat("  Only", n_cells, "cells - skipping\n")
    next
  }
  
  # Get species observations at this site
  sp_obs <- inat %>%
    filter(common_name == sp, site_id == sid)
  
  # Convert to spatial
  sp_sf <- st_as_sf(sp_obs, coords = c("longitude", "latitude"), crs = 4326)
  sp_sf <- st_transform(sp_sf, crs = 32614)
  
  # Assign to cells
  obs_in_cells <- st_intersects(sp_sf, hex_filtered)
  sp_sf$cell_id <- sapply(obs_in_cells, function(x) {
    if (length(x) > 0) return(x[1]) else return(NA)
  })
  
  cells_detected <- unique(sp_sf$cell_id[!is.na(sp_sf$cell_id)])
  n_cells_detected <- length(cells_detected)
  
  # Get ALL observations at site for effort assessment
  all_site_obs <- inat %>% filter(site_id == sid)
  all_sf <- st_as_sf(all_site_obs, coords = c("longitude", "latitude"), crs = 4326)
  all_sf <- st_transform(all_sf, crs = 32614)
  
  all_in_cells <- st_intersects(all_sf, hex_filtered)
  all_sf$cell_id <- sapply(all_in_cells, function(x) {
    if (length(x) > 0) return(x[1]) else return(NA)
  })
  
  cells_with_effort <- unique(all_sf$cell_id[!is.na(all_sf$cell_id)])
  n_cells_effort <- length(cells_with_effort)
  
  naive_occ <- n_cells_detected / n_cells_effort
  
  # Extract land cover covariates for cells with effort
  hex_effort <- hex_filtered %>% filter(cell_id %in% cells_with_effort)
  
  # Quick covariate extraction
  cov_data <- data.frame(cell_id = hex_effort$cell_id)
  
  for (j in 1:nrow(hex_effort)) {
    cell <- hex_effort[j, ]
    lc_vals <- terra::extract(landcover, vect(cell), df = TRUE)
    
    if (nrow(lc_vals) == 0 || all(is.na(lc_vals[,2]))) {
      cov_data$trees_pct[j] <- NA
      cov_data$grass_pct[j] <- NA
      cov_data$shrub_pct[j] <- NA
      cov_data$flooded_veg_pct[j] <- NA
      cov_data$crops_pct[j] <- NA
      cov_data$water_pct[j] <- NA
      cov_data$habitat_div[j] <- NA
    } else {
      lc_col <- names(lc_vals)[2]
      total_px <- sum(!is.na(lc_vals[[lc_col]]))
      
      lc_sum <- lc_vals %>%
        filter(!is.na(.data[[lc_col]])) %>%
        count(.data[[lc_col]]) %>%
        mutate(pct = n / total_px * 100)
      names(lc_sum)[1] <- "class"
      
      get_pct <- function(cls) {
        ifelse(cls %in% lc_sum$class, lc_sum$pct[lc_sum$class == cls], 0)
      }
      
      cov_data$water_pct[j] <- get_pct(0)
      cov_data$trees_pct[j] <- get_pct(1)
      cov_data$grass_pct[j] <- get_pct(2)
      cov_data$flooded_veg_pct[j] <- get_pct(3)
      cov_data$crops_pct[j] <- get_pct(4)
      cov_data$shrub_pct[j] <- get_pct(5)
      
      # Habitat diversity
      props <- lc_sum$pct / 100
      props <- props[props > 0]
      cov_data$habitat_div[j] <- -sum(props * log(props), na.rm = TRUE)
    }
  }
  
  # Calculate variances
  trees_var <- var(cov_data$trees_pct, na.rm = TRUE)
  grass_var <- var(cov_data$grass_pct, na.rm = TRUE)
  shrub_var <- var(cov_data$shrub_pct, na.rm = TRUE)
  flooded_var <- var(cov_data$flooded_veg_pct, na.rm = TRUE)
  crops_var <- var(cov_data$crops_pct, na.rm = TRUE)
  water_var <- var(cov_data$water_pct, na.rm = TRUE)
  hab_div_var <- var(cov_data$habitat_div, na.rm = TRUE)
  
  # Count viable covariates (variance > 0.01)
  var_threshold <- 0.01
  vars <- c(trees = trees_var, grass = grass_var, shrub = shrub_var,
            flooded_veg = flooded_var, crops = crops_var, 
            water = water_var, habitat_div = hab_div_var)
  viable_covs <- names(vars)[!is.na(vars) & vars > var_threshold]
  n_viable <- length(viable_covs)
  
  # Store results
  result <- data.frame(
    species = sp,
    guild = guild,
    site_id = sid,
    site_area_ha = round(site_area_ha, 2),
    n_obs = n_obs,
    n_cells_total = n_cells,
    n_cells_with_effort = n_cells_effort,
    n_cells_detected = n_cells_detected,
    naive_occupancy = round(naive_occ, 3),
    trees_var = round(trees_var, 2),
    grass_var = round(grass_var, 2),
    shrub_var = round(shrub_var, 2),
    flooded_veg_var = round(flooded_var, 2),
    crops_var = round(crops_var, 2),
    water_var = round(water_var, 2),
    habitat_div_var = round(hab_div_var, 4),
    n_viable_covariates = n_viable,
    viable_covariates = paste(viable_covs, collapse = ", ")
  )
  
  viability_results <- rbind(viability_results, result)
}

####################
#   Summary
####################

cat("\n=== SUMMARY ===\n\n")

cat("Total site-species combinations assessed:", nrow(viability_results), "\n\n")

cat("By guild:\n")
print(table(viability_results$guild))

cat("\nCells with detections summary:\n")
print(summary(viability_results$n_cells_detected))

cat("\nViable covariates summary:\n")
print(summary(viability_results$n_viable_covariates))

cat("\nNaive occupancy summary:\n")
print(summary(viability_results$naive_occupancy))

####################
#   Save Results
####################

cat("\n=== SAVING ===\n\n")

write.csv(viability_results, "within_site_viability_assessment.csv", row.names = FALSE)
cat("Saved: within_site_viability_assessment.csv\n")

# Also save a quick summary by species
species_summary <- viability_results %>%
  group_by(species, guild) %>%
  summarise(
    n_viable_sites = n(),
    total_obs = sum(n_obs),
    mean_cells_detected = round(mean(n_cells_detected), 1),
    mean_naive_occ = round(mean(naive_occupancy), 3),
    .groups = "drop"
  ) %>%
  arrange(guild, desc(n_viable_sites))

print(species_summary)

write.csv(species_summary, "within_site_species_summary.csv", row.names = FALSE)
cat("\nSaved: within_site_species_summary.csv\n")

cat("\n=== COMPLETE ===\n")
