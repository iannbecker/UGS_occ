##############################
#
# Occupancy Modelling Figures
# 8/29/2025
# Ian Becker 
#
##############################

library(ggplot2)
library(dplyr)
library(viridis)
library(sf)
library(terra)
library(dplyr)
library(igraph)
library(tidygraph)
library(ggraph)
library(ggpubr)
library(RColorBrewer)

load("multi_species_occupancy_results_2025-08-29.RData")
coef_data <- read.csv("occupancy_coefficients_summary.csv")
study_grid <- readRDS("study_grid.rds") 
landcover <- rast("urban_hotspots_landcover.tif")

###################
#
#  Heatmap
#
##################

# prepping heatmap data

coef_heatmap_data <- coef_data %>%
  filter(covariate != "(Intercept)") %>%
  mutate(
    covariate_clean = case_when(
      covariate == "forest_pct" ~ "Forest %",
      covariate == "grass_pct" ~ "Grassland %", 
      covariate == "urban_pct" ~ "Urban %",
      covariate == "crops_pct" ~ "Crops %",           
      covariate == "shrub_pct" ~ "Shrub/Scrub %",    
      covariate == "flooded_veg_pct" ~ "Wetlands %", 
      covariate == "habitat_diversity" ~ "Habitat Diversity",
      TRUE ~ covariate
    ),
    effect_strength = case_when(
      !significant ~ 0,
      significant & mean > 0 ~ mean,
      significant & mean < 0 ~ mean,
      TRUE ~ 0
    )
  )

# Multi-panel figure combining heatmap + individual forest plots

effect_range <- range(coef_heatmap_data$effect_strength[coef_heatmap_data$effect_strength != 0])
max_abs_effect <- max(abs(effect_range))

p1 <- ggplot(coef_heatmap_data, aes(x = covariate_clean, y = reorder(species, mean))) +
  geom_tile(aes(fill = effect_strength), color = "white", size = 0.5) +
  scale_fill_gradient2(
    low = "darkblue", 
    mid = "white", 
    high = "darkred",
    midpoint = 0,
    limits = c(-max_abs_effect, max_abs_effect),  # Symmetric limits
    name = "Effect Size"
  ) +
  theme_minimal() +
  labs(title = "Species-Habitat Associations Across Urban Green Spaces",
       x = "Habitat Variables", y = "Species") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p1

ggsave("species_habitat_heatmap.png", plot = p1, width = 10, height = 8, dpi = 300)

###################
#
#  conceptual map
#
##################

# Step 1: Load your single site shapefile

single_site <- st_read("utrgv_shapefile")
single_site <- st_transform(single_site, crs = 32614)

# Step 2: Create grid correctly
cell_size <- 150
site_grid_geom <- st_make_grid(single_site, cellsize = cell_size, square = FALSE)

# Step 3: Convert to proper sf object FIRST
site_grid <- st_sf(
  site_id = 1:length(site_grid_geom), 
  geometry = site_grid_geom
)

# Step 4: Set CRS explicitly
st_crs(site_grid) <- st_crs(single_site)

# Step 5: Now do intersection
grid_intersects <- st_intersects(site_grid, single_site, sparse = FALSE)[,1]

# Step 6: Filter to intersecting cells
site_grid <- site_grid[grid_intersects, ]

# Step 7: Renumber site_ids after filtering
site_grid$site_id <- 1:nrow(site_grid)

cat("Created", nrow(site_grid), "grid cells for single site\n")

# Step 4: Extract landcover data for this site's grid
# (You'll need to run the same landcover extraction as in your original script)

# Initialize the landcover variables dataframe

site_landcover_vars <- data.frame(
  site_id = integer(),
  forest_pct = numeric(), 
  grass_pct = numeric(), 
  crops_pct = numeric(),
  urban_pct = numeric(),
  shrub_pct = numeric(),
  flooded_veg_pct = numeric(),
  habitat_diversity = numeric(),
  stringsAsFactors = FALSE
)

# Extract landcover for each grid cell
for (i in 1:nrow(site_grid)) {
  
  if (i %% 10 == 0) cat("Processing cell", i, "of", nrow(site_grid), "\n")
  
  grid_cell <- site_grid[i,]
  lc_values <- terra::extract(landcover, vect(grid_cell), df = TRUE)
  
  if (nrow(lc_values) == 0) {
    # No data for this cell
    cell_results <- data.frame(
      site_id = grid_cell$site_id,
      forest_pct = 0, grass_pct = 0, crops_pct = 0, urban_pct = 0, 
      shrub_pct = 0, habitat_diversity = 0
    )
  } else {
    
    # Calculate land cover percentages
    lc_column <- names(lc_values)[2]
    total_pixels <- nrow(lc_values)
    
    lc_summary <- lc_values %>%
      count(.data[[lc_column]]) %>%
      mutate(percentage = n / total_pixels * 100)
    
    names(lc_summary)[1] <- "class"
    
    # Initialize variables
    forest_pct <- ifelse(1 %in% lc_summary$class, 
                         lc_summary$percentage[lc_summary$class == 1], 0)
    grass_pct <- ifelse(2 %in% lc_summary$class, 
                        lc_summary$percentage[lc_summary$class == 2], 0)
    crops_pct <- ifelse(4 %in% lc_summary$class, 
                        lc_summary$percentage[lc_summary$class == 4], 0)
    shrub_pct <- ifelse(5 %in% lc_summary$class,
                        lc_summary$percentage[lc_summary$class == 5], 0)
    urban_pct <- ifelse(6 %in% lc_summary$class, 
                        lc_summary$percentage[lc_summary$class == 6], 0)
    flooded_veg_pct <- ifelse(3 %in% lc_summary$class,
                              lc_summary$percentage[lc_summary$class == 3], 0)
    
    # Habitat diversity (Shannon index)
    proportions <- lc_summary$percentage / 100
    proportions <- proportions[proportions > 0]
    habitat_diversity <- -sum(proportions * log(proportions), na.rm = TRUE)
    
    cell_results <- data.frame(
      site_id = grid_cell$site_id,
      forest_pct = forest_pct,
      grass_pct = grass_pct, 
      crops_pct = crops_pct,
      urban_pct = urban_pct,
      shrub_pct = shrub_pct,
      flooded_veg_pct = flooded_veg_pct,
      habitat_diversity = habitat_diversity
    )
  }
  
  site_landcover_vars <- rbind(site_landcover_vars, cell_results)
}

cat("Extracted landcover data for", nrow(site_landcover_vars), "grid cells\n")

# Step 5: Apply your occupancy models to predict for this site
# Use your existing models to predict occupancy for each grid cell

# Initialize prediction list
prediction_list <- list()

for (species_name in names(occupancy_results)) {
  if (occupancy_results[[species_name]]$success) {
    
    model <- occupancy_results[[species_name]]$model
    
    # Get covariates
    new_covs <- site_landcover_vars[, c("forest_pct", "grass_pct", "urban_pct", 
                                        "crops_pct", "shrub_pct", "flooded_veg_pct", 
                                        "habitat_diversity")]
    
    # BETTER SCALING: Handle zero-variance columns
    new_covs_scaled <- new_covs
    for(col in names(new_covs)) {
      if(sd(new_covs[[col]], na.rm = TRUE) > 0) {
        new_covs_scaled[[col]] <- scale(new_covs[[col]])[,1]
      } else {
        new_covs_scaled[[col]] <- 0  # Set constant columns to 0
      }
    }
    
    # Check for NaN values
    if(any(is.na(new_covs_scaled)) || any(is.nan(as.matrix(new_covs_scaled)))) {
      cat("NaN values found in covariates for", species_name, "- skipping\n")
      next
    }
    
    # Get model coefficients
    beta_mean <- apply(model$beta.samples, 2, mean)
    
    # Create design matrix
    design_matrix <- as.matrix(cbind(1, new_covs_scaled))
    
    # Calculate predictions
    linear_pred <- design_matrix %*% beta_mean
    occupancy_pred <- plogis(linear_pred)
    
    prediction_list[[species_name]] <- data.frame(
      site_id = site_landcover_vars$site_id,
      species = species_name,
      occupancy_prob = as.vector(occupancy_pred)
    )
    
    cat("Predictions for", species_name, "- range:", 
        round(range(occupancy_pred), 3), "\n")
  }
}

# Step 6: Create dominant species map

all_predictions <- do.call(rbind, prediction_list)

dominant_species_data <- all_predictions %>%
  group_by(site_id) %>%
  slice_max(occupancy_prob, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  filter(occupancy_prob > 0.1)  # Adjust threshold as needed

# Join with spatial grid
final_map_data <- site_grid %>%
  filter(site_id %in% dominant_species_data$site_id) %>%
  left_join(dominant_species_data, by = "site_id")

# Create the hexagonal-style dominant species map
hexagonal_dominant_map <- ggplot(final_map_data) +
  geom_sf(aes(fill = species), color = "white", size = 0.3) +
  scale_fill_viridis_d(name = "Dominant\nSpecies", option = "plasma", 
                       labels = function(x) substr(x, 1, 15)) +  # Shorten long names
  theme_void() +
  theme(
    legend.position = "right",
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5, 
                              margin = margin(b = 20)),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10, face = "italic"),
    legend.key.size = unit(1, "cm"),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  ) +
  labs(title = "Habitat Partitioning Within Urban Green Space") +
  coord_sf(expand = FALSE)

print(hexagonal_dominant_map)

# Alternative: Show dominant species with occupancy probability as transparency
hexagonal_with_confidence <- ggplot(final_map_data) +
  geom_sf(aes(fill = species, alpha = occupancy_prob), color = "white", size = 0.3) +
  scale_fill_viridis_d(name = "Dominant\nSpecies", option = "plasma") +
  scale_alpha_continuous(name = "Occupancy\nProbability", 
                         range = c(0.4, 1.0),
                         breaks = c(0.2, 0.4, 0.6, 0.8)) +
  theme_void() +
  theme(
    legend.position = "right",
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 9, face = "italic")
  ) +
  labs(title = "Species Partitioning with Occupancy Confidence") +
  coord_sf(expand = FALSE)

print(hexagonal_with_confidence)

# Save both versions
ggsave("hexagonal_dominant_species.png", plot = hexagonal_dominant_map, 
       width = 12, height = 8, dpi = 300, bg = "white")

ggsave("hexagonal_with_confidence.png", plot = hexagonal_with_confidence, 
       width = 12, height = 8, dpi = 300, bg = "white")

################ TESTING

# Check 1: Do you have predictions?
cat("Number of species with predictions:", length(prediction_list), "\n")
if(length(prediction_list) > 0) {
  cat("First few predictions:\n")
  print(head(prediction_list[[1]]))
}

# Check 2: Do you have combined predictions?
if(exists("all_predictions")) {
  cat("All predictions rows:", nrow(all_predictions), "\n")
  cat("Occupancy range:", range(all_predictions$occupancy_prob), "\n")
} else {
  cat("all_predictions doesn't exist - combining predictions now\n")
  all_predictions <- do.call(rbind, prediction_list)
}

# Check 3: How many meet the threshold?
if(exists("all_predictions")) {
  high_occ <- all_predictions %>% filter(occupancy_prob > 0.3)
  cat("Cells with occupancy > 0.3:", nrow(high_occ), "\n")
  
  # Try lower threshold
  medium_occ <- all_predictions %>% filter(occupancy_prob > 0.1)
  cat("Cells with occupancy > 0.1:", nrow(medium_occ), "\n")
}

# Check 4: Lower the threshold and try again
dominant_species_data <- all_predictions %>%
  group_by(site_id) %>%
  slice_max(occupancy_prob, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  filter(occupancy_prob > 0.1)  # Lower threshold

cat("Dominant species data rows:", nrow(dominant_species_data), "\n")

############## Multi Site Occupancy Maps

# Load in shapefile

multi_site <- st_read("/Users/ianbecker/Library/CloudStorage/OneDrive-TheUniversityofTexas-RioGrandeValley/CampusBirds/ugs/urban_hotspot_shapefile")
multi_site <- st_transform(multi_site, crs = 32614)

# Separate polygons and add unique id column

individual_polygons <- st_cast(multi_site, "POLYGON")
individual_polygons$site_name <- paste0("Site_", 1:nrow(individual_polygons))

# Extracting species to create single color legend

all_species <- names(occupancy_results)[sapply(occupancy_results, function(x) x$success)]
cat("Total species to map:", length(all_species), "\n")

# Creating species palette

species_colors <- setNames(
  viridis(length(all_species), option = "turbo"),  
  all_species
)

# Loop through polygons and create maps by site

all_maps <- list()
all_map_data <- list()

for(polygon_idx in 1:nrow(individual_polygons)) {
  
  cat("\n======================================\n")
  cat("Processing polygon", polygon_idx, "of", nrow(individual_polygons), "\n")
  cat("======================================\n")
  
  # Extract single polygon
  single_site <- individual_polygons[polygon_idx, ]
  site_name <- single_site$site_name
  
  # Step 4: Create grid for this polygon
  cell_size <- 75
  site_grid_geom <- st_make_grid(single_site, cellsize = cell_size, square = FALSE)
  
  site_grid <- st_sf(
    site_id = 1:length(site_grid_geom), 
    geometry = site_grid_geom
  )
  
  st_crs(site_grid) <- st_crs(single_site)
  
  grid_intersects <- st_intersects(site_grid, single_site, sparse = FALSE)[,1]
  site_grid <- site_grid[grid_intersects, ]
  site_grid$site_id <- 1:nrow(site_grid)
  
  cat("Created", nrow(site_grid), "grid cells for", site_name, "\n")
  
  # Skip if no grid cells (polygon too small)
  if(nrow(site_grid) == 0) {
    cat("No grid cells for", site_name, "- skipping\n")
    next
  }
  
  # Step 5: Extract landcover data
  site_landcover_vars <- data.frame(
    site_id = integer(),
    forest_pct = numeric(), 
    grass_pct = numeric(), 
    crops_pct = numeric(),
    urban_pct = numeric(),
    shrub_pct = numeric(),
    flooded_veg_pct = numeric(),
    habitat_diversity = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (i in 1:nrow(site_grid)) {
    
    if (i %% 10 == 0) cat("  Extracting landcover for cell", i, "of", nrow(site_grid), "\n")
    
    grid_cell <- site_grid[i,]
    lc_values <- terra::extract(landcover, vect(grid_cell), df = TRUE)
    
    if (nrow(lc_values) == 0) {
      cell_results <- data.frame(
        site_id = grid_cell$site_id,
        forest_pct = 0, grass_pct = 0, crops_pct = 0, urban_pct = 0, 
        shrub_pct = 0, flooded_veg_pct = 0, habitat_diversity = 0
      )
    } else {
      lc_column <- names(lc_values)[2]
      total_pixels <- nrow(lc_values)
      
      lc_summary <- lc_values %>%
        count(.data[[lc_column]]) %>%
        mutate(percentage = n / total_pixels * 100)
      
      names(lc_summary)[1] <- "class"
      
      forest_pct <- ifelse(1 %in% lc_summary$class, 
                           lc_summary$percentage[lc_summary$class == 1], 0)
      grass_pct <- ifelse(2 %in% lc_summary$class, 
                          lc_summary$percentage[lc_summary$class == 2], 0)
      crops_pct <- ifelse(4 %in% lc_summary$class, 
                          lc_summary$percentage[lc_summary$class == 4], 0)
      shrub_pct <- ifelse(5 %in% lc_summary$class,
                          lc_summary$percentage[lc_summary$class == 5], 0)
      urban_pct <- ifelse(6 %in% lc_summary$class, 
                          lc_summary$percentage[lc_summary$class == 6], 0)
      flooded_veg_pct <- ifelse(3 %in% lc_summary$class,
                                lc_summary$percentage[lc_summary$class == 3], 0)
      
      proportions <- lc_summary$percentage / 100
      proportions <- proportions[proportions > 0]
      habitat_diversity <- -sum(proportions * log(proportions), na.rm = TRUE)
      
      cell_results <- data.frame(
        site_id = grid_cell$site_id,
        forest_pct = forest_pct,
        grass_pct = grass_pct, 
        crops_pct = crops_pct,
        urban_pct = urban_pct,
        shrub_pct = shrub_pct,
        flooded_veg_pct = flooded_veg_pct,
        habitat_diversity = habitat_diversity
      )
    }
    
    site_landcover_vars <- rbind(site_landcover_vars, cell_results)
  }
  
  # Step 6: Generate predictions
  prediction_list <- list()
  
  for (species_name in names(occupancy_results)) {
    if (occupancy_results[[species_name]]$success) {
      
      model <- occupancy_results[[species_name]]$model
      
      new_covs <- site_landcover_vars[, c("forest_pct", "grass_pct", "urban_pct", 
                                          "crops_pct", "shrub_pct", "flooded_veg_pct", 
                                          "habitat_diversity")]
      
      new_covs_scaled <- new_covs
      for(col in names(new_covs)) {
        if(sd(new_covs[[col]], na.rm = TRUE) > 0) {
          new_covs_scaled[[col]] <- scale(new_covs[[col]])[,1]
        } else {
          new_covs_scaled[[col]] <- 0
        }
      }
      
      if(any(is.na(new_covs_scaled)) || any(is.nan(as.matrix(new_covs_scaled)))) {
        cat("NaN values found in covariates for", species_name, "- skipping\n")
        next
      }
      
      beta_mean <- apply(model$beta.samples, 2, mean)
      design_matrix <- as.matrix(cbind(1, new_covs_scaled))
      linear_pred <- design_matrix %*% beta_mean
      occupancy_pred <- plogis(linear_pred)
      
      prediction_list[[species_name]] <- data.frame(
        site_id = site_landcover_vars$site_id,
        species = species_name,
        occupancy_prob = as.vector(occupancy_pred)
      )
    }
  }
  
  if(length(prediction_list) == 0) {
    cat("No predictions generated for", site_name, "- skipping\n")
    next
  }
  
  # Step 7: Create dominant species map for this polygon
  all_predictions <- do.call(rbind, prediction_list)
  
  dominant_species_data <- all_predictions %>%
    group_by(site_id) %>%
    slice_max(occupancy_prob, n = 1, with_ties = FALSE) %>%
    ungroup() 
  
  if(nrow(dominant_species_data) == 0) {
    cat("No dominant species above threshold for", site_name, "- skipping\n")
    next
  }
  
  final_map_data <- site_grid %>%
    filter(site_id %in% dominant_species_data$site_id) %>%
    left_join(dominant_species_data, by = "site_id")
  
  # Create map
  hexagonal_dominant_map <- ggplot(final_map_data) +
    geom_sf(aes(fill = species), color = "white", size = 0.3) +
    scale_fill_manual(                                          # <-- USE THIS INSTEAD
      name = "Dominant\nSpecies",
      values = species_colors,                                  # <-- Uses your palette
      labels = function(x) substr(x, 1, 15),
      drop = FALSE
    ) +
    theme_void() +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5, 
                                margin = margin(b = 20)),
      legend.title = element_text(size = 12, face = "bold"),
      legend.text = element_text(size = 10, face = "italic"),
      legend.key.size = unit(1, "cm"),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    ) +
    labs(title = paste("Occupancy Map -", site_name)) +
    coord_sf(expand = FALSE)
  
  print(hexagonal_dominant_map)
  
  # Save map
  ggsave(paste0("hexagonal_dominant_species_", site_name, ".png"), 
         plot = hexagonal_dominant_map, 
         width = 12, height = 8, dpi = 300, bg = "white")
  
  # Store for later use
  all_maps[[site_name]] <- hexagonal_dominant_map
  all_map_data[[site_name]] <- final_map_data
  
  cat("Completed map for", site_name, "\n")
}

cat("\n======================================\n")
cat("Processing complete! Generated", length(all_maps), "maps\n")
cat("======================================\n")

# Optional: Create a multi-panel figure with all sites
if(length(all_maps) > 1 && length(all_maps) <= 6) {
  library(patchwork)
  
  combined_plot <- wrap_plots(all_maps, ncol = 2)
  
  ggsave("all_sites_combined.png", plot = combined_plot, 
         width = 20, height = 10 * ceiling(length(all_maps)/2), dpi = 300, bg = "white")
  
  cat("Combined multi-panel figure saved as 'all_sites_combined.png'\n")
}

# 

dummy_data <- data.frame(
  x = rep(1, length(all_species)),
  y = 1:length(all_species),
  species = factor(all_species, levels = all_species)
)

# Create a plot just to generate the legend
legend_plot <- ggplot(dummy_data, aes(x = x, y = y, fill = species)) +
  geom_tile() +  # or geom_point(shape = 22, size = 10)
  scale_fill_manual(
    name = "Dominant Species",
    values = species_colors,
    labels = function(x) substr(x, 1, 25)  # Adjust character limit as needed
  ) +
  theme_void() +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 16, face = "bold"),
    legend.text = element_text(size = 12, face = "italic"),
    legend.key.size = unit(1.2, "cm"),
    legend.key = element_rect(color = "white", linewidth = 0.5)
  )

# Extract just the legend
library(ggpubr)
standalone_legend <- get_legend(legend_plot)

# Convert to ggplot object and save
legend_as_plot <- as_ggplot(standalone_legend)

ggsave("species_legend_only.png", 
       plot = legend_as_plot, 
       width = 4,    # Adjust width as needed
       height = 8,   # Adjust height based on number of species
       dpi = 300, 
       bg = "white")

cat("Standalone legend saved as 'species_legend_only.png'\n")
