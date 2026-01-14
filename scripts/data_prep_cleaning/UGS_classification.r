##############################
#
# Spatial Data prep
# 12/15/2025
# Ian Becker
#
##############################

# In this script I pull urban areas for the LRGV
# I then utilize land cover data to classify urban green spaces

library(tigris)
library(dplyr)
library(sf)
library(terra)

options(tigris_use_cache = TRUE)

##############################
#   LRGV Urban Areas
##############################

# LRGV counties in Texas

lrgv_counties <- c("Starr", "Hidalgo", "Cameron", "Willacy")

# Get Texas counties
tx_counties <- counties(state = "TX", cb = TRUE, year = 2020)

# Filter to LRGV

lrgv <- tx_counties %>%
  filter(NAME %in% lrgv_counties)

# Check what we got

print(lrgv$NAME)

# Save county boundaries

st_write(lrgv, "lrgv_counties.shp", delete_dsn = TRUE)

###################### Download urban areas

cat("\nLoading your existing urban areas shapefile...\n")

# Urban areas shapefile

urban_areas_path <- "tl_2020_us_uac20"

# Load urban areas

urban_areas_all <- st_read(urban_areas_path)

cat("Loaded", nrow(urban_areas_all), "urban areas from shapefile\n")

# Filter to those intersecting LRGV counties

cat("Filtering to LRGV...\n")
lrgv_urban <- urban_areas_all[lrgv, ]

# Transform to UTM

lrgv_urban_utm <- st_transform(lrgv_urban, crs = 32614)

cat("Found", nrow(lrgv_urban), "urban areas in LRGV\n")

# Check urban areas

print(lrgv_urban$NAME20)

# Save urban areas

st_write(lrgv_urban, "lrgv_urban_areas.shp", delete_dsn = TRUE)

##############################
#   Urban Green Spaces
##############################

cat("\nStarting Urban Green Space Classification...\n")

# Load Dynamic World Cover landcover data

landcover <- rast("lrgv_dynamic_world_cover.tif")

# Project landcover to match urban areas CRS (UTM Zone 14N)

landcover_utm <- project(landcover, "EPSG:32614", method = "near")
lrgv_urban_utm <- st_transform(lrgv_urban_areas, crs = 32614)

cat("Landcover reprojected\n")

# Create 150m hexagonal grid over urban areas

cat("\nCreating 150m hexagonal grid cells...\n")

cell_size <- 150  # 150 meters
grid_hexagons <- st_make_grid(lrgv_urban_utm, 
                              cellsize = cell_size, 
                              square = FALSE)  # FALSE = hexagons

cat("Created", length(grid_hexagons), "total hexagonal cells\n")

# Filter to only cells that intersect urban areas

cat("Filtering to cells within urban boundaries...\n")
grid_intersects <- st_intersects(grid_hexagons, lrgv_urban_utm, sparse = FALSE)
grid_urban <- grid_hexagons[apply(grid_intersects, 1, any)]

cat("Retained", length(grid_urban), "cells within urban areas\n")

# Convert to sf object with IDs

grid_urban_sf <- st_sf(
  cell_id = 1:length(grid_urban),
  geometry = grid_urban
)

# Calculate vegetation percentage for each grid cell

cat("\nCalculating vegetation percentage per cell...\n")

# Initialize results dataframe

vegetation_results <- data.frame(
  cell_id = integer(),
  trees_pct = numeric(),
  grass_pct = numeric(),
  shrub_pct = numeric(),
  flooded_veg_pct = numeric(),
  total_veg_pct = numeric(),
  urban_green_space = logical(),
  stringsAsFactors = FALSE
)

# Extract landcover for each cell

for (i in 1:nrow(grid_urban_sf)) {
  
  if (i %% 100 == 0) {
    cat("Processing cell", i, "of", nrow(grid_urban_sf), 
        "(", round(i/nrow(grid_urban_sf)*100, 1), "%)\n")
  }
  
  grid_cell <- grid_urban_sf[i, ]
  
  # Extract landcover values within this cell
  lc_values <- terra::extract(landcover_utm, vect(grid_cell), df = TRUE)
  
  if (nrow(lc_values) == 0) {
    # No data - assume not green space
    cell_result <- data.frame(
      cell_id = grid_cell$cell_id,
      trees_pct = 0,
      grass_pct = 0,
      shrub_pct = 0,
      flooded_veg_pct = 0,
      total_veg_pct = 0,
      urban_green_space = FALSE
    )
  } else {
    
    # Get the landcover column (should be the second column after ID)
    lc_column <- names(lc_values)[2]
    total_pixels <- nrow(lc_values)
    
    # Count pixels by class
    lc_summary <- lc_values %>%
      count(.data[[lc_column]]) %>%
      mutate(percentage = n / total_pixels * 100)
    
    names(lc_summary)[1] <- "class"
    
    # Dynamic World Cover classes:
    # 1 = Trees/Forest
    # 2 = Grass
    # 5 = Shrub/Scrub  
    # 3 = Flooded vegetation
    
    trees_pct <- ifelse(1 %in% lc_summary$class, 
                        lc_summary$percentage[lc_summary$class == 1], 0)
    grass_pct <- ifelse(2 %in% lc_summary$class, 
                        lc_summary$percentage[lc_summary$class == 2], 0)
    shrub_pct <- ifelse(5 %in% lc_summary$class,
                        lc_summary$percentage[lc_summary$class == 5], 0)
    flooded_veg_pct <- ifelse(3 %in% lc_summary$class,
                              lc_summary$percentage[lc_summary$class == 3], 0)
    
    # Calculate total vegetation percentage
    total_veg_pct <- trees_pct + grass_pct + shrub_pct + flooded_veg_pct
    
    # Classify as urban green space if ≥50% vegetation
    is_green_space <- total_veg_pct >= 50
    
    cell_result <- data.frame(
      cell_id = grid_cell$cell_id,
      trees_pct = trees_pct,
      grass_pct = grass_pct,
      shrub_pct = shrub_pct,
      flooded_veg_pct = flooded_veg_pct,
      total_veg_pct = total_veg_pct,
      urban_green_space = is_green_space
    )
  }
  
  vegetation_results <- rbind(vegetation_results, cell_result)
}

cat("\nVegetation analysis complete!\n")

# Join vegetation results back to spatial grid
grid_classified <- grid_urban_sf %>%
  left_join(vegetation_results, by = "cell_id")

# Summary statistics
cat("\n=== URBAN GREEN SPACE SUMMARY ===\n")
cat("Total grid cells analyzed:", nrow(grid_classified), "\n")
cat("Cells classified as urban green space (≥50% veg):", 
    sum(grid_classified$urban_green_space), "\n")
cat("Percentage of cells that are green space:", 
    round(sum(grid_classified$urban_green_space) / nrow(grid_classified) * 100, 1), "%\n")

# Vegetation composition summary
cat("\n=== VEGETATION COMPOSITION (all cells) ===\n")
cat("Mean % Trees:", round(mean(grid_classified$trees_pct), 1), "\n")
cat("Mean % Grass:", round(mean(grid_classified$grass_pct), 1), "\n")
cat("Mean % Shrub/Scrub:", round(mean(grid_classified$shrub_pct), 1), "\n")
cat("Mean % Flooded Veg:", round(mean(grid_classified$flooded_veg_pct), 1), "\n")
cat("Mean % Total Vegetation:", round(mean(grid_classified$total_veg_pct), 1), "\n")

# Green space cells only
green_cells <- grid_classified %>% filter(urban_green_space == TRUE)

if (nrow(green_cells) > 0) {
  cat("\n=== VEGETATION COMPOSITION (green space cells only) ===\n")
  cat("Mean % Trees:", round(mean(green_cells$trees_pct), 1), "\n")
  cat("Mean % Grass:", round(mean(green_cells$grass_pct), 1), "\n")
  cat("Mean % Shrub/Scrub:", round(mean(green_cells$shrub_pct), 1), "\n")
  cat("Mean % Flooded Veg:", round(mean(green_cells$flooded_veg_pct), 1), "\n")
  cat("Mean % Total Vegetation:", round(mean(green_cells$total_veg_pct), 1), "\n")
}

# Save results
cat("\nSaving results...\n")

# Save all classified cells
st_write(grid_classified, "lrgv_urban_grid_150m_classified.shp", delete_dsn = TRUE)
cat("Saved: lrgv_urban_grid_150m_classified.shp\n")

# Save only urban green space cells
st_write(green_cells, "lrgv_urban_green_spaces_150m.shp", delete_dsn = TRUE)
cat("Saved: lrgv_urban_green_spaces_150m.shp\n")

# Save vegetation results as CSV
write.csv(vegetation_results, "lrgv_grid_vegetation_data.csv", row.names = FALSE)
cat("Saved: lrgv_grid_vegetation_data.csv\n")

cat("\n=== PROCESSING COMPLETE ===\n")
cat("You now have:\n")
cat("1. All grid cells with vegetation data (lrgv_urban_grid_150m_classified.shp)\n")
cat("2. Only urban green space cells (lrgv_urban_green_spaces_150m.shp)\n")
cat("3. Vegetation data table (lrgv_grid_vegetation_data.csv)\n")






