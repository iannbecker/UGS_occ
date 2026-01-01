##############################
#
# iNat Data Pull 
# 12/19/2025
# Ian Becker
#
##############################

# this script is used for pulling and filtering inat data for select species
# in our urban green spaces (OLD)

library(rinat)
library(dplyr)
library(sf)

####################
#   iNat Data Pull
####################

################## Data Prep

# Reading in shapefile

urban_hotspot <- st_read("lrgv_osm_green_spaces_filtered")

cat("Original CRS:", st_crs(urban_hotspot)$input, "\n")

# Transform to WGS84 for iNaturalist API (expects lat/lon)
urban_hotspot_wgs84 <- st_transform(urban_hotspot, crs = 4326)

cat("Transformed to WGS84 for API calls\n")

# Shapefile conversion to bounding box

urban_hotspot_bbox <- st_bbox(urban_hotspot_wgs84)
bounds <- c(urban_hotspot_bbox$ymin, urban_hotspot_bbox$xmin, urban_hotspot_bbox$ymax, urban_hotspot_bbox$xmax)

cat("Bounding box (lat/lon):\n")
cat("  South:", bounds[1], "\n")
cat("  West:", bounds[2], "\n")
cat("  North:", bounds[3], "\n")
cat("  East:", bounds[4], "\n\n")

# Species list from CSV

species_df <- read.csv("exploratory_species.csv")
species_list <- species_df$scientific_name

cat("Loading", length(species_list), "species:\n")
print(species_list)
cat("\n")

# Prepping for data pull

all_rinat_data <- list()
inat_max_records <- 10000

################## Data Pull

for (i in seq_along(species_list)) {
  
  species_name <- species_list[i]
  cat(paste0("[", i, "/", length(species_list), "] Extracting data for: ", species_name, "\n"))
  
  # Try to get the data with error handling
  tryCatch({
    
    # Get observations from iNaturalist directly
    species_data <- get_inat_obs(
      taxon_name = species_name,
      bounds = bounds,
      maxresults = inat_max_records,
      quality = "research"  # research grade only
    )
    
    if (!is.null(species_data) && nrow(species_data) > 0) {
      
      # Filter to actual polygon (not just bounding box)
      # Convert to spatial points
      species_sf <- st_as_sf(species_data, 
                             coords = c("longitude", "latitude"), 
                             crs = 4326)
      
      # Check which points intersect with ANY of the green space cells
      # Use st_intersects and check if any cell contains the point
      intersections <- st_intersects(species_sf, urban_hotspot_wgs84)
      inside_polygon <- lengths(intersections) > 0
      
      # Filter to only points inside polygon
      species_filtered <- species_data[inside_polygon, ]
      
      if (nrow(species_filtered) > 0) {
        species_filtered$query_species <- species_name
        all_rinat_data[[species_name]] <- species_filtered
        
        cat("  Found", nrow(species_data), "in bounding box,", 
            nrow(species_filtered), "inside polygon\n")
      } else {
        cat("  Found", nrow(species_data), "in bounding box, but 0 inside polygon\n")
      }
      
    } else {
      cat("  No records found\n")
    }
    
  }, error = function(e) {
    cat("  Error occurred:", e$message, "\n")
  })
  
  # Be nice to iNaturalist API
  Sys.sleep(1)
  cat("\n")
}


####################
#   Data Filtering
####################

filtered_inat_data <- list()
filtering_summary <- data.frame(
  species = character(),
  raw_records = numeric(),
  after_captive_filter = numeric(),
  after_obscured_filter = numeric(),
  after_accuracy_filter = numeric(),
  final_records = numeric(),
  meets_25_threshold = logical(),
  stringsAsFactors = FALSE
)

for (species_name in names(all_rinat_data)) {
  
  cat("Filtering", species_name, ":\n")
  species_data <- all_rinat_data[[species_name]]
  
  # Track numbers through each filter
  raw_count <- nrow(species_data)
  cat("  Raw records:", raw_count, "\n")
  
  # (1) Research grade entries - already filtered during extraction
  
  # (2) Non-captive and non-cultivated individuals
  species_data <- species_data %>%
    filter(
      captive_cultivated != "true" | is.na(captive_cultivated)
    )
  
  captive_count <- nrow(species_data)
  cat("  After captive/cultivated filter:", captive_count, "\n")
  
  # (3) Spatially unobscured records
  species_data <- species_data %>%
    filter(
      coordinates_obscured != "true" | is.na(coordinates_obscured)
    )
  
  obscured_count <- nrow(species_data)
  cat("  After obscured coordinates filter:", obscured_count, "\n")
  
  # (4) Maximum inaccuracy of 100 meters
  species_data <- species_data %>%
    filter(
      positional_accuracy <= 100 | is.na(positional_accuracy)
    )
  
  accuracy_count <- nrow(species_data)
  cat("  After 100m accuracy filter:", accuracy_count, "\n")
  
  # (5) Species with 25 or more observations (check but don't filter yet)
  meets_threshold <- accuracy_count >= 25
  
  if (meets_threshold) {
    filtered_inat_data[[species_name]] <- species_data
    cat("  ✓ Meets 25+ observation threshold - INCLUDED\n")
  } else {
    cat("  ✗ Below 25 observation threshold - EXCLUDED\n")
  }
  
  # Record summary
  filtering_summary <- rbind(filtering_summary, data.frame(
    species = species_name,
    raw_records = raw_count,
    after_captive_filter = captive_count,
    after_obscured_filter = obscured_count,
    after_accuracy_filter = accuracy_count,
    final_records = accuracy_count,
    meets_25_threshold = meets_threshold
  ))
  
  cat("\n")
}

####################
#   Save Results
####################

saveRDS(filtered_inat_data, "lrgv_flycatcher_inat_data.rds")
cat("Saved filtered data to: lrgv_flycatcher_inat_data.rds\n")

# Save filtering summary
write.csv(filtering_summary, "lrgv_flycatcher_filtering_summary.csv", row.names = FALSE)
cat("Saved filtering summary to: lrgv_flycatcher_filtering_summary.csv\n")

cat("\n=== DATA PULL COMPLETE ===\n")
cat("Species with sufficient data (≥25 obs):", sum(filtering_summary$meets_25_threshold), "\n")
cat("Total filtered observations:", sum(filtering_summary$final_records), "\n")

