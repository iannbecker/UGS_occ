##############################
#
# iNat Data Pull 
# 8/6/2025
# Ian Becker
#
##############################

library(rgbif)
library(rinat)
library(dplyr)
library(sf)

####################
#   GBIF
####################

################## Data Prep

# Reading in shapefile

urban_hotspot <- st_read("C:/Users/ianbe/OneDrive - The University of Texas-Rio Grande Valley/CampusBirds/ugs/urban_hotspot_shapefile")

# Converting to wkt format for rgbif

urban_hotspot_wkt <- st_as_text(st_geometry(urban_hotspot)[1])

# Species list

species <- read.csv("patch_occupancy_candidates.csv")
species_list <- unique(species$scientific_name)

# Prepping for data pull 

inat_dataset_key <- "50c9509d-22c7-4a22-a47d-8c48425ef4a7"

all_gbif_data <- list()
gbif_max_records <- 5000 

################## Data Pull

for (i in seq_along(species_list)) {
  
  species <- species_list[i]
  cat(paste0("[", i, "/", length(species_list), "] Extracting data for: ", species, "\n"))
  
  # Try to get the data with error handling
  tryCatch({
    
    # Search for the species
    species_data <- occ_search(
      scientificName = species,
      datasetKey = inat_dataset_key,
      geometry = urban_hotspot_wkt,
      hasCoordinate = TRUE,
      hasGeospatialIssue = FALSE,
      limit = max_records
    )
    
    # Extract and process the data frame
    if (!is.null(species_data$data) && nrow(species_data$data) > 0) {
      
      df <- species_data$data
      df$query_species <- species  # Add species name for tracking
      
      # Store the data
      all_gbif_data[[species]] <- df
      
      cat("  Found", nrow(df), "records\n")
      
    } else {
      cat("  No records found\n")
    }
    
  }, error = function(e) {
    cat("  Error occurred:", e$message, "\n")
  })
  
  # Be nice to the API - add a small delay
  Sys.sleep(1)
  cat("\n")
}


####################
#   iNat
####################

################## Data Prep

# Shapefile conversion to bounding box

urban_hotspot_bbox <- st_bbox(urban_hotspot)
bounds <- c(urban_hotspot_bbox$ymin, urban_hotspot_bbox$xmin, urban_hotspot_bbox$ymax, urban_hotspot_bbox$xmax)

# Prepping for data pull

all_rinat_data <- list()
inat_max_records <- 10000

################## Data Pull

for (i in seq_along(species_list)) {
  
  species <- species_list[i]
  cat(paste0("[", i, "/", length(species_list), "] Extracting data for: ", species, "\n"))
  
  # Try to get the data with error handling
  tryCatch({
    
    # Get observations from iNaturalist directly
    species_data <- get_inat_obs(
      taxon_name = species,
      bounds = bounds,
      maxresults = inat_max_records,
      quality = "research"  # research grade only, or use "any" for all
    )
    
    if (!is.null(species_data) && nrow(species_data) > 0) {
      
      # Filter to actual polygon (not just bounding box)
      # Convert to spatial points
      species_sf <- st_as_sf(species_data, 
                             coords = c("longitude", "latitude"), 
                             crs = 4326)
      
      # Check which points are inside the polygon
      inside_polygon <- st_within(species_sf, urban_hotspot, sparse = FALSE)[,1]
      
      # Filter to only points inside polygon
      species_filtered <- species_data[inside_polygon, ]
      
      if (nrow(species_filtered) > 0) {
        species_filtered$query_species <- species
        all_rinat_data[[species]] <- species_filtered
        
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
#   Data Filtering (inat)
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
  # (assuming you used quality = "research" in get_inat_obs)
  
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

########## Saving 

saveRDS(filtered_inat_data, "urban_hotspot_species_inat.rds")
