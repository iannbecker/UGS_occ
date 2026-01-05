##############################
#
# iNat Data Processing and Detection Matrix Creation
# Ian Becker
#
##############################

# This script:
# 1. Loads raw iNat data
# 2. Applies quality filters
# 3. Assigns observations to sites
# 4. Creates detection matrices for occupancy modeling

library(sf)
library(dplyr)
library(lubridate)

####################
#   Setup Output Directory
####################

# Create output directory for detection matrices
if (!dir.exists("detection_matrices")) {
  dir.create("detection_matrices")
  cat("Created detection_matrices/ directory\n\n")
} else {
  cat("Using existing detection_matrices/ directory\n\n")
}

####################
#   Part 1: Load and Filter iNat Data
####################

cat("=== PART 1: LOADING AND FILTERING iNAT DATA ===\n\n")

# Load raw iNat download
cat("Loading raw iNaturalist data...\n")
inat_raw <- read.csv(file.path("observations-660979.csv/observations-660979.csv"))  # CHANGE THIS to your actual file

cat("Loaded", nrow(inat_raw), "raw observations\n\n")

# Check column names (adjust as needed)
cat("Available columns:\n")
print(names(inat_raw))
cat("\n")

####################
#   Apply Quality Filters
####################

cat("Applying quality filters...\n")

inat_filtered <- inat_raw %>%
  filter(
    # Research grade only
    quality_grade == "research",
    
    # Non-captive/cultivated
    captive_cultivated == "false" | is.na(captive_cultivated),
    
    # Unobscured coordinates  
    coordinates_obscured == "false" | is.na(coordinates_obscured),
    
    # Good spatial accuracy (≤100m)
    positional_accuracy <= 100 | is.na(positional_accuracy),
    
    # Has valid coordinates
    !is.na(longitude) & !is.na(latitude) & 
      longitude != "" & latitude != ""
  )

cat("After quality filtering:", nrow(inat_filtered), "observations\n")
cat("Removed:", nrow(inat_raw) - nrow(inat_filtered), "observations\n\n")

####################
#   Part 2: Load Sites and Assign Observations
####################

cat("=== PART 2: ASSIGNING OBSERVATIONS TO SITES ===\n\n")

# Load urban green space sites
sites <- st_read("lrgv_green_spaces_detection_filtered")
n_sites <- nrow(sites)

cat("Loaded", n_sites, "urban green space sites\n\n")

# Convert iNat data to spatial
cat("Converting iNat data to spatial format...\n")

# Ensure coordinates are numeric
inat_filtered$longitude <- as.numeric(inat_filtered$longitude)
inat_filtered$latitude <- as.numeric(inat_filtered$latitude)

# Remove any remaining invalid coordinates
inat_filtered <- inat_filtered %>%
  filter(!is.na(longitude) & !is.na(latitude))

# Create sf object
inat_sf <- st_as_sf(inat_filtered,
                    coords = c("longitude", "latitude"),
                    crs = 4326)

# Transform to match sites CRS
if (st_crs(inat_sf) != st_crs(sites)) {
  inat_sf <- st_transform(inat_sf, st_crs(sites))
}

cat("Converted to spatial format\n\n")

####################
#   Spatial Assignment to Sites
####################

cat("Assigning observations to sites (this may take a minute)...\n")

# Find which site each observation falls in
obs_in_sites <- st_intersects(inat_sf, sites)

# Add site_id to observations
inat_sf$site_id <- sapply(obs_in_sites, function(x) {
  if (length(x) > 0) return(x[1])
  else return(NA)
})

# Filter to observations within sites
inat_in_sites <- inat_sf %>%
  filter(!is.na(site_id))

cat("\n")
cat("Total observations within sites:", nrow(inat_in_sites), "\n")
cat("Observations outside sites:", nrow(inat_sf) - nrow(inat_in_sites), "\n")

# How many sites have data?
sites_with_data <- n_distinct(inat_in_sites$site_id)
cat("Sites with at least 1 observation:", sites_with_data, "of", n_sites,
    "(", round(sites_with_data/n_sites*100, 1), "%)\n\n")

####################
#   Save Filtered and Assigned Data
####################

cat("Saving filtered iNat data with site assignments...\n")

# Convert back to regular dataframe for easier use
inat_with_sites <- inat_in_sites %>%
  st_drop_geometry()

write.csv(inat_with_sites, "inat_observations_with_sites.csv", row.names = FALSE)

cat("Saved: inat_observations_with_sites.csv\n\n")

####################
#   Species Summary
####################

cat("=== SPECIES SUMMARY ===\n\n")

species_summary <- inat_with_sites %>%
  group_by(common_name) %>%
  summarise(
    n_observations = n(),
    n_sites = n_distinct(site_id),
    .groups = "drop"
  ) %>%
  filter(n_observations >= 25) %>%  # Occupancy modeling threshold
  arrange(desc(n_observations))

cat("Species with ≥25 observations:\n")
cat("Total:", nrow(species_summary), "species\n\n")

# Show top 20
print(head(species_summary, 20))

# Save species summary
write.csv(species_summary, "species_summary_in_sites.csv", row.names = FALSE)
cat("\nSaved: species_summary_in_sites.csv\n\n")

####################
#   Part 3: Build Detection Matrices for All Species
####################

cat("=== PART 3: BUILDING DETECTION MATRICES ===\n\n")

# Define all species for analysis (15 total)

# Urban Waterbirds (6 species)
waterbirds <- c(
  "Black-bellied Whistling-Duck",
  "Green Heron",
  "Neotropic Cormorant",
  "Great Blue Heron",
  "Great Egret"
)

# Urban Generalists (5 species)
generalists <- c(
  "Great-tailed Grackle",
  "Northern Mockingbird",
  "White-winged Dove",
  "Inca Dove",
  "Mourning Dove"
)

# Urban Specialists (4 species)
specialists <- c(
  "Great Kiskadee",
  "Tropical Kingbird",
  "Green Jay",
  "Long-billed Thrasher",
  "Golden-fronted Woodpecker"
)

# Combine all species
test_species_list <- c(waterbirds, generalists, specialists)

cat("Building detection matrices for", length(test_species_list), "species:\n\n")
cat("URBAN WATERBIRDS (n=5):\n")
for (sp in waterbirds) cat("  -", sp, "\n")
cat("\nURBAN GENERALISTS (n=5):\n")
for (sp in generalists) cat("  -", sp, "\n")
cat("\nURBAN SPECIALISTS (n=5):\n")
for (sp in specialists) cat("  -", sp, "\n")
cat("\n")

####################
#   Prepare Temporal Data
####################

# Add date information
inat_with_sites$date <- as.Date(inat_with_sites$observed_on)
inat_with_sites$year <- year(inat_with_sites$date)
inat_with_sites$month <- month(inat_with_sites$date)

# Filter to study period
inat_with_sites <- inat_with_sites %>%
  filter(year >= 2015, year <= 2025)

cat("Filtered to 2015-2025 study period\n")
cat("Total observations:", nrow(inat_with_sites), "\n\n")

####################
#   Define Temporal Structure
####################

# Primary periods (years)
years <- 2015:2025
n_years <- length(years)

# Secondary occasions (months)
months <- 1:12
n_months <- length(months)

cat("Temporal structure:\n")
cat("  Primary periods (years):", n_years, "\n")
cat("  Secondary occasions (months):", n_months, "\n\n")

####################
#   Loop Through Species
####################

# Store all matrices and metadata
all_detection_matrices <- list()
all_metadata <- list()

for (test_species in test_species_list) {
  
  # Wrap everything in error handling
  tryCatch({
    
    cat("\n========================================\n")
    cat("PROCESSING:", test_species, "\n")
    cat("========================================\n\n")
    
    ####################
    #   Filter to Species
    ####################
    
    species_data <- inat_with_sites %>%
      filter(common_name == test_species)
    
    cat("Total observations:", nrow(species_data), "\n")
    cat("Sites with observations:", n_distinct(species_data$site_id), "\n\n")
    
    if (nrow(species_data) == 0) {
      cat("No observations found - skipping\n")
      next
    }
    
    ####################
    #   Create Detection Matrix
    ####################
    
    cat("Creating detection matrix [", n_sites, "sites ×", n_years, "years ×", 
        n_months, "months]...\n")
    
    # Initialize 3D array
    detection_matrix <- array(0, 
                              dim = c(n_sites, n_years, n_months),
                              dimnames = list(
                                site = 1:n_sites,
                                year = years,
                                month = months
                              ))
    
    # Fill in detections
    for (i in 1:nrow(species_data)) {
      
      if (i %% 100 == 0) cat("  Processing observation", i, "of", nrow(species_data), "\n")
      
      obs <- species_data[i, ]
      
      # Get indices
      site_idx <- obs$site_id
      year_idx <- which(years == obs$year)
      month_idx <- obs$month
      
      # Skip if invalid
      if (is.na(site_idx) || site_idx < 1 || site_idx > n_sites) next
      if (length(year_idx) == 0) next
      if (is.na(month_idx) || month_idx < 1 || month_idx > 12) next
      
      # Mark as detected
      detection_matrix[site_idx, year_idx, month_idx] <- 1
    }
    
    cat("\nDetection matrix complete!\n\n")
    
    ####################
    #   Matrix Summary
    ####################
    
    cat("=== DETECTION MATRIX SUMMARY ===\n\n")
    
    # Sites with detections
    sites_with_detections <- apply(detection_matrix, 1, function(x) any(x == 1))
    n_sites_detected <- sum(sites_with_detections)
    
    cat("Sites with detections:", n_sites_detected, "of", n_sites, 
        "(", round(n_sites_detected/n_sites*100, 1), "%)\n\n")
    
    # Detections by year
    detections_by_year <- apply(detection_matrix, 2, sum)
    cat("Detections by year (site-months):\n")
    for (i in 1:length(years)) {
      cat("  ", years[i], ":", detections_by_year[i], "\n")
    }
    
    # Naive occupancy
    naive_occupancy <- n_sites_detected / n_sites
    cat("\nNaive occupancy:", round(naive_occupancy, 3), "\n\n")
    
    ####################
    #   Example Detection Histories
    ####################
    
    cat("=== EXAMPLE DETECTION HISTORIES ===\n\n")
    
    # Show first 2 sites with detections
    sites_to_show <- which(sites_with_detections)[1:min(2, n_sites_detected)]
    
    for (site_idx in sites_to_show) {
      
      cat("Site", site_idx, ":\n")
      site_history <- detection_matrix[site_idx, , ]
      
      for (y in 1:min(3, n_years)) {  # Just show first 3 years
        detections <- paste(site_history[y, ], collapse = "")
        n_detections <- sum(site_history[y, ])
        cat("  ", years[y], ":", detections, "(", n_detections, "months)\n")
      }
      cat("  ...\n")
    }
    cat("\n")
    
    ####################
    #   Save Results
    ####################
    
    cat("Saving results...\n")
    
    # Save detection matrix
    species_filename <- gsub(" ", "_", test_species)
    saveRDS(detection_matrix, 
            paste0("detection_matrices/detection_matrix_", species_filename, ".rds"))
    cat("Saved: detection_matrices/detection_matrix_", species_filename, ".rds\n")
    
    # Save metadata
    metadata <- list(
      species = test_species,
      n_sites = n_sites,
      n_years = n_years,
      n_months = n_months,
      years = years,
      months = months,
      n_observations = nrow(species_data),
      n_sites_detected = n_sites_detected,
      naive_occupancy = naive_occupancy,
      date_created = Sys.Date()
    )
    
    saveRDS(metadata, 
            paste0("detection_matrices/detection_metadata_", species_filename, ".rds"))
    cat("Saved: detection_matrices/detection_metadata_", species_filename, ".rds\n")
    
    # Store in lists
    all_detection_matrices[[test_species]] <- detection_matrix
    all_metadata[[test_species]] <- metadata
    
    cat("\n", test_species, "COMPLETE\n")
    
    # Clean up to free memory
    rm(detection_matrix, species_data, site_history, sites_with_detections)
    gc()  # Force garbage collection
    
  }, error = function(e) {
    cat("\nERROR processing", test_species, ":", e$message, "\n")
    cat("Continuing to next species...\n")
  })
  
}

cat("\n========================================\n")
cat("ALL DETECTION MATRICES COMPLETE\n")
cat("========================================\n\n")

cat("Created matrices for", length(all_detection_matrices), "species:\n")
for (sp in names(all_detection_matrices)) {
  cat("  ✓", sp, "\n")
}