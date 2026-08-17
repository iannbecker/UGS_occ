##############################
#
# iNat Detection Matrix Creation
# Ian Becker
# May 2026
#
##############################

# Builds detection matrices from GBIF (iNat) data for 
# landscape level occupancy modelling 

library(sf)
library(dplyr)
library(lubridate)

# ============================================================================
# 1. LOAD DATA AND PREP DIRECTORIES
# ============================================================================

setwd("~/Desktop/project_code/UGS_occ/data")

# create detection matrix output directory

if (!dir.exists("detection_matrices_gbif")) {
  dir.create("detection_matrices_gbif")
  cat("Created detection_matrices_gbif/ directory\n\n")
} else {
  cat("Using existing detection_matrices_gbif/ directory\n\n")
}

# Load in filtered iNat data

gbif_with_sites <- read.csv("gbif_in_sites.csv")

# check numbers 

nrow(gbif_with_sites)
n_distinct(gbif_with_sites$site_id)

# load in site covariates 

site_covs_gbif <- readRDS("site_covariates_gbif.rds")
gbif_site_ids  <- sort(unique(site_covs_gbif$site_id))

#  check that numbers match

n_sites <- length(gbif_site_ids)

# load in species list

species_csv       <- read.csv("lrgv_ugs_species_FILTERED.csv")
test_species_list <- species_csv$common_name

# Check species list 

length(test_species_list)
for (sp in test_species_list) cat("  -", sp, "\n")

# ============================================================================
# 2. PREP TEMPORAL DATA
# ============================================================================

# Not necessary but sanity check filter to 2015-2025 in case there are any outliers

gbif_with_sites <- gbif_with_sites %>%
  filter(year >= 2015, year <= 2025,
         !is.na(year), !is.na(month))


nrow(gbif_with_sites)

# Setup temporal grain for detection matrices

years    <- 2015:2025
n_years  <- length(years)
months   <- 1:12
n_months <- length(months)

# ============================================================================
# 3. LOOP THROUGH SPECIES CREATING DETECTION MATRICES
# ============================================================================

all_detection_matrices <- list()
all_metadata           <- list()

for (test_species in test_species_list) {
  
  tryCatch({
    
    cat("\n========================================\n")
    cat("PROCESSING:", test_species, "\n")
    cat("========================================\n\n")
    
    # Filter to species
    
    species_data <- gbif_with_sites %>%
      filter(common_name == test_species)
    
    cat("Total observations:", nrow(species_data), "\n")
    cat("Sites with observations:", n_distinct(species_data$site_id), "\n\n")
    
    if (nrow(species_data) == 0) {
      cat("No observations found - skipping\n")
      next
    }
    
  #### Create detection matrix
    
    cat("Creating detection matrix [", n_sites, "sites x",
        n_years, "years x", n_months, "months]...\n")
    
    detection_matrix <- array(0,
                              dim = c(n_sites, n_years, n_months),
                              dimnames = list(
                                site  = gbif_site_ids,
                                year  = years,
                                month = months
                              ))
    
    for (i in 1:nrow(species_data)) {
      
      if (i %% 100 == 0) cat("  Processing observation", i,
                             "of", nrow(species_data), "\n")
      
      obs <- species_data[i, ]
      
      # Map original site_id to 1:n_sites index
      
      site_idx  <- which(gbif_site_ids == obs$site_id)
      year_idx  <- which(years == obs$year)
      month_idx <- obs$month
      
      if (length(site_idx) == 0) next
      if (length(year_idx) == 0) next
      if (is.na(month_idx) || month_idx < 1 || month_idx > 12) next
      
      detection_matrix[site_idx, year_idx, month_idx] <- 1
    }
    
    cat("\nDetection matrix complete!\n\n")
    
  ### Summarize
    
    sites_with_detections <- apply(detection_matrix, 1, function(x) any(x == 1))
    n_sites_detected      <- sum(sites_with_detections)
    naive_occupancy       <- n_sites_detected / n_sites
    
    cat("Sites with detections:", n_sites_detected, "of", n_sites,
        "(", round(n_sites_detected / n_sites * 100, 1), "%)\n")
    cat("Naive occupancy:", round(naive_occupancy, 3), "\n\n")
    
    detections_by_year <- apply(detection_matrix, 2, sum)
    cat("Detections by year (site-months):\n")
    for (i in 1:length(years)) {
      cat("  ", years[i], ":", detections_by_year[i], "\n")
    }
    cat("\n")
    
  ### Save results 
    
    species_filename <- gsub(" ", "_", test_species)
    
  # Save detection matrix 
    
    saveRDS(detection_matrix,
            paste0("detection_matrices_gbif/detection_matrix_",
                   species_filename, ".rds"))
    cat("Saved: detection_matrices_gbif/detection_matrix_",
        species_filename, ".rds\n")
    
  # Save associated metadata
    
    metadata <- list(
      species          = test_species,
      n_sites          = n_sites,
      site_ids         = gbif_site_ids,
      n_years          = n_years,
      n_months         = n_months,
      years            = years,
      months           = months,
      n_observations   = nrow(species_data),
      n_sites_detected = n_sites_detected,
      naive_occupancy  = naive_occupancy,
      data_source      = "GBIF",
      date_created     = Sys.Date()
    )
    
    saveRDS(metadata,
            paste0("detection_matrices_gbif/detection_metadata_",
                   species_filename, ".rds"))
    cat("Saved: detection_matrices_gbif/detection_metadata_",
        species_filename, ".rds\n")
 
  # Add to all results list
    
    all_detection_matrices[[test_species]] <- detection_matrix
    all_metadata[[test_species]]           <- metadata
    
    cat("\n", test_species, "COMPLETE\n")
  
  # Cleanup to save space
    
    rm(detection_matrix, species_data, sites_with_detections)
    gc(verbose = FALSE)
    
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
  cat("  v", sp, "\n")
}