##############################
#
# GBIF Detection Matrix Creation
# Ian Becker
# May 2026
#
##############################

# Builds detection matrices from GBIF data for occupancy modeling
# Uses same species list as original iNat analysis
# All outputs saved to _gbif directories to avoid overwriting iNat results

library(sf)
library(dplyr)
library(lubridate)

setwd("~/Desktop/project_code/UGS_occ/data")

####################
#   Setup Output Directory
####################

if (!dir.exists("detection_matrices_gbif")) {
  dir.create("detection_matrices_gbif")
  cat("Created detection_matrices_gbif/ directory\n\n")
} else {
  cat("Using existing detection_matrices_gbif/ directory\n\n")
}

####################
#   Load GBIF Data
####################

cat("Loading GBIF filtered data...\n")

gbif_with_sites <- read.csv("gbif_in_sites.csv")

cat("Loaded", nrow(gbif_with_sites), "observations\n")
cat("Sites represented:", n_distinct(gbif_with_sites$site_id), "\n\n")

####################
#   Load Sites — Filter to GBIF Observed Sites
####################

cat("Loading site covariates to identify GBIF sites...\n")

site_covs_gbif <- readRDS("site_covariates_gbif.rds")
gbif_site_ids  <- sort(unique(site_covs_gbif$site_id))
n_sites        <- length(gbif_site_ids)

cat("Modeling over", n_sites, "GBIF-observed sites\n")
cat("Site IDs:", paste(gbif_site_ids, collapse = ", "), "\n\n")

####################
#   Load Species List
####################

cat("Loading species list...\n")

species_csv       <- read.csv("lrgv_ugs_species_FILTERED.csv")
test_species_list <- species_csv$common_name

cat("Species to process:", length(test_species_list), "\n")
for (sp in test_species_list) cat("  -", sp, "\n")
cat("\n")

####################
#   Prepare Temporal Data
####################

# GBIF uses separate year/month/day columns — confirm and filter
gbif_with_sites <- gbif_with_sites %>%
  filter(year >= 2015, year <= 2025,
         !is.na(year), !is.na(month))

cat("Filtered to 2015-2025 study period\n")
cat("Total observations:", nrow(gbif_with_sites), "\n\n")

####################
#   Define Temporal Structure
####################

years    <- 2015:2025
n_years  <- length(years)
months   <- 1:12
n_months <- length(months)

cat("Temporal structure:\n")
cat("  Primary periods (years):", n_years, "\n")
cat("  Secondary occasions (months):", n_months, "\n\n")

####################
#   Loop Through Species
####################

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
    
    ####################
    #   Create Detection Matrix
    ####################
    
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
    
    ####################
    #   Matrix Summary
    ####################
    
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
    
    ####################
    #   Save Results
    ####################
    
    species_filename <- gsub(" ", "_", test_species)
    
    saveRDS(detection_matrix,
            paste0("detection_matrices_gbif/detection_matrix_",
                   species_filename, ".rds"))
    cat("Saved: detection_matrices_gbif/detection_matrix_",
        species_filename, ".rds\n")
    
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
    
    all_detection_matrices[[test_species]] <- detection_matrix
    all_metadata[[test_species]]           <- metadata
    
    cat("\n", test_species, "COMPLETE\n")
    
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