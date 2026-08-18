##############################
#
# eBird Detection Matrix Creation — Checklist Subsample
# Ian Becker
# August 2026
#
##############################

# Builds detection matrices for ebird checklist data for landscape-level occupancy modeling

library(dplyr)

setwd("~/Desktop/project_code/UGS_occ/data")

# ============================================================================
# 1. SETUP AND LOAD DATA
# ============================================================================

n_sample <- 10  # checklists subsampled per site-year

# Create output folder

dir.create("detection_matrices_ebird_checklist")

# Load in eBird data

ebird <- read.csv("ebird_filtered_lrgv_with_sites.csv") %>%
  mutate(
    observation_date = as.Date(observation_date),
    year             = as.integer(format(observation_date, "%Y")),
    month            = as.integer(format(observation_date, "%m"))
  ) %>%
  filter(year >= 2015, year <= 2025,
         !is.na(year), !is.na(month),
         !is.na(site_id))

# Number check 

nrow(ebird)
n_distinct(ebird$site_id)

# Keep track of site ids and setup for loop

ebird_site_ids <- sort(unique(ebird$site_id))
n_sites <- length(ebird_site_ids)

# Load in species list

species_csv       <- read.csv("lrgv_ugs_species_FILTERED.csv")
test_species_list <- species_csv$common_name

# Check species list

print(test_species_list)

# Setup temporal grain

years <- 2015:2025
n_years  <- length(years)
n_occ <- n_sample  # secondary occasions = subsampled checklists per site-year

# ============================================================================
# 2. SUBSAMPLE CHECKLISTS PER SITE-YEAR
# ============================================================================

set.seed(42)  

# Get unique checklists 

checklists_all <- ebird %>%
  distinct(checklist_id, site_id, year, month,
           duration_minutes, effort_distance_km, number_observers,
           species_observed, common_name)

# Subsample unique checklist IDs per site-year

sampled_ids <- checklists_all %>%
  distinct(checklist_id, site_id, year) %>%
  group_by(site_id, year) %>%
  group_modify(~ slice_sample(.x, n = min(n_sample, nrow(.x)),
                              replace = FALSE)) %>%
  ungroup()


# Filter full data to sampled checklist IDs only

ebird_sampled <- checklists_all %>%
  semi_join(sampled_ids, by = c("checklist_id", "site_id", "year"))

# ============================================================================
# 3. BUILD EFFORT VARIABLES FOR MODELLING
# ============================================================================

# Assign occasion index (1:n_sample) within each site-year

ebird_sampled <- ebird_sampled %>%
  group_by(site_id, year, checklist_id) %>%
  mutate(occ_idx = cur_group_id()) %>%
  ungroup()

# Map occasion index to 1:n_sample within each site-year

ebird_sampled <- ebird_sampled %>%
  group_by(site_id, year) %>%
  mutate(occ_idx = as.integer(factor(checklist_id,
                                     levels = unique(checklist_id)))) %>%
  ungroup()

# Save effort lookup for modeling script

effort_lookup <- ebird_sampled %>%
  distinct(checklist_id, site_id, year, month, occ_idx,
           duration_minutes, effort_distance_km, number_observers)

saveRDS(effort_lookup,
        "detection_matrices_ebird_checklist/effort_lookup.rds")


# ============================================================================
# 4. LOOP TO BUILD SPECIES DETECTION MATRIX 
# ============================================================================

# Setup output lists

all_detection_matrices <- list()
all_metadata           <- list()

for (test_species in test_species_list) {
  
  tryCatch({
    
    cat("\n========================================\n")
    cat("PROCESSING:", test_species, "\n")
    cat("========================================\n\n")
    
    species_data <- ebird_sampled %>%
      filter(common_name == test_species,
             site_id %in% ebird_site_ids)
    
    if (nrow(species_data) == 0) {
      cat("No records found - skipping\n")
      next
    }
    
    cat("Checklist-level records:", nrow(species_data), "\n")
    cat("Sites with any detection:",
        n_distinct(species_data$site_id[species_data$species_observed]), "\n\n")
    
    cat("Creating detection matrix [", n_sites, "sites x",
        n_years, "years x", n_occ, "occasions]...\n")
    
    # NA = no checklist submitted for that occasion slot
    # 0 = checklist submitted, species not detected
    # 1 = species detected on checklist
    
    detection_matrix <- array(NA,
                              dim = c(n_sites, n_years, n_occ),
                              dimnames = list(
                                site      = ebird_site_ids,
                                year      = years,
                                occasion  = 1:n_occ
                              ))
    
    for (i in 1:nrow(species_data)) {
      
      obs <- species_data[i, ]
      
      site_idx <- which(ebird_site_ids == obs$site_id)
      year_idx <- which(years == obs$year)
      occ_idx  <- obs$occ_idx
      
      if (length(site_idx) == 0) next
      if (length(year_idx) == 0) next
      if (is.na(occ_idx) || occ_idx < 1 || occ_idx > n_occ) next
      
      detection_matrix[site_idx, year_idx, occ_idx] <- as.integer(obs$species_observed)
    }
    
    cat("Detection matrix complete!\n\n")
    
    sites_with_detections <- apply(detection_matrix, 1,
                                   function(x) any(x == 1, na.rm = TRUE))
    n_sites_detected      <- sum(sites_with_detections)
    naive_occupancy       <- n_sites_detected / n_sites
    
    cat("Sites with detections:", n_sites_detected, "of", n_sites,
        "(", round(n_sites_detected / n_sites * 100, 1), "%)\n")
    cat("Naive occupancy:", round(naive_occupancy, 3), "\n\n")
    
    species_filename <- gsub(" ", "_", test_species)
    
    # Save detection matrix
    
    saveRDS(detection_matrix,
            paste0("detection_matrices_ebird_checklist/detection_matrix_",
                   species_filename, ".rds"))
    cat("Saved: detection_matrix_", species_filename, ".rds\n")
    
    # Generate detection matrix metadata
    
    metadata <- list(
      species          = test_species,
      n_sites          = n_sites,
      site_ids         = ebird_site_ids,
      n_years          = n_years,
      n_occ            = n_occ,
      years            = years,
      n_sample         = n_sample,
      n_sites_detected = n_sites_detected,
      naive_occupancy  = naive_occupancy,
      data_source      = "eBird_checklist",
      date_created     = Sys.Date()
    )
    
    # Save metadata
    
    saveRDS(metadata,
            paste0("detection_matrices_ebird_checklist/detection_metadata_",
                   species_filename, ".rds"))
    cat("Saved: detection_metadata_", species_filename, ".rds\n")
    
    all_detection_matrices[[test_species]] <- detection_matrix
    all_metadata[[test_species]]           <- metadata
    
    cat("\n", test_species, "COMPLETE\n")
    
    rm(detection_matrix, species_data, sites_with_detections)
    gc(verbose = FALSE)
    
  }, error = function(e) {
    cat("\nERROR processing", test_species, ":", e$message, "\n")
    cat("Continuing to next species...\n")
  })
  
  if(i == length(test_species_list)) {
    cat("\nAll species processed.\n")
  }
}
