##############################
#
# eBird Detection Matrix Creation
# Ian Becker
# August 2026
#
##############################

# Builds detection matrices from eBird zerofilled data for occupancy modeling
# Key difference from GBIF: eBird has explicit non-detections via zerofill
# Monthly detection occasion = 1 if any checklist at site detected species,
#                              0 if all checklists recorded non-detection,
#                              NA if no checklists at that site-month
# All outputs saved to detection_matrices_ebird/

library(dplyr)

setwd("~/Desktop/project_code/UGS_occ/data")

# ============================================================================
# 1. SETUP OUTPUT DIRECTORY
# ============================================================================

if (!dir.exists("detection_matrices_ebird")) {
  dir.create("detection_matrices_ebird")
  cat("Created detection_matrices_ebird/ directory\n\n")
} else {
  cat("Using existing detection_matrices_ebird/ directory\n\n")
}

# ============================================================================
# 2. LOAD EBIRD DATA
# ============================================================================

cat("Loading eBird filtered data...\n")

ebird <- read.csv("ebird_filtered_lrgv_with_sites.csv")

cat("Loaded", nrow(ebird), "checklist-species records\n")
cat("Sites represented:", n_distinct(ebird$site_id), "\n\n")

# ============================================================================
# 3. IDENTIFY EBIRD SITES
# ============================================================================

cat("Identifying eBird sites...\n")

# Use sites with at least one checklist
ebird_site_ids <- sort(unique(ebird$site_id[!is.na(ebird$site_id)]))
n_sites        <- length(ebird_site_ids)

cat("Modeling over", n_sites, "eBird-observed sites\n")
cat("Site IDs:", paste(ebird_site_ids, collapse = ", "), "\n\n")

# ============================================================================
# 4. LOAD SPECIES LIST
# ============================================================================

cat("Loading species list...\n")

species_csv       <- read.csv("lrgv_ugs_species_FILTERED.csv")
test_species_list <- species_csv$common_name

cat("Species to process:", length(test_species_list), "\n")
for (sp in test_species_list) cat("  -", sp, "\n")
cat("\n")

# ============================================================================
# 5. PREPARE TEMPORAL DATA
# ============================================================================

ebird <- ebird %>%
  mutate(
    observation_date = as.Date(observation_date),
    year             = as.integer(format(observation_date, "%Y")),
    month            = as.integer(format(observation_date, "%m"))
  ) %>%
  filter(year >= 2015, year <= 2025,
         !is.na(year), !is.na(month))

cat("Filtered to 2015-2025 study period\n")
cat("Total records:", nrow(ebird), "\n\n")

years    <- 2015:2025
n_years  <- length(years)
months   <- 1:12
n_months <- length(months)

cat("Temporal structure:\n")
cat("  Primary periods (years):", n_years, "\n")
cat("  Secondary occasions (months):", n_months, "\n\n")

# ============================================================================
# 6. PRE-COLLAPSE TO SITE-YEAR-MONTH LEVEL
# ============================================================================

# For each species, site, year, month:
# detected = 1 if any checklist detected species
# detected = 0 if checklists exist but none detected species
# NA if no checklists at that site-month

cat("Pre-collapsing to site-year-month level...\n")

ebird_collapsed <- ebird %>%
  group_by(common_name, site_id, year, month) %>%
  summarise(
    detected         = as.integer(any(species_observed)),
    n_checklists     = n(),
    mean_duration    = mean(duration_minutes,   na.rm = TRUE),
    mean_distance    = mean(effort_distance_km, na.rm = TRUE),
    mean_observers   = mean(number_observers,   na.rm = TRUE),
    .groups          = "drop"
  )

cat("Collapsed records:", nrow(ebird_collapsed), "\n\n")

# ============================================================================
# 7. BUILD DETECTION MATRICES — LOOP THROUGH SPECIES
# ============================================================================

all_detection_matrices <- list()
all_metadata           <- list()

for (test_species in test_species_list) {
  
  tryCatch({
    
    cat("\n========================================\n")
    cat("PROCESSING:", test_species, "\n")
    cat("========================================\n\n")
    
    species_data <- ebird_collapsed %>%
      filter(common_name == test_species,
             site_id %in% ebird_site_ids)
    
    cat("Site-month records:", nrow(species_data), "\n")
    cat("Sites with any detection:",
        n_distinct(species_data$site_id[species_data$detected == 1]), "\n\n")
    
    if (nrow(species_data) == 0) {
      cat("No records found - skipping\n")
      next
    }
    
    cat("Creating detection matrix [", n_sites, "sites x",
        n_years, "years x", n_months, "months]...\n")
    
    # Initialize with NA — NA means no checklists submitted at that site-month
    # 0 = checklists submitted but species not detected
    # 1 = species detected on at least one checklist
    detection_matrix <- array(NA,
                              dim = c(n_sites, n_years, n_months),
                              dimnames = list(
                                site  = ebird_site_ids,
                                year  = years,
                                month = months
                              ))
    
    for (i in 1:nrow(species_data)) {
      
      obs <- species_data[i, ]
      
      site_idx  <- which(ebird_site_ids == obs$site_id)
      year_idx  <- which(years == obs$year)
      month_idx <- obs$month
      
      if (length(site_idx) == 0) next
      if (length(year_idx) == 0) next
      if (is.na(month_idx) || month_idx < 1 || month_idx > 12) next
      
      detection_matrix[site_idx, year_idx, month_idx] <- obs$detected
    }
    
    cat("\nDetection matrix complete!\n\n")
    
    sites_with_detections <- apply(detection_matrix, 1,
                                   function(x) any(x == 1, na.rm = TRUE))
    n_sites_detected      <- sum(sites_with_detections)
    naive_occupancy       <- n_sites_detected / n_sites
    
    cat("Sites with detections:", n_sites_detected, "of", n_sites,
        "(", round(n_sites_detected / n_sites * 100, 1), "%)\n")
    cat("Naive occupancy:", round(naive_occupancy, 3), "\n\n")
    
    detections_by_year <- apply(detection_matrix, 2,
                                function(x) sum(x == 1, na.rm = TRUE))
    cat("Detections by year (site-months):\n")
    for (i in 1:length(years)) {
      cat("  ", years[i], ":", detections_by_year[i], "\n")
    }
    cat("\n")
    
    species_filename <- gsub(" ", "_", test_species)
    
    saveRDS(detection_matrix,
            paste0("detection_matrices_ebird/detection_matrix_",
                   species_filename, ".rds"))
    cat("Saved: detection_matrices_ebird/detection_matrix_",
        species_filename, ".rds\n")
    
    # Save effort dataframe for use in modeling script
    saveRDS(species_data,
            paste0("detection_matrices_ebird/effort_",
                   species_filename, ".rds"))
    cat("Saved: detection_matrices_ebird/effort_",
        species_filename, ".rds\n")
    
    metadata <- list(
      species          = test_species,
      n_sites          = n_sites,
      site_ids         = ebird_site_ids,
      n_years          = n_years,
      n_months         = n_months,
      years            = years,
      months           = months,
      n_sites_detected = n_sites_detected,
      naive_occupancy  = naive_occupancy,
      data_source      = "eBird",
      date_created     = Sys.Date()
    )
    
    saveRDS(metadata,
            paste0("detection_matrices_ebird/detection_metadata_",
                   species_filename, ".rds"))
    cat("Saved: detection_matrices_ebird/detection_metadata_",
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

# ============================================================================
# 8. FINAL SUMMARY
# ============================================================================

cat("\n========================================\n")
cat("ALL DETECTION MATRICES COMPLETE\n")
cat("========================================\n\n")

cat("Created matrices for", length(all_detection_matrices), "species:\n")
for (sp in names(all_detection_matrices)) {
  cat("  v", sp, "\n")
}
