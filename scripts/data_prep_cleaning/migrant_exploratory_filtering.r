##############################
#
# Migrants exploratory data prep
# Ian Becker
# February 2026
#
##############################

# This script checks the migration season window for selected species then 
# creates detection histories to be used for occupancy modelling

library(dplyr)
library(lubridate)
library(purrr)
library(tidyverse)
library(sf)

setwd("~/Desktop/project_code/UGS_occ/data")

#===================================================================
# CALCULATING MIGRATION WINDOWS (95% CI AROUND DATA)
#===================================================================

# load site data

inat_with_sites <- read.csv("inat_observations_with_sites.csv")

# Function to calculate migration windows

calculate_migration_windows <- function(species_name, inat_data) {
  
  species_obs <- inat_data %>%
    filter(common_name == species_name) %>%
    mutate(
      # Ensure date is properly formatted
      date = as.Date(observed_on),
      year = year(date),
      month = month(date),  # Extract month as numeric
      doy = yday(date),     # Day of year (1-365)
      # Define seasons
      season = case_when(
        month %in% c(3, 4, 5, 6) ~ "spring",
        month %in% c(8, 9, 10, 11) ~ "fall",
        TRUE ~ "other"
      )
    ) %>%
    filter(season %in% c("spring", "fall"))
  
  # Check if we have any data
  if (nrow(species_obs) == 0) {
    cat("Warning: No spring/fall observations for", species_name, "\n")
    return(NULL)
  }
  
  # Calculate 95% quantiles for each season
  windows <- species_obs %>%
    group_by(season) %>%
    summarise(
      n_obs = n(),
      start_doy = round(quantile(doy, 0.025)),  # 2.5th percentile
      end_doy = round(quantile(doy, 0.975)),    # 97.5th percentile
      peak_doy = round(median(doy)),
      n_sites = n_distinct(site_id),
      .groups = "drop"
    ) %>%
    mutate(
      duration_days = end_doy - start_doy,
      # Convert DOY to actual dates for easier interpretation
      start_date = as.Date(start_doy - 1, origin = "2020-01-01"),
      end_date = as.Date(end_doy - 1, origin = "2020-01-01"),
      peak_date = as.Date(peak_doy - 1, origin = "2020-01-01"),
      start_month = month(start_date),
      end_month = month(end_date)
    )
  
  return(list(
    species = species_name,
    windows = windows,
    raw_obs_count = nrow(species_obs)
  ))
}

# Test with passage migrants

passage_migrants <- c(
  "Black-throated Green Warbler",
  "Tennessee Warbler",
  "Summer Tanager",
  "Baltimore Oriole",
  "Black-and-white Warbler"
)

cat("=== CALCULATING MIGRATION WINDOWS ===\n\n")

migration_windows <- list()
for (sp in passage_migrants) {
  cat("Processing:", sp, "... ")
  
  result <- calculate_migration_windows(sp, inat_with_sites)
  
  if (!is.null(result)) {
    migration_windows[[sp]] <- result
    cat("✓\n")
  } else {
    cat("✗ (no data)\n")
  }
}

cat("\n=== MIGRATION WINDOW RESULTS ===\n\n")

# View results in a nice format

for (sp in names(migration_windows)) {
  cat("\n", sp, ":\n")
  cat(strrep("-", nchar(sp) + 3), "\n")
  
  windows <- migration_windows[[sp]]$windows
  
  for (i in 1:nrow(windows)) {
    row <- windows[i, ]
    cat(sprintf("  %s: %s to %s (peak: %s)\n",
                toupper(row$season),
                format(row$start_date, "%b %d"),
                format(row$end_date, "%b %d"),
                format(row$peak_date, "%b %d")))
    cat(sprintf("    Duration: %d days | n = %d obs | %d sites\n",
                row$duration_days, row$n_obs, row$n_sites))
  }
}

# Summary table

cat("\n\n=== SAMPLE SIZE SUMMARY ===\n\n")

summary_table <- map_dfr(migration_windows, function(sp_data) {
  sp_data$windows %>%
    mutate(
      species = sp_data$species,
      months = paste(month.abb[start_month], "-", month.abb[end_month]),
      viable = n_obs >= 25 & n_sites >= 10
    ) %>%
    select(species, season, months, n_obs, n_sites, duration_days, viable)
})

print(summary_table)

cat("\n")
cat("Viable for modeling (≥25 obs, ≥10 sites):\n")
viable <- summary_table %>% filter(viable)
if (nrow(viable) > 0) {
  for (i in 1:nrow(viable)) {
    cat(sprintf("  ✓ %s (%s): %d obs across %d sites\n",
                viable$species[i], viable$season[i], 
                viable$n_obs[i], viable$n_sites[i]))
  }
} else {
  cat("  None - consider lowering thresholds or focusing on winterers\n")
}

#===================================================================
# CREATING DETECTION HISTORIES FOR SELECTED MIGRANTS/SEASONS
#===================================================================

# Ensure temporal variables exist

if (!"month" %in% names(inat_with_sites)) {
  cat("Adding temporal variables...\n")
  inat_with_sites <- inat_with_sites %>%
    mutate(
      date = as.Date(observed_on),
      year = year(date),
      month = month(date)
    )
}

# Load sites

sites <- st_read("lrgv_green_spaces_detection_filtered")
n_sites <- nrow(sites)

# Species-specific windows based on 95% quantiles

spring_migrant_windows <- list(
  "Summer Tanager" = list(
    months = c(3, 4, 5),  # March-May
    window = "Mar 15 - May 9",
    n_obs = 173,
    n_sites = 10
  ),
  "Baltimore Oriole" = list(
    months = c(4, 5),  # April-May ONLY
    window = "Apr 12 - May 4", 
    n_obs = 244,
    n_sites = 10
  ),
  "Black-and-white Warbler" = list(
    months = c(3, 4, 5),  # March-May
    window = "Mar 15 - May 8",
    n_obs = 147,
    n_sites = 10
  )
)

cat("Species-specific windows:\n")
for (sp in names(spring_migrant_windows)) {
  info <- spring_migrant_windows[[sp]]
  cat(sprintf("  %-25s: %s (%d months)\n", 
              sp, info$window, length(info$months)))
}
cat("\n")

####################
#   Build Detection Matrices
####################

for (test_species in names(spring_migrant_windows)) {
  
  cat("----------------------------------------\n")
  cat("PROCESSING:", test_species, "\n")
  cat("----------------------------------------\n")
  
  # Get species-specific info
  sp_info <- spring_migrant_windows[[test_species]]
  spring_months <- sp_info$months
  n_months <- length(spring_months)
  
  cat("  Window:", sp_info$window, "\n")
  cat("  Months:", paste(month.abb[spring_months], collapse = ", "), "\n\n")
  
  tryCatch({
    
    # Filter to species-specific months
    species_data <- inat_with_sites %>%
      filter(
        common_name == test_species,
        month %in% spring_months,
        year >= 2015, year <= 2025
      )
    
    cat("  Observations in window:", nrow(species_data), "\n")
    
    if (nrow(species_data) == 0) {
      cat("  No observations - skipping\n\n")
      next
    }
    
    # Define temporal structure
    years <- 2015:2025
    n_years <- length(years)
    
    # Initialize detection matrix
    detection_matrix <- array(0, 
                              dim = c(n_sites, n_years, n_months),
                              dimnames = list(
                                site = 1:n_sites,
                                year = years,
                                month = month.abb[spring_months]
                              ))
    
    # Fill in detections
    for (i in 1:nrow(species_data)) {
      
      if (i %% 50 == 0) {
        cat("    Processing observation", i, "of", nrow(species_data), "\n")
      }
      
      obs <- species_data[i, ]
      
      site_idx <- obs$site_id
      year_idx <- which(years == obs$year)
      month_idx <- which(spring_months == obs$month)
      
      # Validate indices
      if (!is.na(site_idx) && length(year_idx) > 0 && length(month_idx) > 0) {
        if (site_idx >= 1 && site_idx <= n_sites) {
          detection_matrix[site_idx, year_idx, month_idx] <- 1
        }
      }
    }
    
    # Calculate summary stats
    sites_with_detections <- apply(detection_matrix, 1, function(x) any(x == 1))
    n_sites_detected <- sum(sites_with_detections)
    naive_occupancy <- n_sites_detected / n_sites
    
    # Detections by year
    detections_by_year <- apply(detection_matrix, 2, sum)
    
    cat("\n  === SUMMARY ===\n")
    cat("  Sites with detections:", n_sites_detected, "of", n_sites, 
        "(", round(n_sites_detected/n_sites*100, 1), "%)\n")
    cat("  Naive occupancy:", round(naive_occupancy, 3), "\n\n")
    
    cat("  Detections by year (site-months):\n")
    for (i in 1:length(years)) {
      cat("    ", years[i], ":", detections_by_year[i], "\n")
    }
    
    # Show example detection histories
    cat("\n  Example detection histories:\n")
    sites_to_show <- which(sites_with_detections)[1:min(2, n_sites_detected)]
    
    for (site_idx in sites_to_show) {
      cat("    Site", site_idx, ":\n")
      site_history <- detection_matrix[site_idx, , ]
      
      for (y in 1:min(3, n_years)) {
        detections <- paste(site_history[y, ], collapse = "")
        n_det <- sum(site_history[y, ])
        cat("      ", years[y], ":", detections, "(", n_det, "months)\n")
      }
      if (n_years > 3) cat("      ...\n")
    }
    
    # Metadata
    metadata <- list(
      species = test_species,
      season = "spring",
      migration_type = "passage",
      months = spring_months,
      month_names = month.abb[spring_months],
      phenology_window = sp_info$window,
      years = years,
      n_sites = n_sites,
      n_years = n_years,
      n_months = n_months,
      n_observations = nrow(species_data),
      n_sites_detected = n_sites_detected,
      naive_occupancy = naive_occupancy,
      date_created = Sys.Date()
    )
    
    # Save
    species_filename <- gsub(" ", "_", test_species)
    saveRDS(detection_matrix, 
            paste0("detection_matrices/detection_matrix_spring_", 
                   species_filename, ".rds"))
    saveRDS(metadata, 
            paste0("detection_matrices/detection_metadata_spring_", 
                   species_filename, ".rds"))
    
    cat("\n  Saved: detection_matrix_spring_", species_filename, ".rds\n")
    cat("  Saved: detection_metadata_spring_", species_filename, ".rds\n\n")
    
    # Clean up
    rm(detection_matrix, species_data, site_history, sites_with_detections)
    gc(verbose = FALSE)
    
  }, error = function(e) {
    cat("  ✗ ERROR:", e$message, "\n\n")
  })
}

cat("========================================\n")
cat("SPRING MIGRANT MATRICES COMPLETE\n")
cat("========================================\n\n")
