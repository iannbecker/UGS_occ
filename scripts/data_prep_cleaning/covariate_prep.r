##############################
#
# Extract Site-Level Habitat Covariates
# Ian Becker
# December 2025
#
##############################

# Extract land cover data for each urban green space site
# and creates covariate matrix for occupancy modeling

library(sf)
library(terra)
library(dplyr)

# ============================================================================
# 1. LOAD DATA
# ============================================================================

# Load sites

sites <- st_read("lrgv_green_spaces_ebird")
n_sites <- nrow(sites)

cat("Loaded", n_sites, "urban green space sites\n")

# Load land cover raster (Dynamic World Cover)

landcover <- rast("lrgv_dynamic_world_cover.tif")

# Ensure same CRS

if (st_crs(sites)$input != crs(landcover, describe = TRUE)$code) {
  cat("Reprojecting sites to match land cover...\n")
  sites <- st_transform(sites, crs(landcover))
}

# ============================================================================
# 2. CALCULATE SITE AREA 
# ============================================================================

# Calculate areas (meters squared and hectares)

sites$area_m2 <- st_area(sites)
sites$area_ha <- as.numeric(sites$area_m2) / 10000  # Convert m^2 to hectares

# Check for any outliers

n_small <- sum(sites$area_ha < 1)
n_large <- sum(sites$area_ha > 100)

# ============================================================================
# 3. EXTRACT LAND COVER FOR EACH SITE
# ============================================================================

# Set up covariate dataframe

site_covariates <- data.frame(
  site_id = integer(),
  area_ha = numeric(),
  log_area = numeric(),
  trees_pct = numeric(),
  grass_pct = numeric(),
  shrub_pct = numeric(),
  flooded_veg_pct = numeric(),
  urban_pct = numeric(),
  crops_pct = numeric(),
  bare_pct = numeric(),
  water_pct = numeric(),
  habitat_diversity = numeric(),
  stringsAsFactors = FALSE
)

# Extract for each site

for (i in 1:n_sites) {
  
  if (i %% 10 == 0) {
    cat("Processing site", i, "of", n_sites, 
        "(", round(i/n_sites*100, 1), "%)\n")
  }
  
  site <- sites[i, ]
  
  # Extract land cover values within this site polygon
  
  lc_values <- terra::extract(landcover, vect(site), df = TRUE)
  
  if (nrow(lc_values) == 0) {
    
    # No land cover data for this site
    
    cat("  Warning: No land cover data for site", i, "\n")
    
    site_result <- data.frame(
      site_id = site$site_id,
      area_ha = site$area_ha,
      log_area = log(site$area_ha + 0.01),
      trees_pct = 0,
      grass_pct = 0,
      shrub_pct = 0,
      flooded_veg_pct = 0,
      urban_pct = 0,
      crops_pct = 0,
      bare_pct = 0,
      water_pct = 0,
      habitat_diversity = 0
    )
    
  } else {
    
    # Get land cover column
    
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
    # 3 = Flooded vegetation
    # 4 = Crops
    # 5 = Shrub/Scrub
    # 6 = Built area/Urban
    # 7 = Bare ground
    # 8 = Water
    
    trees_pct <- ifelse(1 %in% lc_summary$class, 
                        lc_summary$percentage[lc_summary$class == 1], 0)
    grass_pct <- ifelse(2 %in% lc_summary$class, 
                        lc_summary$percentage[lc_summary$class == 2], 0)
    flooded_veg_pct <- ifelse(3 %in% lc_summary$class,
                              lc_summary$percentage[lc_summary$class == 3], 0)
    crops_pct <- ifelse(4 %in% lc_summary$class, 
                        lc_summary$percentage[lc_summary$class == 4], 0)
    shrub_pct <- ifelse(5 %in% lc_summary$class,
                        lc_summary$percentage[lc_summary$class == 5], 0)
    urban_pct <- ifelse(6 %in% lc_summary$class, 
                        lc_summary$percentage[lc_summary$class == 6], 0)
    bare_pct <- ifelse(7 %in% lc_summary$class, 
                       lc_summary$percentage[lc_summary$class == 7], 0)
    water_pct <- ifelse(0 %in% lc_summary$class, 
                        lc_summary$percentage[lc_summary$class == 0], 0)
    
    # Calculate habitat diversity (Shannon index)
    
    proportions <- lc_summary$percentage / 100
    proportions <- proportions[proportions > 0]
    habitat_diversity <- -sum(proportions * log(proportions), na.rm = TRUE)
    
    site_result <- data.frame(
      site_id = site$site_id,
      area_ha = site$area_ha,
      log_area = log(site$area_ha + 0.01),
      trees_pct = trees_pct,
      grass_pct = grass_pct,
      shrub_pct = shrub_pct,
      flooded_veg_pct = flooded_veg_pct,
      urban_pct = urban_pct,
      crops_pct = crops_pct,
      bare_pct = bare_pct,
      water_pct = water_pct,
      habitat_diversity = habitat_diversity
    )
  }
  
  site_covariates <- rbind(site_covariates, site_result)
  cat("\n=== EXTRACTION COMPLETE ===\n\n")
}

# Calculate summary stats

covar_summary <- site_covariates %>%
  select(-site_id) %>%
  summarise(across(everything(), 
                   list(mean = ~mean(., na.rm = TRUE),
                        sd = ~sd(., na.rm = TRUE),
                        min = ~min(., na.rm = TRUE),
                        max = ~max(., na.rm = TRUE))))

# ============================================================================
# 4. DATA AND CORRELATION CHECKS 
# ============================================================================

# Sites with missing data

sites_no_data <- sum(rowSums(site_covariates[, -1]) == 0)
cat("Sites with no land cover data:", sites_no_data, "\n")

# Check total percentages

total_pct <- rowSums(site_covariates[, c("trees_pct", "grass_pct", "shrub_pct", 
                                         "flooded_veg_pct", "urban_pct", "water_pct",
                                         "crops_pct", "bare_pct")])

round(mean(total_pct), 1)
round(min(total_pct), 1)
round(max(total_pct), 1)

# Sites with outlier totals (not close to 100%)

outlier <- which(total_pct < 90 | total_pct > 110)

# Check for correlations between land cover classes

cor_matrix <- cor(site_covariates[, -1], use = "complete.obs")
print(round(cor_matrix, 2))

# High correlations

high_cor <- which(abs(cor_matrix) > 0.6 & abs(cor_matrix) < 1, arr.ind = TRUE)
if (nrow(high_cor) > 0) {
  cat("High correlations (|r| > 0.6) detected:\n")
  for (i in 1:nrow(high_cor)) {
    var1 <- rownames(cor_matrix)[high_cor[i, 1]]
    var2 <- colnames(cor_matrix)[high_cor[i, 2]]
    cor_val <- cor_matrix[high_cor[i, 1], high_cor[i, 2]]
    if (var1 < var2) {  
      cat("  ", var1, "vs", var2, ":", round(cor_val, 2), "\n")
    }
  }
  cat("\nConsider removing one of each highly correlated pair\n")
}

cat("\n")

# ============================================================================
# 5. SAVE RESULTS
# ============================================================================

# Save as CSV

write.csv(site_covariates, "site_covariates_ebird.csv", row.names = FALSE)


# Save as RDS

saveRDS(site_covariates, "site_covariates_ebird.rds")








