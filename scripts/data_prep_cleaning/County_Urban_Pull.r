##############################
#
# County and Urban Areas data pull
# 12/15/2025
# Ian Becker
#
##############################

library(tigris)
library(dplyr)
library(sf)

options(tigris_use_cache = TRUE)

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
