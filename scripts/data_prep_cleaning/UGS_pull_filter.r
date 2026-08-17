##############################
#
# Green Space Pull and Filtering 
# February 2026
# Ian Becker
#
##############################

# This script manually pulls green spaces within the valley from OpenStreetMap
# and then manually filters through all pulled green space polygons after mapping.

library(sf)
library(osmdata)
library(tigris)
library(mapview)
library(dplyr)

options(tigris_use_cache = TRUE)
setwd("PATH HERE")

#===================================================================
# DEFINE STUDY AREA
#===================================================================

# Define LRGV counties 

lrgv_counties <- c("Starr", "Hidalgo", "Cameron", "Willacy")

# Get Texas counties

tx_counties <- counties(state = "TX", cb = TRUE, year = 2020)

# Filter to LRGV

lrgv <- tx_counties %>%
  filter(NAME %in% lrgv_counties)

cat("Loaded", nrow(lrgv), "LRGV counties:", paste(lrgv$NAME, collapse = ", "), "\n")

# Load Census urban areas (2020)

urban_areas_all <- st_read("tl_2020_us_uac20")

# Filter to urban areas intersecting LRGV

lrgv_urban <- urban_areas_all[lrgv, ]

# Transform to WGS84 for OSM query

lrgv_urban_wgs84 <- st_transform(lrgv_urban, crs = 4326)

# Get bounding box for OSM query

bbox <- st_bbox(lrgv_urban_wgs84)

cat("\nBounding box for OSM query:\n")
cat("  West:", bbox["xmin"], "\n")
cat("  South:", bbox["ymin"], "\n")
cat("  East:", bbox["xmax"], "\n")
cat("  North:", bbox["ymax"], "\n\n")

#===================================================================
# 1. PULL OSM GREEN SPACES
#===================================================================

# Query OSM for green spaces
# Multiple queries for different green space types

tryCatch({
  
  # Parks and gardens
  
  cat("Pulling parks and gardens...\n")
  parks <- opq(bbox = bbox) %>%
    add_osm_feature(key = "leisure", value = c("park", "garden")) %>%
    osmdata_sf()
  
  # Golf courses
  
  cat("Pulling golf courses...\n")
  golf <- opq(bbox = bbox) %>%
    add_osm_feature(key = "leisure", value = "golf_course") %>%
    osmdata_sf()
  
  # Cemeteries
  
  cat("Pulling cemeteries...\n")
  cemeteries <- opq(bbox = bbox) %>%
    add_osm_feature(key = "landuse", value = "cemetery") %>%
    osmdata_sf()
  
  # Recreation grounds
  
  cat("Pulling recreation grounds...\n")
  recreation <- opq(bbox = bbox) %>%
    add_osm_feature(key = "leisure", value = "recreation_ground") %>%
    osmdata_sf()
  
  cat("\nOSM queries complete!\n\n")
  
}, error = function(e) {
  cat("ERROR querying OSM:", e$message, "\n")
  stop("OSM query failed. Check internet connection or try again later.")
})

#===================================================================
# 2. PROCESS PULLED OSM GREEN SPACES
#===================================================================

# Extract polygons only (multipolygons and polygons)

parks_poly <- parks$osm_polygons
golf_poly <- golf$osm_polygons
cemeteries_poly <- cemeteries$osm_polygons
recreation_poly <- recreation$osm_polygons

# Remove NULL results

all_polys <- list(parks_poly, golf_poly, cemeteries_poly, recreation_poly)
all_polys <- all_polys[!sapply(all_polys, is.null)]

if (length(all_polys) == 0) {
  stop("ERROR: No green space polygons found in OSM query!")
}

# Function to standardize columns

standardize_osm <- function(osm_sf, type_label) {
  if (is.null(osm_sf)) return(NULL)
  
  osm_sf %>%
    mutate(
      type = type_label,
      osm_id = osm_id,
      name = if("name" %in% names(.)) name else NA_character_
    ) %>%
    select(osm_id, name, type, geometry)
}

# Standardize and combine

parks_std <- standardize_osm(parks_poly, "park")
golf_std <- standardize_osm(golf_poly, "golf_course")
cemeteries_std <- standardize_osm(cemeteries_poly, "cemetery")
recreation_std <- standardize_osm(recreation_poly, "recreation_ground")

# Combine all

all_green_spaces <- bind_rows(
  parks_std,
  golf_std,
  cemeteries_std,
  recreation_std
)

cat("Combined", nrow(all_green_spaces), "green space polygons\n\n")

#===================================================================
# FILTER OSM GREEN SPACES
#===================================================================

# Transform to common CRS (UTM Zone 14N)

all_green_spaces <- st_transform(all_green_spaces, crs = 32614)
lrgv_urban_utm <- st_transform(lrgv_urban, crs = 32614)

# Filter to sites within urban areas

green_spaces_urban <- st_intersection(all_green_spaces, 
                                      st_union(lrgv_urban_utm))

cat("  Retained", nrow(green_spaces_urban), "sites within urban boundaries\n")

# Calculate area

green_spaces_urban$area_m2 <- st_area(green_spaces_urban)
green_spaces_urban$area_ha <- as.numeric(green_spaces_urban$area_m2) / 10000

# Filter by size (1-100 ha)

green_spaces_filtered <- green_spaces_urban %>%
  filter(area_ha >= 1, area_ha <= 100)

cat("  Retained", nrow(green_spaces_filtered), "sites between 1-100 ha\n")

# Filter to named sites only (excludes small unlabeled patches)

green_spaces_named <- green_spaces_filtered %>%
  filter(!is.na(name) & name != "")

cat("  Retained", nrow(green_spaces_named), "named sites\n")

# Add site_id

osm_sites <- green_spaces_named %>%
  mutate(site_id = row_number()) %>%
  select(site_id, osm_id, name, type, area_ha, geometry)

# Save initial OSM extract

st_write(osm_sites, "lrgv_osm_green_spaces.shp", delete_dsn = TRUE)
cat("\nSaved: lrgv_osm_green_spaces.shp\n")

# Load in green spaces

osm_sites <- st_read("lrgv_osm_green_spaces")

# Interactive map with labels

m <- mapview(osm_sites, 
             zcol = "type",
             label = paste0("ID:", osm_sites$site_id, " - ", 
                            osm_sites$name, " (", 
                            round(osm_sites$area_ha, 1), " ha)"),
             layer.name = "Site Type")

# Add sites to remove

sites_to_remove <- c(132, 134, 154, 187, 197, 157, 155, 156, 202, 198, 201, 200, 203, 199)  

osm_final <- osm_sites %>%
  filter(!site_id %in% sites_to_remove)

st_write(osm_final, "lrgv_osm_green_spaces_filtered.shp", delete_dsn = TRUE)

