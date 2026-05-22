##############################
#
# Study Area Map
# Ian Becker
# December 2025
#
##############################

library(tigris)
library(ggplot2)
library(sf)
library(dplyr)
library(ggspatial)
options(tigris_use_cache = TRUE)

setwd("~/Desktop/project_code/UGS_occ/data")

output_dir <- "/Users/ianbecker/Desktop/project_code/UGS_occ/figures_tables"

####################
#   Load Data
####################

# load in green space shapefile
ugs <- st_read("lrgv_green_spaces_combined")

# load urban areas boundary
census_urban <- st_read("tl_2020_us_uac20")

# Load Texas county boundaries
tx_counties <- counties(state = "TX", cb = TRUE) %>%
  st_transform(st_crs(ugs))

# LRGV counties in Texas
lrgv_counties <- c("Starr", "Hidalgo", "Cameron", "Willacy")

# Filter to LRGV
lrgv <- tx_counties %>%
  filter(NAME %in% lrgv_counties)

####################
#   Prep data
####################

# Transform census urban to match CRS
census_urban <- st_transform(census_urban, st_crs(ugs))

# Clip census urban areas to LRGV boundary
lrgv_boundary <- st_union(lrgv)
census_urban_lrgv <- st_intersection(census_urban, lrgv_boundary)

# Get centroids of green spaces for point locations
ugs_centroids <- st_centroid(ugs)

####################
#   Plot
####################

# Create the map
study_area_map <- ggplot() +
  
  # County boundaries
  geom_sf(data = lrgv, fill = "gray95", color = "gray40", linewidth = 0.5) +
  
  # Census urban areas (background polygons)
  geom_sf(data = census_urban_lrgv, fill = "#D4B996", color = NA, alpha = 0.6) +
  
  # Study sites as points
  geom_sf(data = ugs_centroids, color = "darkgreen", size = 1.5, alpha = 0.7) +
  
  # County labels
  geom_sf_text(data = lrgv, aes(label = NAME), size = 3.5, color = "gray30") +
  
  # Scale bar
  annotation_scale(
    location   = "bl",
    width_hint = 0.25,
    style      = "ticks",
    text_cex   = 0.8
  ) +
  
  # North arrow
  annotation_north_arrow(
    location    = "tr",
    which_north = "true",
    height      = unit(1, "cm"),
    width       = unit(1, "cm"),
    style       = north_arrow_fancy_orienteering()
  ) +
  
  # Clean theme
  theme_minimal() +
  theme(
    panel.grid  = element_blank(),
    axis.text   = element_text(size = 13),
    axis.title = element_blank()
  ) 

print(study_area_map)

# Save
ggsave(
  file.path(output_dir, "study_area_map.png"),
  study_area_map, width = 10, height = 8, dpi = 300, bg = "white"
)

####################
#   Filter for any species
####################

# ── USER INPUT — change species name here ─────────────────────────────────────
focal_species <- "Roseate Spoonbill"

# Load iNat observations
inat <- read.csv("inat_observations_with_sites.csv")

# Get site IDs with focal species detections
species_sites <- inat %>%
  filter(common_name == focal_species) %>%
  pull(site_id) %>%
  unique()

cat("Sites with", focal_species, "detections:", length(species_sites), "\n")

# Flag detected vs non-detected sites
ugs_centroids_sp <- ugs_centroids %>%
  mutate(detected = site_id %in% species_sites)

# Create species map — same as base map but with detected sites highlighted
species_map <- ggplot() +
  
  geom_sf(data = lrgv, fill = "gray95", color = "gray40", linewidth = 0.5) +
  
  geom_sf(data = census_urban_lrgv, fill = "#D4B996", color = NA, alpha = 0.6) +
  
  # Non-detected sites — grey
 # geom_sf(
  #  data  = ugs_centroids_sp %>% filter(!detected),
   # color = "gray60", size = 1.5, alpha = 0.6
#  ) +
  
  # Detected sites — yellow circle with black border matching exemplar map
  geom_sf(
    data   = ugs_centroids_sp %>% filter(detected),
    shape  = 21,
    fill   = "yellow",
    color  = "black",
    size   = 3,
    stroke = 0.6,
    alpha  = 0.85
  ) +
  
  geom_sf_text(data = lrgv, aes(label = NAME), size = 3.5, color = "gray30") +
  
  annotation_scale(
    location   = "bl",
    width_hint = 0.25,
    style      = "ticks",
    text_cex   = 0.8
  ) +
  
  annotation_north_arrow(
    location    = "tr",
    which_north = "true",
    height      = unit(1, "cm"),
    width       = unit(1, "cm"),
    style       = north_arrow_fancy_orienteering()
  ) +
  
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text  = element_text(size = 8),
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5)
  ) +
  
  labs(x = "Longitude", y = "Latitude")

print(species_map)

# Save
species_filename <- gsub(" ", "_", focal_species)

ggsave(
  file.path(output_dir, paste0("study_area_map_", species_filename, ".png")),
  species_map, width = 10, height = 8, dpi = 300, bg = "white"
)

cat("Saved: study_area_map_", species_filename, ".png\n")
