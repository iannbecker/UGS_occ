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
    location = "bl",           # bottom-left
    width_hint = 0.25,         # proportion of plot width
    style = "ticks",           # "bar" or "ticks"
    text_cex = 0.8             # text size
  ) +
  
  # North arrow
  annotation_north_arrow(
    location = "tr",           # top-right
    which_north = "true",
    height = unit(1, "cm"),
    width = unit(1, "cm"),
    style = north_arrow_fancy_orienteering()  # or north_arrow_minimal()
  ) +
  
  
  # Clean theme
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text = element_text(size = 8),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5, color = "gray50")
  ) +
  labs(
    x = "Longitude", y = "Latitude"
  )

print(study_area_map)

# Save
ggsave("study_area_map_pre_filter.png", plot = study_area_map, 
       width = 10, height = 8, dpi = 300, bg = "white")



