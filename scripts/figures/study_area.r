##############################
#
# Study Area Map with Insets
# Ian Becker
# May 2026
#
##############################

library(tigris)
library(ggplot2)
library(sf)
library(dplyr)
library(ggspatial)
library(cowplot)

options(tigris_use_cache = TRUE)

setwd("~/Desktop/project_code/UGS_occ/data")

output_dir <- "/Users/ianbecker/Desktop/project_code/UGS_occ/figures_tables"

####################
#   Load Data
####################

ugs          <- st_read("lrgv_green_spaces_detection_filtered")
census_urban <- st_read("tl_2020_us_uac20")

tx_counties <- counties(state = "TX", cb = TRUE) %>%
  st_transform(st_crs(ugs))

lrgv_counties <- c("Starr", "Hidalgo", "Cameron", "Willacy")

lrgv <- tx_counties %>%
  filter(NAME %in% lrgv_counties)

####################
#   Prep Data
####################

census_urban      <- st_transform(census_urban, st_crs(ugs))
lrgv_boundary     <- st_union(lrgv)
census_urban_lrgv <- st_intersection(census_urban, lrgv_boundary)
ugs_centroids     <- st_centroid(ugs)

####################
#   Inset 1: US Map with Texas Highlighted
####################

us_states <- states(cb = TRUE) %>%
  st_transform(st_crs(ugs)) %>%
  filter(!STUSPS %in% c("AK", "HI", "PR", "GU", "VI", "MP", "AS"))

texas <- us_states %>% filter(STUSPS == "TX")

p_us_inset <- ggplot() +
  geom_sf(data = us_states, fill = "gray90", color = "gray60", linewidth = 0.2) +
  geom_sf(data = texas, fill = "#2d6a2d", color = "gray40", linewidth = 0.3) +
  theme_void()

####################
#   Inset 2: Texas Map with LRGV Counties Highlighted
####################

p_tx_inset <- ggplot() +
  geom_sf(data = tx_counties, fill = "gray90", color = "gray60", linewidth = 0.2) +
  geom_sf(data = lrgv, fill = "#2d6a2d", color = "gray40", linewidth = 0.3) +
  theme_void()

####################
#   Save Inset Maps Separately
####################

ggsave(
  file.path(output_dir, "inset_us.png"),
  p_us_inset, width = 8, height = 5, dpi = 300, bg = "white"
)
cat("Saved: inset_us.png\n")

ggsave(
  file.path(output_dir, "inset_texas.png"),
  p_tx_inset, width = 8, height = 5, dpi = 300, bg = "white"
)
cat("Saved: inset_texas.png\n\n")

####################
#   Main Study Area Map
####################

study_area_map <- ggplot() +
  
  geom_sf(data = lrgv, fill = "gray95", color = "gray40", linewidth = 0.5) +
  
  geom_sf(data = census_urban_lrgv,
          aes(fill = "Urban Areas"), color = NA, alpha = 0.6) +
  
  geom_sf(data = ugs_centroids,
          aes(color = "Study Sites (n = 69)"), size = 1.5, alpha = 0.7) +
  
  scale_fill_manual(
    name   = NULL,
    values = c("Urban Areas" = "#D4B996")
  ) +
  
  scale_color_manual(
    name   = NULL,
    values = c("Study Sites (n = 69)" = "darkgreen"),
    guide  = guide_legend(override.aes = list(size = 3, alpha = 1))
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
    panel.grid           = element_blank(),
    axis.text            = element_text(size = 13),
    axis.title           = element_blank(),
    legend.position      = c(0.02, 0.05),
    legend.justification = c(0, 0),
    legend.background    = element_rect(fill = "white", color = NA),
    legend.text          = element_text(size = 10)
  )

print(study_area_map)

ggsave(
  file.path(output_dir, "study_area_map.png"),
  study_area_map, width = 10, height = 8, dpi = 300, bg = "white"
)

cat("Saved: study_area_map.png\n")

####################
#   Filter for any species
####################

focal_species <- "Roseate Spoonbill"

inat <- read.csv("inat_observations_with_sites.csv")

species_sites <- inat %>%
  filter(common_name == focal_species) %>%
  pull(site_id) %>%
  unique()

cat("Sites with", focal_species, "detections:", length(species_sites), "\n")

ugs_centroids_sp <- ugs_centroids %>%
  mutate(detected = site_id %in% species_sites)

species_map <- ggplot() +
  
  geom_sf(data = lrgv, fill = "gray95", color = "gray40", linewidth = 0.5) +
  
  geom_sf(data = census_urban_lrgv, fill = "#D4B996", color = NA, alpha = 0.6) +
  
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

species_filename <- gsub(" ", "_", focal_species)

ggsave(
  file.path(output_dir, paste0("study_area_map_", species_filename, ".png")),
  species_map, width = 10, height = 8, dpi = 300, bg = "white"
)

cat("Saved: study_area_map_", species_filename, ".png\n")