##############################
#
# Exemplar Site Map with Detections
# Ian Becker
# May 2026
#
##############################

library(sf)
library(terra)
library(ggplot2)
library(tidyterra)
library(dplyr)

####################
#   USER INPUTS — change these
####################

exemplar_site_id <- 25        # site ID
cell_size        <- 100        # meters
target_species   <- "Golden-fronted Woodpecker"  # common name as in iNat data

####################
#   Paths
####################

input_dir  <- "/Users/ianbecker/Desktop/project_code/UGS_occ/data"
output_dir <- "/Users/ianbecker/Desktop/project_code/UGS_occ/figures_tables/site_maps"

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

####################
#   Load Data
####################

cat("Loading data...\n")

sites     <- st_read(file.path(input_dir, "lrgv_green_spaces_detection_filtered"), quiet = TRUE)
sites_utm <- st_transform(sites, crs = 32614)
site      <- sites_utm %>% filter(site_id == exemplar_site_id)

if (nrow(site) == 0) stop(paste("Site", exemplar_site_id, "not found!"))

site_area <- round(as.numeric(st_area(site)) / 10000, 2)
cat("Site", exemplar_site_id, ":", site_area, "ha\n")

# iNat observations
inat <- read.csv(file.path(input_dir, "inat_observations_with_sites.csv"))

# Land cover
landcover <- rast(file.path(input_dir, "lrgv_dynamic_world_cover.tif"))
if (crs(landcover, describe = TRUE)$code != "32614") {
  landcover <- project(landcover, "EPSG:32614", method = "near")
}

cat("Data loaded\n\n")

####################
#   Crop Land Cover
####################

site_buffered <- st_buffer(site, dist = 100)
lc_cropped    <- crop(landcover, vect(site_buffered))
lc_masked     <- mask(lc_cropped, vect(site_buffered))
lc_factor     <- as.factor(lc_masked)

####################
#   Create Hex Grid
####################

hex_grid     <- st_make_grid(site, cellsize = cell_size, square = FALSE)
hex_sf       <- st_sf(geometry = hex_grid)
hex_filtered <- hex_sf[lengths(st_intersects(hex_sf, site)) > 0, ]
hex_filtered$cell_id <- 1:nrow(hex_filtered)

cat("Created", nrow(hex_filtered), "hexagonal cells\n")

####################
#   Filter Observations
####################

# All observations at this site (any species) — for effort context
all_site_obs <- inat %>% filter(site_id == exemplar_site_id)

# Target species observations at this site
species_obs <- inat %>%
  filter(site_id == exemplar_site_id,
         common_name == target_species)

cat("All observations at site:", nrow(all_site_obs), "\n")
cat(target_species, "observations at site:", nrow(species_obs), "\n\n")

if (nrow(species_obs) == 0) {
  warning(paste("No observations of", target_species, "at site", exemplar_site_id))
}

# Convert species observations to spatial
if (nrow(species_obs) > 0) {
  species_sf <- st_as_sf(species_obs,
                         coords = c("longitude", "latitude"),
                         crs = 4326)
  species_sf <- st_transform(species_sf, crs = 32614)
}

####################
#   Land Cover Colors and Labels
####################

lc_colors <- c(
  "0" = "#419BDF",  # Water
  "1" = "#397D49",  # Trees
  "2" = "#88B053",  # Grass
  "3" = "#7A87C6",  # Flooded Veg
  "4" = "#E49635",  # Crops
  "5" = "#DFC35A",  # Shrub
  "6" = "#C4281B",  # Built
  "7" = "#A59B8F",  # Bare
  "8" = "#FFFFFF"   # Snow/Ice
)

lc_labels <- c(
  "0" = "Water",
  "1" = "Trees",
  "2" = "Grass",
  "3" = "Flooded Veg",
  "4" = "Crops",
  "5" = "Shrub",
  "6" = "Built",
  "7" = "Bare",
  "8" = "Snow/Ice"
)

####################
#   Bounding Box
####################

bbox <- st_bbox(site_buffered)

####################
#   Plot 1: Land Cover + Grid Only
####################

cat("Creating base map...\n")

p_base <- ggplot() +
  geom_spatraster(data = lc_factor) +
  scale_fill_manual(
    values     = lc_colors,
    labels     = lc_labels,
    name       = "Land Cover",
    na.value   = "transparent",
    na.translate = FALSE
  ) +
  geom_sf(data = hex_filtered, fill = NA, color = "white",  linewidth = 0.35, alpha = 0.6) +
  coord_sf(xlim = c(bbox["xmin"], bbox["xmax"]),
           ylim = c(bbox["ymin"], bbox["ymax"]),
           expand = FALSE) +
  labs(title = paste0("Site ", exemplar_site_id)) +
  theme_void() +
  theme(
    plot.title      = element_text(size = 14, face = "bold", hjust = 0.5),
    legend.position = "right",
    plot.margin     = margin(10, 10, 10, 10)
  )

ggsave(
  file.path(output_dir, paste0("site_", exemplar_site_id, "_landcover_grid.png")),
  p_base, width = 10, height = 8, dpi = 300, bg = "white"
)
cat("Saved base map\n")

####################
#   Plot 2: Land Cover + Grid + Species Detections
####################

if (nrow(species_obs) > 0) {
  
  cat("Creating detection map...\n")
  
  species_label <- gsub(" ", "_", target_species)
  
  # Detection label for legend including n
  det_label <- paste0("Detection (n = ", nrow(species_obs), ")")
  
  p_detections <- p_base +
    
    # Detection points — mapped to shape aesthetic so they appear in legend
    geom_sf(
      data  = species_sf,
      aes(shape = det_label),
      fill  = "yellow",
      color = "black",
      size  = 2.5,
      stroke = 0.6,
      alpha = 0.85
    ) +
    
    scale_shape_manual(
      name   = NULL,
      values = setNames(21, det_label)
    ) +
    
    guides(
      fill  = guide_legend(order = 1, title = "Land Cover"),
      shape = guide_legend(order = 2,
                           override.aes = list(fill = "yellow", color = "black",
                                               size = 3, stroke = 0.6))
    ) +
    
    labs(title = paste0("Site ", exemplar_site_id, " \u2014 ", target_species)) +
    
    theme(plot.subtitle = element_blank())
  
  ggsave(
    file.path(output_dir, 
              paste0("site_", exemplar_site_id, "_", species_label, "_detections.png")),
    p_detections, width = 10, height = 8, dpi = 300, bg = "white"
  )
  cat("Saved detection map\n\n")
  
} else {
  cat("No detections to plot for", target_species, "at site", exemplar_site_id, "\n\n")
}

####################
#   Console Summary
####################

cat("=== SITE SUMMARY ===\n")
cat("Site:", exemplar_site_id, "\n")
cat("Area:", site_area, "ha\n")
cat("Cells:", nrow(hex_filtered), "\n")
cat("Total observations (all species):", nrow(all_site_obs), "\n")
cat(target_species, "observations:", nrow(species_obs), "\n")

if (nrow(species_obs) > 0) {
  cat("\nDetections by year:\n")
  species_obs$year <- as.integer(substr(species_obs$observed_on, 1, 4))
  year_summary <- species_obs %>% count(year) %>% arrange(year)
  print(year_summary)
}

cat("\n=== DONE ===\n")
cat("Outputs saved to:", output_dir, "\n")

