##############################
#
# Exemplar Site Map with Land Cover and Grid
# Ian Becker
# January 2026
#
##############################

library(sf)
library(terra)
library(ggplot2)
library(tidyterra)
library(dplyr)

# Set paths
input_dir <- "/Users/ianbecker/Desktop/project_code/UGS_occ/data"
output_dir <- "/Users/ianbecker/Desktop/project_code/UGS_occ/figures_tables"

####################
#   Settings
####################

exemplar_site_id <- 191  # Change this to try different sites
cell_size <- 200  # meters

####################
#   Load Data
####################

cat("Loading data...\n")

# Sites
sites <- st_read(file.path(input_dir, "lrgv_green_spaces_detection_filtered"), quiet = TRUE)
sites_utm <- st_transform(sites, crs = 32614)

# Get exemplar site
site <- sites_utm %>% filter(site_id == exemplar_site_id)

if (nrow(site) == 0) {
  stop(paste("Site", exemplar_site_id, "not found!"))
}

site_area <- round(as.numeric(st_area(site)) / 10000, 2)
cat("Site", exemplar_site_id, ":", site_area, "ha\n")

# Land cover
landcover <- rast(file.path(input_dir, "lrgv_dynamic_world_cover.tif"))

# Ensure same CRS
if (crs(landcover, describe = TRUE)$code != "32614") {
  landcover <- project(landcover, "EPSG:32614", method = "near")
}

####################
#   Crop Land Cover to Site (with buffer for context)
####################

cat("Cropping land cover...\n")

# Add buffer around site for context
site_buffered <- st_buffer(site, dist = 100)

# Crop raster
lc_cropped <- crop(landcover, vect(site_buffered))
lc_masked <- mask(lc_cropped, vect(site_buffered))

####################
#   Create Hex Grid
####################

cat("Creating hexagonal grid...\n")

hex_grid <- st_make_grid(site, cellsize = cell_size, square = FALSE)
hex_sf <- st_sf(geometry = hex_grid)
hex_filtered <- hex_sf[lengths(st_intersects(hex_sf, site)) > 0, ]
hex_filtered$cell_id <- 1:nrow(hex_filtered)

cat("Created", nrow(hex_filtered), "hexagonal cells\n")

####################
#   Create Land Cover Categories for Plotting
####################

# Dynamic World classes
lc_classes <- data.frame(
  value = 0:8,
  class = c("Water", "Trees", "Grass", "Flooded Veg", "Crops", 
            "Shrub", "Built", "Bare", "Snow/Ice")
)

# Create factor raster for plotting
lc_factor <- as.factor(lc_masked)

# Define colors (Dynamic World palette-ish)
lc_colors <- c(
  "Water" = "#419BDF",
  "Trees" = "#397D49", 
  "Grass" = "#88B053",
  "Flooded Veg" = "#7A87C6",
  "Crops" = "#E49635",
  "Shrub" = "#DFC35A",
  "Built" = "#C4281B",
  "Bare" = "#A59B8F",
  "Snow/Ice" = "#FFFFFF"
)

####################
#   Load iNat observations for this site (optional overlay)
####################

inat <- read.csv(file.path(input_dir, "inat_observations_with_sites.csv"))

# Filter to this site and study species
waterbirds <- c("Black-bellied Whistling-Duck", "Green Heron", 
                "Neotropic Cormorant", "Great Blue Heron", "Great Egret")

site_obs <- inat %>%
  filter(site_id == exemplar_site_id,
         common_name %in% waterbirds)

cat("Observations at site:", nrow(site_obs), "\n")

# Convert to spatial
if (nrow(site_obs) > 0) {
  obs_sf <- st_as_sf(site_obs, coords = c("longitude", "latitude"), crs = 4326)
  obs_sf <- st_transform(obs_sf, crs = 32614)
}

####################
#   Create Map
####################

cat("Creating map...\n")

# Get bounding box for plot limits
bbox <- st_bbox(site_buffered)

# Convert raster to factor for discrete colors
lc_factor <- as.factor(lc_masked)

# Define colors
lc_colors <- c(
  "0" = "#419BDF",
  "1" = "#397D49", 
  "2" = "#88B053",
  "3" = "#7A87C6",
  "4" = "#E49635",
  "5" = "#DFC35A",
  "6" = "#C4281B",
  "7" = "#A59B8F",
  "8" = "#FFFFFF"
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

# Base plot with land cover
p <- ggplot() +
  # Land cover raster
  geom_spatraster(data = lc_factor) +
  scale_fill_manual(
    values = lc_colors,
    labels = lc_labels,
    name = "Land Cover",
    na.value = "transparent",
    na.translate = FALSE
  ) +
  
  # Site boundary
  geom_sf(data = site, fill = NA, color = "black", linewidth = 1.2) +
  
  # Hex grid
  geom_sf(data = hex_filtered, fill = NA, color = "white", linewidth = 0.6) +
  
  # Coordinate limits
  coord_sf(xlim = c(bbox["xmin"], bbox["xmax"]),
           ylim = c(bbox["ymin"], bbox["ymax"]),
           expand = FALSE) +
  
  # Labels
  labs(
    title = paste0("Site ", exemplar_site_id)
  ) +
  
  # Theme - no axes
  theme_void() +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5, color = "gray40"),
    legend.position = "right",
    plot.margin = margin(10, 10, 10, 10)
  )

print(p)
# Save
ggsave(file.path(output_dir, paste0("site_", exemplar_site_id, "_landcover_grid.png")),
       p, width = 10, height = 8, dpi = 300)

cat("Saved map!\n")

####################
#   Version with observations overlaid
####################

if (nrow(site_obs) > 0) {
  
  p_obs <- p +
    geom_sf(data = obs_sf, aes(color = common_name), size = 1.5, alpha = 0.7) +
    scale_color_brewer(palette = "Set1", name = "Species") +
    labs(
      title = paste0("Site ", exemplar_site_id, " - Waterbird Observations"),
      subtitle = paste0("Area: ", site_area, " ha | ", nrow(site_obs), " observations")
    )
  
  print(p_obs)
  
  ggsave(file.path(output_dir, paste0("site_", exemplar_site_id, "_landcover_observations.png")),
         p_obs, width = 12, height = 8, dpi = 300)
  
  cat("Saved map with observations!\n")
}

####################
#   Summary stats for this site
####################

cat("\n=== SITE SUMMARY ===\n")

# Land cover composition
lc_vals <- terra::extract(lc_masked, vect(site), df = TRUE)
lc_summary <- lc_vals %>%
  filter(!is.na(label)) %>%
  count(label) %>%
  mutate(pct = round(n / sum(n) * 100, 1))

cat("\nLand cover composition:\n")
for (i in 1:nrow(lc_summary)) {
  class_name <- lc_classes$class[lc_classes$value == lc_summary$label[i]]
  cat("  ", class_name, ":", lc_summary$pct[i], "%\n")
}

# Species at this site
cat("\nSpecies observations:\n")
sp_summary <- site_obs %>%
  count(common_name) %>%
  arrange(desc(n))
print(sp_summary)

cat("\n=== DONE ===\n")