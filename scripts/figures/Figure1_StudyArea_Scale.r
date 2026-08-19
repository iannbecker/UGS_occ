##############################
#
# Figure 1: Study Area and Scale Framework
# Ian Becker
# May 2026
#
##############################

# This script produces all the panels for Figure 1 in the manuscript.
# This includes:
#   - Study area map (panel A)
#   - US and Texas inset maps (saved separately)
#   - Clean site maps for tilting (saved as PNGs)
#   - Tilted versions of both site maps (panels B and C)

library(tigris)
library(ggplot2)
library(sf)
library(terra)
library(tidyterra)
library(dplyr)
library(ggspatial)
library(magick)

options(tigris_use_cache = TRUE)

# ============================================================================
# 1. SCRIPT PREP
# ============================================================================

exemplar_site_id <- 48     # site ID for tilted maps (see shapefile for more info)
cell_size <- 100    # hex grid cell size in meters

input_dir  <- "/Users/ianbecker/Desktop/project_code/UGS_occ/data"
output_dir <- "/Users/ianbecker/Desktop/project_code/UGS_occ/figures_tables"

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

setwd(input_dir)

# ============================================================================
# 2. LOAD DATA 
# ============================================================================

# Load in urban green space and urban areas shapefile

ugs          <- st_read("lrgv_green_spaces_ebird", quiet = TRUE)
census_urban <- st_read("tl_2020_us_uac20", quiet = TRUE)

# Load in county shapefile

tx_counties <- counties(state = "TX", cb = TRUE) %>%
  st_transform(st_crs(ugs))

# Filter to LRGV

lrgv_counties <- c("Starr", "Hidalgo", "Cameron", "Willacy")
lrgv <- tx_counties %>%
  filter(NAME %in% lrgv_counties)

# Align and clip to LRGV

census_urban      <- st_transform(census_urban, st_crs(ugs))
lrgv_boundary     <- st_union(lrgv)
census_urban_lrgv <- st_intersection(census_urban, lrgv_boundary)
ugs_centroids     <- st_centroid(ugs)

# Exemplar site (see script prep)

sites_utm <- st_transform(ugs, crs = 32614)
site      <- sites_utm %>% filter(site_id == exemplar_site_id)
if (nrow(site) == 0) stop(paste("Site", exemplar_site_id, "not found!"))

# Land cover

landcover <- rast("lrgv_dynamic_world_cover.tif")
if (crs(landcover, describe = TRUE)$code != "32614") {
  landcover <- project(landcover, "EPSG:32614", method = "near")
}

cat("Data loaded\n\n")

# ============================================================================
# 3. PANEL A: STUDY AREA MAP
# ============================================================================

# Load in US shapefile

us_states <- states(cb = TRUE) %>%
  st_transform(st_crs(ugs)) %>%
  filter(!STUSPS %in% c("AK", "HI", "PR", "GU", "VI", "MP", "AS"))

# Filter to Texas

texas <- us_states %>% filter(STUSPS == "TX")

# Create study area map

study_area_map <- ggplot() +
  
  geom_sf(data = lrgv, fill = "gray95", color = "gray40", linewidth = 0.5) +
  
  geom_sf(data = census_urban_lrgv,
          aes(fill = "Urban Areas"), color = NA, alpha = 0.6) +
  
  geom_sf(data = ugs_centroids,
          aes(color = "Study Sites (n = 59)"), size = 1.5, alpha = 0.7) +
  
  scale_fill_manual(
    name   = NULL,
    values = c("Urban Areas" = "#D4B996")
  ) +
  
  scale_color_manual(
    name   = NULL,
    values = c("Study Sites (n = 59)" = "darkgreen"),
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

# Save

ggsave(
  file.path(output_dir, "fig1A_study_area_map.png"),
  study_area_map, width = 10, height = 8, dpi = 300, bg = "white"
)

cat("Saved: fig1A_study_area_map.png\n")

# ============================================================================
# 4. INSET MAPS FOR PANEL A
# ============================================================================

# US inset map

p_us_inset <- ggplot() +
  geom_sf(data = us_states, fill = "gray90", color = "gray60", linewidth = 0.2) +
  geom_sf(data = texas, fill = "#2d6a2d", color = "gray40", linewidth = 0.3) +
  theme_void()

# Texas inset map

p_tx_inset <- ggplot() +
  geom_sf(data = tx_counties, fill = "gray90", color = "gray60", linewidth = 0.2) +
  geom_sf(data = lrgv, fill = "#2d6a2d", color = "gray40", linewidth = 0.3) +
  theme_void()

ggsave(file.path(output_dir, "fig1_inset_us.png"),
       p_us_inset, width = 8, height = 5, dpi = 300, bg = "white")
cat("Saved: fig1_inset_us.png\n")

ggsave(file.path(output_dir, "fig1_inset_texas.png"),
       p_tx_inset, width = 8, height = 5, dpi = 300, bg = "white")
cat("Saved: fig1_inset_texas.png\n\n")

# ============================================================================
# 5. PANEL B: EXEMPLAR SITE MAP
# ============================================================================

# Prep Land Cover data

# Crop and mask to land cover raster

site_buffered <- st_buffer(site, dist = 500)
lc_cropped    <- crop(landcover, vect(site_buffered))
lc_masked     <- mask(lc_cropped, vect(st_buffer(site, dist = 100)))
lc_factor     <- as.factor(lc_masked)

# Create grid and filter to hexagons that intersect the site

hex_grid     <- st_make_grid(site, cellsize = cell_size, square = FALSE)
hex_sf       <- st_sf(geometry = hex_grid)
hex_filtered <- hex_sf[lengths(st_intersects(hex_sf, site)) > 0, ]
hex_filtered$cell_id <- 1:nrow(hex_filtered)

cat("Created", nrow(hex_filtered), "hexagonal cells\n\n")

# Classify land cover class by color

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
  "0" = "Water",    "1" = "Trees",
  "2" = "Grass",    "3" = "Flooded Veg",
  "4" = "Crops",    "5" = "Shrub",
  "6" = "Built",    "7" = "Bare",
  "8" = "Snow/Ice"
)

# Set up bounding box by site

bbox <- st_bbox(site_buffered)


# Plot clean land cover site map (no grid)

p_landcover <- ggplot() +
  geom_spatraster(data = lc_factor) +
  scale_fill_manual(
    values       = lc_colors,
    labels       = lc_labels,
    na.value     = "transparent",
    na.translate = FALSE
  ) +
  coord_sf(
    xlim   = c(bbox["xmin"], bbox["xmax"]),
    ylim   = c(bbox["ymin"], bbox["ymax"]),
    expand = FALSE
  ) +
  theme_void() +
  theme(legend.position = "none")

# Save

path_landcover <- file.path(output_dir,
                            paste0("fig1B_site_", exemplar_site_id, "_landcover.png"))

ggsave(path_landcover, p_landcover,
       width = 10, height = 8, dpi = 300, bg = "white")
cat("Saved:", path_landcover, "\n")

# ============================================================================
# 6. PANEL C: EXEMPLAR SITE MAP WITH GRID
# ============================================================================

# Plot gridded site map

p_grid <- ggplot() +
  geom_spatraster(data = lc_factor) +
  scale_fill_manual(
    values       = lc_colors,
    labels       = lc_labels,
    na.value     = "transparent",
    na.translate = FALSE
  ) +
  geom_sf(data = hex_filtered, fill = NA, color = "white",
          linewidth = 0.35, alpha = 0.6) +
  coord_sf(
    xlim   = c(bbox["xmin"], bbox["xmax"]),
    ylim   = c(bbox["ymin"], bbox["ymax"]),
    expand = FALSE
  ) +
  theme_void() +
  theme(legend.position = "none")

path_grid <- file.path(output_dir,
                       paste0("fig1C_site_", exemplar_site_id, "_landcover_grid.png"))

ggsave(path_grid, p_grid,
       width = 10, height = 8, dpi = 300, bg = "white")
cat("Saved:", path_grid, "\n\n")

