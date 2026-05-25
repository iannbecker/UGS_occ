##############################
#
# Tilted Site Map Generator
# Ian Becker
# May 2026
#
##############################

# Step 1: Generate clean site maps (no legend, no title, no axes)
#         - Map A: land cover only
#         - Map B: land cover + hex grid
# Step 2: Apply perspective tilt to both maps using magick

library(sf)
library(terra)
library(ggplot2)
library(tidyterra)
library(dplyr)
library(magick)

####################
#   USER INPUTS
####################

exemplar_site_id <- 48     # site ID
cell_size        <- 100    # meters

input_dir  <- "/Users/ianbecker/Desktop/project_code/UGS_occ/data"
output_dir <- "/Users/ianbecker/Desktop/project_code/UGS_occ/figures_tables/tilted_maps"

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

####################
#   Load Data
####################

cat("Loading data...\n")

sites     <- st_read(file.path(input_dir, "lrgv_green_spaces_detection_filtered"), quiet = TRUE)
sites_utm <- st_transform(sites, crs = 32614)
site      <- sites_utm %>% filter(site_id == exemplar_site_id)

if (nrow(site) == 0) stop(paste("Site", exemplar_site_id, "not found!"))

landcover <- rast(file.path(input_dir, "lrgv_dynamic_world_cover.tif"))
if (crs(landcover, describe = TRUE)$code != "32614") {
  landcover <- project(landcover, "EPSG:32614", method = "near")
}

cat("Data loaded\n\n")

####################
#   Prep Land Cover
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

cat("Created", nrow(hex_filtered), "hexagonal cells\n\n")

####################
#   Land Cover Colors
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
  "0" = "Water",    "1" = "Trees",
  "2" = "Grass",    "3" = "Flooded Veg",
  "4" = "Crops",    "5" = "Shrub",
  "6" = "Built",    "7" = "Bare",
  "8" = "Snow/Ice"
)

bbox <- st_bbox(st_buffer(site, dist = 500))

####################
#   Step 1A: Land Cover Only — No legend, no title, no axes
####################

cat("Generating clean land cover map...\n")

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

path_landcover <- file.path(output_dir,
                            paste0("site_", exemplar_site_id, "_landcover_clean.png"))

ggsave(path_landcover, p_landcover,
       width = 10, height = 8, dpi = 300, bg = "white")
cat("Saved:", path_landcover, "\n")

####################
#   Step 1B: Land Cover + Hex Grid — No legend, no title, no axes
####################

cat("Generating clean land cover + grid map...\n")

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
                       paste0("site_", exemplar_site_id, "_landcover_grid_clean.png"))

ggsave(path_grid, p_grid,
       width = 10, height = 8, dpi = 300, bg = "white")
cat("Saved:", path_grid, "\n\n")

####################
#   Step 2: Apply Perspective Tilt
####################

cat("Applying perspective tilt...\n")

tilt_perspective <- function(img) {
  
  w <- image_info(img)$width
  h <- image_info(img)$height
  
  # Crop a small amount off the top before tilting to remove stray pixels
  img <- image_crop(img, paste0(w, "x", round(h * 0.90), "+0+", round(h * 0.05)))
  
  w <- image_info(img)$width
  h <- image_info(img)$height
  
  control_points <- c(
    0,   0,        round(w * 0.18),  round(h * 0.20),
    w,   0,        round(w * 0.91),  round(h * 0.05),
    0,   h,        round(w * 0.06),  round(h * 0.57),
    w,   h,        round(w * 0.80),  round(h * 0.47)
  )
  
  image_distort(img, "Perspective", control_points, bestfit = TRUE)
}
# Load and tilt
img_landcover <- image_read(path_landcover)
img_grid      <- image_read(path_grid)

img_landcover_tilted <- tilt_perspective(img_landcover) %>% image_trim()
img_grid_tilted      <- tilt_perspective(img_grid)      %>% image_trim()

# Save tilted versions
path_landcover_tilted <- file.path(output_dir,
                                   paste0("site_", exemplar_site_id, "_landcover_tilted.png"))
path_grid_tilted      <- file.path(output_dir,
                                   paste0("site_", exemplar_site_id, "_grid_tilted.png"))

image_write(img_landcover_tilted, path_landcover_tilted)
cat("Saved:", path_landcover_tilted, "\n")

image_write(img_grid_tilted, path_grid_tilted)
cat("Saved:", path_grid_tilted, "\n\n")

cat("=== DONE ===\n")
cat("Tilted maps saved to:", output_dir, "\n")