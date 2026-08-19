##############################
#
# Figure 5: Roseate Spoonbill Scale Example
# Ian Becker
# May 2026
#
##############################

# This script generates all panels for Figure 5.
# This includes:
#   Panel A — Combined landscape/site level water cover occupancy curves
#   Panel B — LRGV map with (SPECIES NAME) detection sites highlighted
#   Panel C — Site ## land cover + hex grid + species detections

library(sf)
library(terra)
library(tidyterra)
library(ggplot2)
library(dplyr)
library(tigris)
library(ggspatial)

options(tigris_use_cache = TRUE)

# ============================================================================
# 1. SSCRIPT PREP
# ============================================================================

focal_species    <- "Green Jay"
focal_site_id    <- 40
focal_cov        <- "trees"      # covariate name in within-site models
focal_cov_land   <- "trees_pct"  # covariate name in landscape models
cell_size        <- 100          # hex grid cell size in meters
n_points         <- 200
ci               <- 0.95

input_dir  <- "/Users/ianbecker/Desktop/project_code/UGS_occ/data"
output_dir <- "/Users/ianbecker/Desktop/project_code/UGS_occ/figures_tables/figureS2_greenjay"

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

setwd(input_dir)

# Colors for land cover and spatial scales

# Landscape vs. Site level

land_color <- "#2d6a4f"   # dark green — landscape level
site_color <- "#9b2335"   # dark red — site level

# Land cover colors

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

# Land cover labels

lc_labels <- c(
  "0" = "Water",    "1" = "Trees",
  "2" = "Grass",    "3" = "Flooded Veg",
  "4" = "Crops",    "5" = "Shrub",
  "6" = "Built",    "7" = "Bare",
  "8" = "Snow/Ice"
)

# ============================================================================
# 2. LOAD DATA
# ============================================================================

# iNat data and site covariates

inat <- read.csv("gbif_in_sites.csv")
ebird <- read.csv("ebird_filtered_lrgv_with_sites.csv")
site_covs <- read.csv("site_covariates_ebird.csv")

# Urban green space shapefile + check

ugs       <- st_read("lrgv_green_spaces_ebird", quiet = TRUE)
sites_utm <- st_transform(ugs, crs = 32614)
site      <- sites_utm %>% filter(site_id == focal_site_id)

if (nrow(site) == 0) stop(paste("Site", focal_site_id, "not found!"))

# Land cover raw data

landcover <- rast("lrgv_dynamic_world_cover.tif")
if (crs(landcover, describe = TRUE)$code != "32614") {
  landcover <- project(landcover, "EPSG:32614", method = "near")
}

# Landscape models

all_results <- readRDS("model_results_ebird_checklist/all_results_ebird_checklist_2026-08-11.rds")

# Within site models

site_file <- file.path("within_site_models_gbif",
                       paste0(gsub(" ", "_", focal_species),
                              "_site", focal_site_id, ".rds"))
if (!file.exists(site_file)) stop("Within-site model file not found")
site_result <- readRDS(site_file)

cat("Data loaded\n\n")

# ============================================================================
# 3. PANEL A: OCCUPANCY SCALE COMPARISON
# ============================================================================

# Back transform covariates to original scale

cov_mean <- mean(site_covs$trees_pct, na.rm = TRUE)
cov_sd   <- sd(site_covs$trees_pct,   na.rm = TRUE)
cov_min  <- min(site_covs$trees_pct,  na.rm = TRUE)
cov_max  <- max(site_covs$trees_pct,  na.rm = TRUE)

# Prediction on original scale

x_orig   <- seq(cov_min, cov_max, length.out = n_points)

# Standardize for model predictions

x_scaled <- (x_orig - cov_mean) / cov_sd

### Landscape level curve 

# Find focal species

sp_result <- NULL
for (nm in names(all_results)) {
  if (all_results[[nm]]$species == focal_species && all_results[[nm]]$success) {
    sp_result <- all_results[[nm]]
    break
  }
}
if (is.null(sp_result)) stop("Species not found in landscape models")

# Extract posterior samples from landscape models

model_land   <- sp_result$model
beta_land    <- as.matrix(model_land$beta.samples)

# Clean names and find focal covariate

params_clean <- gsub("[0-9]+$", "", colnames(beta_land))
int_idx      <- which(params_clean == "(Intercept)")
cov_idx      <- which(params_clean == focal_cov_land)

# Effect of covariate holding everything at standardized mean

pred_land <- matrix(NA, nrow = nrow(beta_land), ncol = n_points)
for (i in 1:nrow(beta_land)) {
  pred_land[i, ] <- plogis(beta_land[i, int_idx] + beta_land[i, cov_idx] * x_scaled)
}

# Summarize posterior predictions 

land_df <- data.frame(
  x     = x_orig,
  mean  = apply(pred_land, 2, mean),
  lower = apply(pred_land, 2, quantile, (1 - ci) / 2),
  upper = apply(pred_land, 2, quantile, 1 - (1 - ci) / 2),
  scale = "Landscape Level"
)

### Site level curve 

model_site  <- site_result$model
beta_site   <- as.matrix(model_site$beta.samples)

# Get scaling parameters from original models

if (!is.null(site_result$covariates_unscaled) &&
    focal_cov %in% names(site_result$covariates_unscaled)) {
  site_cov      <- site_result$covariates_unscaled[[focal_cov]]
  site_cov_mean <- mean(site_cov, na.rm = TRUE)
  site_cov_sd   <- sd(site_cov,   na.rm = TRUE)
} else {
  site_cov_mean <- cov_mean
  site_cov_sd   <- cov_sd
}

# Clean (see above)

params_site_clean <- gsub("[0-9]+$", "", colnames(beta_site))
int_site_idx      <- which(params_site_clean == "(Intercept)")
cov_site_idx      <- which(params_site_clean == focal_cov)

# Standardize same as above to share x-axis

x_site_scaled <- (x_orig - site_cov_mean) / site_cov_sd

# Predict within-site occupancy across range (same as above)

pred_site <- matrix(NA, nrow = nrow(beta_site), ncol = n_points)
for (i in 1:nrow(beta_site)) {
  pred_site[i, ] <- plogis(beta_site[i, int_site_idx] +
                             beta_site[i, cov_site_idx] * x_site_scaled)
}

# Summarize within-site posterior predictions

site_df <- data.frame(
  x     = x_orig,
  mean  = apply(pred_site, 2, mean),
  lower = apply(pred_site, 2, quantile, (1 - ci) / 2),
  upper = apply(pred_site, 2, quantile, 1 - (1 - ci) / 2),
  scale = "Site Level"
)

### Raw data points

# Filter to subsampled checklists used in the actual analysis

effort_lookup <- readRDS("detection_matrices_ebird_checklist/effort_lookup.rds")
sampled_ids   <- unique(effort_lookup$checklist_id)

# Panel B detection sites — from subsampled checklists only

sp_sites <- ebird %>%
  filter(checklist_id %in% sampled_ids,
         common_name == focal_species,
         species_observed == TRUE) %>%
  pull(site_id) %>%
  unique()

# Get landscape-level detections 

land_pts <- site_covs %>%
  mutate(detected = site_id %in% sp_sites) %>%
  transmute(x = trees_pct, y = as.numeric(detected), scale = "Landscape Level")

# Get site-level detections

det_matrix <- site_result$detection_matrix
naive_occ  <- apply(det_matrix, 1, function(x) as.numeric(any(x == 1, na.rm = TRUE)))

site_pts <- if (!is.null(site_result$covariates_unscaled) &&
                focal_cov %in% names(site_result$covariates_unscaled)) {
  data.frame(x     = site_result$covariates_unscaled[[focal_cov]],
             y     = naive_occ,
             scale = "Site Level")
} else NULL

# Combine data for plotting

all_pts    <- rbind(land_pts, site_pts)
all_curves <- rbind(land_df, site_df)

all_curves$scale <- factor(all_curves$scale, levels = c("Landscape Level", "Site Level"))
all_pts$scale    <- factor(all_pts$scale,    levels = c("Landscape Level", "Site Level"))

scale_colors <- c("Landscape Level" = land_color, "Site Level" = site_color)

# Plot Panel A

p_curve <- ggplot() +
  geom_ribbon(data = all_curves,
              aes(x = x, ymin = lower, ymax = upper, fill = scale),
              alpha = 0.25) +
  geom_line(data = all_curves,
            aes(x = x, y = mean, color = scale),
            linewidth = 1.2) +
  geom_point(data = all_pts,
             aes(x = x, y = y, color = scale),
             size = 2, alpha = 0.5) +
  scale_color_manual(values = scale_colors, name = NULL) +
  scale_fill_manual( values = scale_colors, name = NULL) +
  scale_x_continuous(limits = c(cov_min, cov_max), n.breaks = 5) +
  scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  labs(x = "cov Cover (%)", y = "Occupancy Probability") +
  theme_classic() +
  theme(
    axis.text            = element_text(size = 15),
    axis.title           = element_blank(),
    legend.position = "bottom",
    legend.text          = element_text(size = 17),
    legend.background    = element_rect(fill = "white", color = NA),
    panel.grid.major.y   = element_line(color = "gray90", linewidth = 0.3)
  )

# Save 

ggsave(file.path(output_dir, "figS2_scale_comparison_curve.png"),
       p_curve, width = 6, height = 5, dpi = 300, bg = "white")
cat("Saved: figS2_scale_comparison_curve.png\n\n")

# ============================================================================
# 4. PANEL B: LANDSCAPE-LEVEL DETECTIONS 
# ============================================================================

# Pull LRGV counties

tx_counties <- counties(state = "TX", cb = TRUE) %>%
  st_transform(st_crs(ugs))

lrgv <- tx_counties %>%
  filter(NAME %in% c("Starr", "Hidalgo", "Cameron", "Willacy"))

# Load in urban areas shapefile

census_urban <- st_read("tl_2020_us_uac20", quiet = TRUE) %>%
  st_transform(st_crs(ugs))

# Align to LRGV boundary

lrgv_boundary     <- st_union(lrgv)
census_urban_lrgv <- st_intersection(census_urban, lrgv_boundary)
ugs_centroids     <- st_centroid(ugs)

# UGS with focal species detections

ugs_centroids_sp <- ugs_centroids %>%
  mutate(detected = site_id %in% sp_sites)

# Plot Panel B

p_map <- ggplot() +
  geom_sf(data = lrgv, fill = "gray95", color = "gray40", linewidth = 0.5) +
  geom_sf(data = census_urban_lrgv, fill = "#D4B996", color = NA, alpha = 0.6) +
  geom_sf(data  = ugs_centroids_sp %>% filter(detected),
          shape = 21, fill = "yellow", color = "black",
          size = 5, stroke = 0.6, alpha = 0.85) +
  geom_sf_text(data = lrgv, aes(label = NAME), size = 3.5, color = "gray30") +
  annotation_scale(location = "bl", width_hint = 0.25,
                   style = "ticks", text_cex = 0.8) +
  annotation_north_arrow(location = "tr", which_north = "true",
                         height = unit(1, "cm"), width = unit(1, "cm"),
                         style = north_arrow_fancy_orienteering()) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text  = element_text(size = 12),
    axis.title = element_blank()
  )

# Save

ggsave(file.path(output_dir, "figS2_species_detection_map.png"),
       p_map, width = 10, height = 8, dpi = 300, bg = "white")
cat("Saved: fig5B_species_detection_map.png\n\n")

# ============================================================================
# 5. PANEL C: SITE-LEVEL DETECTIONS
# ============================================================================

# Round site area 

site_area <- round(as.numeric(st_area(site)) / 10000, 2)

# Create site map buffer

site_buffered <- st_buffer(site, dist = 100)

# Align land cover data with site

lc_cropped    <- crop(landcover, vect(site_buffered))
lc_masked     <- mask(lc_cropped, vect(site_buffered))
lc_factor     <- as.factor(lc_masked)

# Construct hexagonal grid

hex_grid     <- st_make_grid(site, cellsize = cell_size, square = FALSE)
hex_sf       <- st_sf(geometry = hex_grid)
hex_filtered <- hex_sf[lengths(st_intersects(hex_sf, site)) > 0, ]
hex_filtered$cell_id <- 1:nrow(hex_filtered)

# Get species detections at focal site

species_obs <- inat %>%
  filter(site_id == focal_site_id, common_name == focal_species)

cat(focal_species, "observations at site:", nrow(species_obs), "\n\n")

# Convert species observations to sf points

species_sf <- st_as_sf(species_obs,
                       coords = c("decimalLongitude", "decimalLatitude"), crs = 4326) %>%
  st_transform(crs = 32614)

# Set up bounding box

bbox      <- st_bbox(site_buffered)
det_label <- paste0("Detection (n = ", nrow(species_obs), ")")

# Plot Panel C

p_base <- ggplot() +
  geom_spatraster(data = lc_factor) +
  scale_fill_manual(values = lc_colors, labels = lc_labels,
                    name = "Land Cover",
                    na.value = "transparent", na.translate = FALSE) +
  geom_sf(data = hex_filtered, fill = NA, color = "white",
          linewidth = 0.35, alpha = 0.6) +
  coord_sf(xlim = c(bbox["xmin"], bbox["xmax"]),
           ylim = c(bbox["ymin"], bbox["ymax"]),
           expand = FALSE) +
  theme_void() +
  theme(legend.position = "right", plot.margin = margin(10, 10, 10, 10))

p_detections <- p_base +
  geom_sf(data   = species_sf,
          aes(shape = det_label),
          fill   = "yellow", color = "black",
          size   = 5, stroke = 0.6, alpha = 0.85) +
  scale_shape_manual(name = NULL, values = setNames(21, det_label)) +
  guides(
    fill  = guide_legend(order = 1, title = "Land Cover"),
    shape = guide_legend(order = 2,
                         override.aes = list(fill = "yellow", color = "black",
                                             size = 3, stroke = 0.6))
  )

# Save

ggsave(file.path(output_dir, "figS2_site_detection_map.png"),
       p_detections, width = 10, height = 8, dpi = 300, bg = "white")
cat("Saved: fig5C_site_detection_map.png\n\n")
