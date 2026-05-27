##############################
#
# Figure 2: Landscape Occupancy
# Ian Becker
# May 2026
#
##############################

# This script creates Figure 2 in the manuscript,
# which shows occupancy response curves for all 38 species across all 8 landscape-level covariates.

library(dplyr)
library(ggplot2)
library(patchwork)

setwd("~/Desktop/project_code/UGS_occ/data")

output_dir <- "/Users/ianbecker/Desktop/project_code/UGS_occ/figures_tables"

# DATA PREP ------------------------------

# Change label and order of covariates for plotting

covariate_info <- list(
  list(param = "grass_pct",          label = "Grass Cover (%)"),
  list(param = "trees_pct",          label = "Tree Cover (%)"),
  list(param = "shrub_pct",          label = "Shrub Cover (%)"),
  list(param = "flooded_veg_pct",    label = "Flooded Vegetation Cover (%)"),
  list(param = "water_pct",          label = "Water Cover (%)"),
  list(param = "crops_pct",          label = "Crop Cover (%)"),
  list(param = "habitat_diversity",  label = "Habitat Diversity"),
  list(param = "log_area",           label = "Site Area (log)")
)

# Number of points to predict along each curve

n_points <- 100

# Load in data

site_covs <- read.csv("site_covariates.csv")

# Calculate mean and SD for each covariate

cov_scales <- site_covs %>%
  select(-site_id, -bare_pct, -urban_pct, -area_ha) %>%
  pivot_longer(everything(), names_to = "param") %>%
  group_by(param) %>%
  summarise(
    cov_mean = mean(value, na.rm = TRUE),
    cov_sd   = sd(value,   na.rm = TRUE),
    x_min    = min(value,  na.rm = TRUE),
    x_max    = max(value,  na.rm = TRUE),
    .groups = "drop"
  )

# Colors for each covariate

cov_colors <- c(
  "Grass Cover (%)"              = "#a8d08d",
  "Tree Cover (%)"               = "#2d6a2d",
  "Shrub Cover (%)"              = "#8c6d31",
  "Flooded Vegetation Cover (%)" = "#7b2d8b",
  "Water Cover (%)"              = "#4393c3",
  "Crop Cover (%)"               = "#d9ef8b",
  "Habitat Diversity"            = "#d73027",
  "Site Area (log)"              = "#636363"
)

# Load in models

all_results <- readRDS("model_results/all_results_2026-03-19.rds")
successful  <- all_results[sapply(all_results, function(x) x$success)]
cat("Successful models:", length(successful), "\n\n")

# HELPER FUNCTION TO GENERATE OCCUPANCY CURVE FOR EACH SPECIES X COVARIATE  ------------------------------

get_curve <- function(model, focal_param, cov_mean, cov_sd, 
                      x_min, x_max, n_points = 100) {
  
  beta_samples <- as.matrix(model$beta.samples)
  param_names  <- colnames(beta_samples)
  
  param_names_clean <- gsub("[0-9]+$", "", param_names)
  focal_idx <- which(param_names_clean == focal_param)
  int_idx   <- which(param_names_clean == "(Intercept)")
  
  if (length(focal_idx) == 0) return(NULL)
  
  # x sequence on original scale, clipped to observed range
  
  x_original <- seq(x_min, x_max, length.out = n_points)
  
  # Scale for prediction
  
  x_scaled <- (x_original - cov_mean) / cov_sd
  
  n_samples   <- nrow(beta_samples)
  pred_matrix <- matrix(NA, nrow = n_samples, ncol = n_points)
  
  for (i in 1:n_samples) {
    linear_pred       <- beta_samples[i, int_idx] + beta_samples[i, focal_idx] * x_scaled
    pred_matrix[i, ]  <- plogis(linear_pred)
  }
  
  data.frame(
    x    = x_original,
    mean = apply(pred_matrix, 2, mean)
  )
}

# GENERATE ALL OCCUPANCY CURVES  ------------------------------

all_curves <- data.frame()

for (sp_name in names(successful)) {
  
  sp    <- successful[[sp_name]]$species
  model <- successful[[sp_name]]$model
  
  for (cov in covariate_info) {
    
    cov_row <- cov_scales[cov_scales$param == cov$param, ]
    if (nrow(cov_row) == 0) next
    
    curve <- get_curve(model, cov$param,
                       cov_mean = cov_row$cov_mean,
                       cov_sd   = cov_row$cov_sd,
                       x_min    = cov_row$x_min,
                       x_max    = cov_row$x_max,
                       n_points = n_points)
    
    if (!is.null(curve)) {
      curve$species   <- sp
      curve$parameter <- cov$param
      curve$label     <- cov$label
      all_curves <- rbind(all_curves, curve)
    }
  }
}

cat("Curves generated for", length(unique(all_curves$species)), "species\n\n")

# CREATE PLOT ------------------------------

# Calulate mean curve for each covariate across species

mean_curves <- all_curves %>%
  group_by(parameter, label, x) %>%
  summarise(mean_occ = mean(mean, na.rm = TRUE), .groups = "drop")

# Set factor levels for plotting

param_order <- sapply(covariate_info, function(x) x$param)
label_order <- sapply(covariate_info, function(x) x$label)

all_curves$label  <- factor(all_curves$label,  levels = label_order)
mean_curves$label <- factor(mean_curves$label, levels = label_order)


# Plot

p <- ggplot() +
  
  # Individual species curves 
  
  geom_line(
    data = all_curves,
    aes(x = x, y = mean, group = species),
    color     = "gray80",
    linewidth = 0.3,
    alpha     = 0.7
  ) +
  
  # Mean curve — colored by covariate
  
  geom_line(
    data = mean_curves,
    aes(x = x, y = mean_occ, color = label),
    linewidth = 1.3
  ) +
  
  scale_color_manual(values = cov_colors, guide = "none") +
  
  # Facet by covariate
  
  facet_wrap(~ label, ncol = 4, scales = "free_x", strip.position = "bottom") +
  
  scale_y_continuous(
    limits = c(0, 1),
    breaks = c(0, 0.5, 1),
    labels = c("0", "0.5", "1")
  ) +
  
  scale_x_continuous(n.breaks = 4) +
  
  labs(
    x = NULL,
    y = "Occupancy Probability"
  ) +
  
  theme_classic() +
  theme(
    strip.text         = element_text(size = 12, face = "bold"),
    strip.background   = element_blank(),
    strip.placement    = "outside",
    axis.text          = element_text(size = 9),
    axis.title.y       = element_text(size = 10),
    panel.grid.major.y = element_line(color = "gray90", linewidth = 0.3),
    panel.spacing      = unit(0.8, "lines")
  )

ggsave(
  file.path(output_dir, "Figure1.png"),
  p, width = 12, height = 6, dpi = 200, bg = "white"
)

cat("Saved: Figure1.png\n")
