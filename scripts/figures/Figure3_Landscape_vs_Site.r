##############################
#
# Figure 3: Landscape vs. Site-level Coefficient Comparison
# Ian Becker
# May 2026
#
##############################

# This script is used to make Figure 3 in the manuscript;
# Plot comparing site-level and landscape-level covariates for each species-covariate combination.

library(dplyr)
library(ggplot2)

setwd("~/Desktop/project_code/UGS_occ/data")

output_dir <- "/Users/ianbecker/Desktop/project_code/UGS_occ/figures_tables"

# ============================================================================
# 1. SETUP AND LOAD DATA
# ============================================================================

# Load in landscape-level results

all_results <- readRDS("model_results_ebird_checklist/all_results_ebird_checklist_2026-08-11.rds")

# Extract landscape-level coefficients and credible intervals

site_coefs <- data.frame()
for (sp in names(all_results)) {
  if (!all_results[[sp]]$success) next
  model <- all_results[[sp]]$model
  beta  <- as.matrix(model$beta.samples)
  for (p in colnames(beta)) {
    vals <- beta[, p]
    site_coefs <- rbind(site_coefs, data.frame(
      species    = all_results[[sp]]$species,
      parameter  = p,
      land_mean  = mean(vals),
      land_lower = quantile(vals, 0.025),
      land_upper = quantile(vals, 0.975),
      stringsAsFactors = FALSE
    ))
  }
}

# Clean up parameter names 

site_coefs <- site_coefs %>%
  filter(parameter != "(Intercept)") %>%
  mutate(
    parameter = gsub("[0-9]+$", "", parameter),
    parameter = gsub("_pct", "", parameter)
  )

# Load in site-level results

within_coefs <- read.csv("within_site_models_gbif/within_site_results_gbif.csv")

# Clean up parameter names to match landscape level

within_coefs <- within_coefs %>%
  filter(parameter != "(Intercept)") %>%
  mutate(
    parameter = gsub("[0-9]+$", "", parameter),
    parameter = gsub("_pct", "", parameter)
  )

# Average within-site coefficients across sites for each species-covariate combination

within_avg <- within_coefs %>%
  group_by(species, parameter) %>%
  summarise(
    within_mean = mean(Mean, na.rm = TRUE),
    n_sites     = n(),
    .groups = "drop"
  )

# Combine landscape and site-level data

combined <- site_coefs %>%
  inner_join(within_avg, by = c("species", "parameter")) %>%
  filter(!is.na(within_mean))


# Fix labels

cov_labels <- c(
  "trees"             = "Tree Cover",
  "grass"             = "Grass Cover",
  "shrub"             = "Shrub Cover",
  "flooded_veg"       = "Flooded Veg",
  "crops"             = "Crop Cover",
  "water"             = "Water Cover",
  "habitat_diversity" = "Habitat Diversity",
  "log_area"          = "Site Area"
)

combined <- combined %>%
  mutate(
    covariate_label = ifelse(parameter %in% names(cov_labels),
                             cov_labels[parameter], parameter))

# ============================================================================
# 2. PLOTTING
# ============================================================================

# Setup colors for each covariate

cov_colors <- c(
  "Tree Cover"        = "#2d6a2d",
  "Grass Cover"       = "#a8d08d",
  "Shrub Cover"       = "#8c6d31",
  "Flooded Veg"       = "#7b2d8b",
  "Crop Cover"        = "#d9ef8b",
  "Water Cover"       = "#4393c3",
  "Habitat Diversity" = "#d73027"
)

# Plot

p <- ggplot(combined,
            aes(x = land_mean, y = within_mean,
                size = n_sites, color = covariate_label)) +
  
  # Reference lines
  
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.5) +
  
  # 1:1 relationship line
  
  geom_abline(slope = 1, intercept = 0, linetype = "dotted",
              color = "gray60", linewidth = 0.5) +
  
  # Points
  
  geom_point(alpha = 0.7, stroke = 0.3) +
  
  scale_color_manual(
    values = cov_colors,
    name   = "Covariate",
    guide  = guide_legend(override.aes = list(size = 6))
  ) +
  
  scale_size_continuous(
    name   = "Viable sites (n)",
    range  = c(2, 8),
    breaks = c(1, 3, 5, 10)
  ) +
  
  labs(
    x = "Landscape-level coefficient (standardized)",
    y = "Site-level coefficient (standardized)"
  ) +
  
  scale_x_continuous(limits = c(-3, 3)) +
  scale_y_continuous(limits = c(-3, 3)) +
  
  theme_classic() +
  theme(
    axis.title       = element_text(size = 15),
    axis.text        = element_text(size = 13),
    legend.position  = "right",
    legend.title     = element_text(size = 11),
    legend.text      = element_text(size = 10),
    panel.background = element_rect(fill = "white"),
    plot.background  = element_rect(fill = "white", color = NA)
  )

ggsave(
  file.path(output_dir, "Figure3.png"),
  p, width = 10, height = 6, dpi = 200, bg = "white"
)

cat("Saved: Figure3.png\n")
