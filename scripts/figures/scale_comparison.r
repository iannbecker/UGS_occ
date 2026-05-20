##############################
#
# Scale Comparison Coefficient Plot
# All Within-Site Points (not averaged)
# Ian Becker
# January 2026
#
##############################

library(dplyr)
library(ggplot2)
library(tidyr)

# Set paths
input_dir <- "/Users/ianbecker/Desktop/project_code/UGS_occ/data"
output_dir <- "/Users/ianbecker/Desktop/project_code/UGS_occ/figures_tables"

####################
#   Load Within-Site Results (keep all individual site estimates)
####################

within <- read.csv(file.path(input_dir, "within_site_models/within_site_results.csv"))

# Clean parameter names (remove trailing numbers)
within$parameter <- gsub("[0-9]+$", "", within$parameter)

# Filter to non-intercepts and add significance flag
within_all <- within %>%
  filter(parameter != "(Intercept)") %>%
  mutate(
    within_is_sig = lower > 0 | upper < 0
  )

####################
#   Load and Prep Site-Level Results
####################

site_results <- readRDS(file.path(input_dir, "model_results/all_results_2026-03-19.rds"))

# Extract site-level coefficients
site_coefs <- data.frame()

for (sp in names(site_results)) {
  if (!site_results[[sp]]$success) next
  
  model <- site_results[[sp]]$model
  beta <- as.matrix(model$beta.samples)
  
  for (p in colnames(beta)) {
    vals <- beta[, p]
    site_coefs <- rbind(site_coefs, data.frame(
      species = site_results[[sp]]$species,
      parameter = p,
      site_mean = mean(vals),
      site_lower = quantile(vals, 0.025),
      site_upper = quantile(vals, 0.975)
    ))
  }
}

# Clean parameter names to match within-site
site_coefs$parameter <- gsub("_pct", "", site_coefs$parameter)

# Add significance
site_coefs <- site_coefs %>%
  filter(parameter != "(Intercept)") %>%
  mutate(
    site_is_sig = site_lower > 0 | site_upper < 0
  )

####################
#   Prep Plot Data
####################

# Focus on key covariates
key_covs <- c("water", "trees", "grass", "flooded_veg", "habitat_diversity", "shrub")

# Filter within-site data
within_plot <- within_all %>%
  filter(parameter %in% key_covs) %>%
  mutate(
    parameter = factor(parameter, 
                       levels = c("water", "flooded_veg", "trees", "shrub", "grass", "habitat_diversity"),
                       labels = c("Water", "Flooded Veg", "Trees", "Shrub", "Grass", "Habitat Diversity"))
  )

# Filter site-level data
site_plot <- site_coefs %>%
  filter(parameter %in% key_covs) %>%
  mutate(
    parameter = factor(parameter, 
                       levels = c("water", "flooded_veg", "trees", "shrub", "grass", "habitat_diversity"),
                       labels = c("Water", "Flooded Veg", "Trees", "Shrub", "Grass", "Habitat Diversity"))
  )

####################
#   Create Plot
####################

p <- ggplot() +
  
  # Zero line
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  
  # Within-site points (triangles)
  geom_point(data = within_plot,
             aes(x = Mean, y = reorder(species, species), fill = within_is_sig),
             shape = 24, color = "#ff7f00", size = 2, stroke = 0.8, alpha = 0.7,
             position = position_jitter(height = 0.15, width = 0, seed = 123)) +
  
  # Site-level points (circles)
  geom_point(data = site_plot,
             aes(x = site_mean, y = species, fill = site_is_sig),
             shape = 21, color = "#1f78b4", size = 3.5, stroke = 1.2) +
  
  # Facet by covariate
  facet_wrap(~parameter, scales = "free_x", ncol = 3) +
  
  # Fill scale
  scale_fill_manual(
    values = c("TRUE" = "black", "FALSE" = "white"),
    name = "Significant",
    labels = c("No", "Yes")
  ) +
  
  # Labels
  labs(
    x = "Coefficient Estimate (standardized)",
    y = NULL,
    title = "Scale-Dependent Habitat Selection",
    subtitle = "Blue circles = Site-level | Orange triangles = Within-site (each site shown) | Filled = Significant"
  ) +
  
  # Theme
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5, color = "gray40"),
    strip.text = element_text(size = 11, face = "bold"),
    axis.text.y = element_text(size = 6),
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    panel.spacing = unit(1, "lines")
  )

print(p)

ggsave(file.path(output_dir, "scale_comparison_all_within_sites.png"),
       p, width = 14, height = 20, dpi = 300)

####################
#   Summary stats
####################

cat("\n=== SUMMARY ===\n")
cat("Species:", length(unique(within_plot$species)), "\n")
cat("Within-site estimates:", nrow(within_plot), "\n")
cat("Site-level estimates:", nrow(site_plot), "\n")

cat("\nDone! Saved to", output_dir, "\n")