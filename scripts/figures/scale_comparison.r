##############################
#
# Scale Comparison Coefficient Plot
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
#   Load and Prep Within-Site Results
####################

within <- read.csv(file.path(input_dir, "within_site_models/within_site_results.csv"))

# Clean parameter names (remove trailing numbers)
within$parameter <- gsub("[0-9]+$", "", within$parameter)

# Summarize within-site by species
within_summary <- within %>%
  filter(parameter != "(Intercept)") %>%
  mutate(sig_pos = lower > 0,
         sig_neg = upper < 0) %>%
  group_by(species, guild, parameter) %>%
  summarise(
    within_mean = mean(Mean),
    within_se = sd(Mean) / sqrt(n()),
    within_n_sites = n(),
    within_n_sig_pos = sum(sig_pos),
    within_n_sig_neg = sum(sig_neg),
    .groups = "drop"
  ) %>%
  mutate(
    within_sig = case_when(
      within_n_sig_pos > 0 & within_n_sig_neg == 0 ~ "positive",
      within_n_sig_neg > 0 & within_n_sig_pos == 0 ~ "negative",
      within_n_sig_pos > 0 & within_n_sig_neg > 0 ~ "mixed",
      TRUE ~ "NS"
    )
  )

####################
#   Load and Prep Site-Level Results
####################

site_results <- readRDS(file.path(input_dir, "model_results/all_results_2026-01-20.rds"))

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
      site_se = sd(vals),
      site_lower = quantile(vals, 0.025),
      site_upper = quantile(vals, 0.975)
    ))
  }
}

# Clean parameter names to match within-site
site_coefs$parameter <- gsub("_pct", "", site_coefs$parameter)

# Add significance
site_coefs <- site_coefs %>%
  mutate(
    site_sig = case_when(
      site_lower > 0 ~ "positive",
      site_upper < 0 ~ "negative",
      TRUE ~ "NS"
    )
  )

####################
#   Merge Site and Within-Site
####################

comparison <- within_summary %>%
  left_join(site_coefs, by = c("species", "parameter")) %>%
  filter(!is.na(site_mean))  # Only keep covariates present in both

# Add guild info
guild_lookup <- within_summary %>% 
  select(species, guild) %>% 
  distinct()

comparison <- comparison %>%
  left_join(guild_lookup, by = "species")

####################
#   Create Comparison Plot
####################

# Focus on key covariates
key_covs <- c("water", "trees", "grass", "flooded_veg", "habitat_diversity", "shrub")

plot_data <- comparison %>%
  filter(parameter %in% key_covs) %>%
  mutate(
    # Clean up labels
    parameter = factor(parameter, 
                       levels = c("water", "flooded_veg", "trees", "shrub", "grass", "habitat_diversity"),
                       labels = c("Water", "Flooded Veg", "Trees", "Shrub", "Grass", "Habitat Diversity")),
    guild = factor(guild.x, 
                   levels = c("waterbird", "specialist", "generalist"),
                   labels = c("Waterbirds", "Specialists", "Generalists"))
  )

plot_data <- plot_data %>%
  mutate(
    site_is_sig = site_sig != "NS",
    within_is_sig = within_sig != "NS" & within_sig != "mixed"
  )

# Main comparison plot
p <- ggplot(plot_data, aes(y = reorder(species, as.numeric(guild)))) +
  
  # Zero line
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  
  # Connecting segment
  geom_segment(aes(x = site_mean, xend = within_mean, yend = species), 
               color = "gray70", linetype = "dotted", linewidth = 0.5) +
  
  # Site-level point (circle)
  geom_point(aes(x = site_mean, fill = site_is_sig), 
             shape = 21, color = "#1f78b4", size = 3.5, stroke = 1.2) +
  
  # Within-site point (triangle)
  geom_point(aes(x = within_mean, fill = within_is_sig), 
             shape = 24, color = "#ff7f00", size = 3.5, stroke = 1.2) +
  
  # Facet by covariate
  facet_wrap(~parameter, scales = "free_x", ncol = 3) +
  
  # Fill scale - significant vs NS
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
    subtitle = "Blue circles = Site-level | Orange triangles = Within-site | Filled = Significant"
  ) +
  
  # Theme
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5, color = "gray40"),
    strip.text = element_text(size = 11, face = "bold"),
    axis.text.y = element_text(size = 8),
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    panel.spacing = unit(1, "lines")
  )

print(p)

ggsave(file.path(output_dir, "scale_comparison_coefficients.png"),
       p, width = 12, height = 10, dpi = 300)

####################
#   Alternative: Grouped by Guild
####################

p_guild <- ggplot(plot_data, aes(y = species)) +
  
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  
  geom_segment(aes(x = site_mean, xend = within_mean, yend = species), 
               color = "gray70", linetype = "dotted", linewidth = 0.5) +
  
  geom_point(aes(x = site_mean), color = "#2d6a4f", size = 3) +
  geom_point(aes(x = within_mean), color = "#ae2012", size = 3) +
  
  facet_grid(guild ~ parameter, scales = "free", space = "free_y") +
  
  labs(
    x = "Coefficient Estimate",
    y = NULL,
    title = "Scale-Dependent Habitat Selection by Guild",
    subtitle = "Green = Site-level | Red = Within-site"
  ) +
  
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5, color = "gray40"),
    strip.text = element_text(size = 9, face = "bold"),
    strip.text.y = element_text(angle = 0),
    axis.text.y = element_text(size = 8),
    panel.grid.minor = element_blank()
  )

print(p_guild)

ggsave(file.path(output_dir, "scale_comparison_by_guild.png"),
       p_guild, width = 14, height = 8, dpi = 300)

cat("\nDone! Saved to", output_dir, "\n")