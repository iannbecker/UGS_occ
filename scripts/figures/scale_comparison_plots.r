##############################
#
# Per-Species Scale Comparison Plots
# Ian Becker
# May 2026
#
##############################

# Saves two folders of PNGs per species:
#   figures_tables/scale_plots_averaged/   <- one triangle per species (mean across sites)
#   figures_tables/scale_plots_allsites/   <- one triangle per viable site (jittered)
#
# Built directly on the confirmed working plot code from March 2026.
# Key design: shape 21 circles (landscape-level, blue) + shape 24 triangles (site-level, orange)
#             fill = black (significant) or white (non-significant)

library(dplyr)
library(ggplot2)

# ── Paths ──────────────────────────────────────────────────────────────────────
input_dir  <- "/Users/ianbecker/Desktop/project_code/UGS_occ/data"
output_dir <- "/Users/ianbecker/Desktop/project_code/UGS_occ/figures_tables"

dir.create(file.path(output_dir, "scale_plots_averaged"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(output_dir, "scale_plots_allsites"), showWarnings = FALSE, recursive = TRUE)

# ── Covariate display name lookup ─────────────────────────────────────────────
# Maps internal parameter names (after _pct stripping) to display labels
covariate_labels <- c(
  "trees"              = "Tree Cover (%)",
  "grass"              = "Grass Cover (%)",
  "shrub"              = "Shrub Cover (%)",
  "flooded_veg"        = "Flooded Vegetation Cover (%)",
  "crops"              = "Crop Cover (%)",
  "habitat_diversity"  = "Habitat Diversity",
  "log_area"           = "Site Area (log)",
  "water"              = "Water Cover (%)"
)

# Display order (internal names — will be mapped to labels for plotting)
covariate_order <- c(
  "trees", "grass", "shrub",
  "flooded_veg", "crops",
  "habitat_diversity", "log_area", "water"
)

# ── Load within-site results ───────────────────────────────────────────────────
cat("Loading within-site results...\n")
within <- read.csv(file.path(input_dir, "within_site_models/within_site_results.csv"))

# Clean parameter names (remove trailing numbers, strip _pct for display)
within$parameter <- gsub("[0-9]+$", "", within$parameter)
within$parameter <- gsub("_pct", "", within$parameter)

within_all <- within %>%
  filter(parameter != "(Intercept)") %>%
  mutate(within_is_sig = lower > 0 | upper < 0)

# ── Load and extract site-level results ───────────────────────────────────────
cat("Extracting site-level coefficients...\n")
site_results <- readRDS(file.path(input_dir, "model_results/all_results_2026-03-19.rds"))

site_coefs <- data.frame()

for (sp in names(site_results)) {
  if (!site_results[[sp]]$success) next
  
  model <- site_results[[sp]]$model
  beta  <- as.matrix(model$beta.samples)
  
  for (p in colnames(beta)) {
    vals <- beta[, p]
    site_coefs <- rbind(site_coefs, data.frame(
      species    = site_results[[sp]]$species,
      parameter  = p,
      site_mean  = mean(vals),
      site_lower = quantile(vals, 0.025),
      site_upper = quantile(vals, 0.975),
      stringsAsFactors = FALSE
    ))
  }
}

# Match parameter name cleaning to within-site
site_coefs$parameter <- gsub("[0-9]+$", "", site_coefs$parameter)
site_coefs$parameter <- gsub("_pct", "", site_coefs$parameter)

site_coefs <- site_coefs %>%
  filter(parameter != "(Intercept)") %>%
  mutate(site_is_sig = site_lower > 0 | site_upper < 0)

cat("Landscape-level species:", length(unique(site_coefs$species)), "\n")
cat("Site-level species:", length(unique(within_all$species)), "\n\n")

# ── Helper: order and label parameters ────────────────────────────────────────
order_params <- function(params) {
  ordered <- intersect(covariate_order, params)
  extras  <- setdiff(params, covariate_order)
  c(ordered, extras)
}

apply_labels <- function(params) {
  # Map internal names to display labels, fall back to raw name if not found
  ifelse(params %in% names(covariate_labels),
         covariate_labels[params],
         params)
}

# ── Shared theme ───────────────────────────────────────────────────────────────
plot_theme <- theme_minimal() +
  theme(
    plot.title       = element_text(size = 13, face = "italic", hjust = 0.5),
    plot.subtitle    = element_blank(),
    strip.text       = element_text(size = 11, face = "bold"),
    axis.text.y      = element_text(size = 9),
    axis.text.x      = element_text(size = 9),
    axis.title       = element_text(size = 10),
    legend.position  = "bottom",
    legend.title     = element_text(size = 9),
    legend.text      = element_text(size = 9),
    legend.spacing.x = unit(0.3, "cm"),
    panel.grid.minor = element_blank()
  )

# ── Helper: build one species plot ────────────────────────────────────────────
make_species_plot <- function(sp, site_df, within_df) {
  
  s <- site_df  %>% filter(species == sp)
  w <- within_df  # already species-filtered before calling
  
  # Get ordered internal parameter names
  all_params   <- union(s$parameter, w$parameter)
  param_levels <- order_params(all_params)
  
  # Apply display labels
  param_labels <- apply_labels(param_levels)
  
  # Recode parameters to display labels
  s <- s %>%
    mutate(parameter = factor(apply_labels(parameter), levels = param_labels)) %>%
    filter(!is.na(parameter))
  w <- w %>%
    mutate(parameter = factor(apply_labels(parameter), levels = param_labels)) %>%
    filter(!is.na(parameter))
  
  ggplot() +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    
    # Site-level triangles (shape 24 = filled upward triangle, orange)
    geom_point(
      data = w,
      aes(x = Mean, y = parameter, fill = within_is_sig,
          shape = "Site-level"),
      color = "#ff7f00", size = 2.8, stroke = 1,
      alpha = 0.75,
      position = position_jitter(height = 0.15, width = 0, seed = 42)
    ) +
    
    # Landscape-level CI bar (blue)
    geom_errorbarh(
      data = s,
      aes(xmin = site_lower, xmax = site_upper, y = parameter),
      color = "#1f78b4", height = 0.15, linewidth = 0.7
    ) +
    
    # Landscape-level circle (shape 21 = filled circle, blue)
    geom_point(
      data = s,
      aes(x = site_mean, y = parameter, fill = site_is_sig,
          shape = "Landscape-level"),
      color = "#1f78b4", size = 4, stroke = 1.2
    ) +
    
    # Shape scale: circle for landscape, triangle for site
    scale_shape_manual(
      name   = NULL,
      values = c("Landscape-level" = 21, "Site-level" = 24),
      guide  = guide_legend(
        override.aes = list(
          color = c("#1f78b4", "#ff7f00"),
          fill  = c("gray50", "gray50"),
          size  = c(4, 3)
        ),
        order = 1
      )
    ) +
    
    # Fill scale: significance
    scale_fill_manual(
      name   = "Significant (95% CI \u2260 0)",
      values = c("TRUE" = "black", "FALSE" = "white"),
      labels = c("No", "Yes"),
      guide  = guide_legend(
        override.aes = list(
          shape = 21,
          color = "gray40",
          size  = 3.5
        ),
        order = 2
      )
    ) +
    
    scale_y_discrete(limits = rev(param_labels)) +
    
    labs(
      title = sp,
      x     = "Coefficient estimate (standardized)",
      y     = NULL
    ) +
    
    plot_theme
}

# ── Helper: landscape-level only plot (no site-level data) ────────────────────
make_siteonly_plot <- function(sp, site_df) {
  
  s <- site_df %>% filter(species == sp)
  param_levels <- order_params(unique(s$parameter))
  param_labels <- apply_labels(param_levels)
  
  s <- s %>%
    mutate(parameter = factor(apply_labels(parameter), levels = param_labels)) %>%
    filter(!is.na(parameter))
  
  ggplot(s) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    geom_errorbarh(
      aes(xmin = site_lower, xmax = site_upper, y = parameter),
      color = "#1f78b4", height = 0.15, linewidth = 0.7
    ) +
    geom_point(
      aes(x = site_mean, y = parameter, fill = site_is_sig,
          shape = "Landscape-level"),
      color = "#1f78b4", size = 4, stroke = 1.2
    ) +
    scale_shape_manual(
      name   = NULL,
      values = c("Landscape-level" = 21),
      guide  = guide_legend(
        override.aes = list(
          color = "#1f78b4",
          fill  = "gray50",
          size  = 4
        ),
        order = 1
      )
    ) +
    scale_fill_manual(
      name   = "Significant (95% CI \u2260 0)",
      values = c("TRUE" = "black", "FALSE" = "white"),
      labels = c("No", "Yes"),
      guide  = guide_legend(
        override.aes = list(
          shape = 21,
          color = "gray40",
          size  = 3.5
        ),
        order = 2
      )
    ) +
    scale_y_discrete(limits = rev(param_labels)) +
    labs(
      title = sp,
      x     = "Coefficient estimate (standardized)",
      y     = NULL
    ) +
    plot_theme
}

# ── Species sets ───────────────────────────────────────────────────────────────
species_site     <- unique(site_coefs$species)
species_within   <- unique(within_all$species)
species_both     <- intersect(species_site, species_within)
species_siteonly <- setdiff(species_site, species_within)

cat("Species with both scales:", length(species_both), "\n")
cat("Species with landscape-level only:", length(species_siteonly), "\n\n")

# ── VERSION 1: Averaged within-site ───────────────────────────────────────────
cat("=== SAVING AVERAGED VERSION ===\n")

within_avg <- within_all %>%
  group_by(species, parameter) %>%
  summarise(
    Mean          = mean(Mean,  na.rm = TRUE),
    lower         = mean(lower, na.rm = TRUE),
    upper         = mean(upper, na.rm = TRUE),
    within_is_sig = mean(lower) > 0 | mean(upper) < 0,
    n_sites       = n(),
    .groups = "drop"
  )

for (sp in species_both) {
  w_sp <- within_avg %>% filter(species == sp)
  p    <- make_species_plot(sp = sp, site_df = site_coefs, within_df = w_sp)
  fname <- paste0(gsub("[^A-Za-z0-9]", "_", sp), "_averaged.png")
  ggsave(file.path(output_dir, "scale_plots_averaged", fname),
         p, width = 6.5, height = 5.5, dpi = 200, bg = "white")
  cat("  Saved:", fname, "\n")
}

for (sp in species_siteonly) {
  p     <- make_siteonly_plot(sp, site_coefs)
  fname <- paste0(gsub("[^A-Za-z0-9]", "_", sp), "_averaged.png")
  ggsave(file.path(output_dir, "scale_plots_averaged", fname),
         p, width = 6.5, height = 5.5, dpi = 200, bg = "white")
  cat("  Saved (landscape-only):", fname, "\n")
}

# ── VERSION 2: All sites shown (jittered) ─────────────────────────────────────
cat("\n=== SAVING ALL-SITES VERSION ===\n")

for (sp in species_both) {
  w_sp <- within_all %>% filter(species == sp)
  p    <- make_species_plot(sp = sp, site_df = site_coefs, within_df = w_sp)
  fname <- paste0(gsub("[^A-Za-z0-9]", "_", sp), "_allsites.png")
  ggsave(file.path(output_dir, "scale_plots_allsites", fname),
         p, width = 6.5, height = 5.5, dpi = 200, bg = "white")
  cat("  Saved:", fname, "\n")
}

for (sp in species_siteonly) {
  fname_src <- file.path(output_dir, "scale_plots_averaged",
                         paste0(gsub("[^A-Za-z0-9]", "_", sp), "_averaged.png"))
  fname_dst <- file.path(output_dir, "scale_plots_allsites",
                         paste0(gsub("[^A-Za-z0-9]", "_", sp), "_allsites.png"))
  file.copy(fname_src, fname_dst, overwrite = TRUE)
  cat("  Copied (landscape-only):", basename(fname_dst), "\n")
}

# ── Summary ───────────────────────────────────────────────────────────────────
cat("\n=== COMPLETE ===\n")
cat("Averaged plots:  scale_plots_averaged/ (", length(species_site), " species)\n")
cat("All-sites plots: scale_plots_allsites/ (", length(species_site), " species)\n")
cat("  With both scales:", length(species_both), "\n")
cat("  Landscape-level only:", length(species_siteonly), "\n")
