##############################
#
# Bubble Plot — Landscape vs Within-Site Coefficients
# Ian Becker
# May 2026
#
##############################

# X = landscape-level coefficient
# Y = within-site mean coefficient
# Bubble size = number of viable sites
# Color = covariate
# Quadrants: I & III = concordant, II & IV = discordant (scale flip)

library(dplyr)
library(ggplot2)

setwd("~/Desktop/project_code/UGS_occ/data")

output_dir <- "/Users/ianbecker/Desktop/project_code/UGS_occ/figures_tables"

# ── Load data ──────────────────────────────────────────────────────────────────
cat("Loading data...\n")

all_results <- readRDS("model_results/all_results_2026-03-19.rds")

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

site_coefs <- site_coefs %>%
  filter(parameter != "(Intercept)") %>%
  mutate(
    parameter = gsub("[0-9]+$", "", parameter),
    parameter = gsub("_pct", "", parameter)
  )

# Within-site — average across sites per species
within_coefs <- read.csv("within_site_models/within_site_results.csv")
within_coefs <- within_coefs %>%
  filter(parameter != "(Intercept)") %>%
  mutate(
    parameter = gsub("[0-9]+$", "", parameter),
    parameter = gsub("_pct", "", parameter)
  )

within_avg <- within_coefs %>%
  group_by(species, parameter) %>%
  summarise(
    within_mean = mean(Mean, na.rm = TRUE),
    n_sites     = n(),
    .groups = "drop"
  )

# ── Join landscape and within-site ────────────────────────────────────────────
combined <- site_coefs %>%
  inner_join(within_avg, by = c("species", "parameter")) %>%
  filter(!is.na(within_mean))

cat("Species-covariate combinations:", nrow(combined), "\n\n")

# ── Covariate display labels ───────────────────────────────────────────────────
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
                             cov_labels[parameter], parameter),
    # Quadrant classification
    quadrant = case_when(
      land_mean > 0 & within_mean > 0 ~ "Concordant positive",
      land_mean < 0 & within_mean < 0 ~ "Concordant negative",
      land_mean > 0 & within_mean < 0 ~ "Discordant (+ to -)",
      land_mean < 0 & within_mean > 0 ~ "Discordant (- to +)",
      TRUE ~ "Near zero"
    )
  )

# ── Color palette for covariates ───────────────────────────────────────────────
cov_colors <- c(
  "Tree Cover"        = "#2d6a2d",
  "Grass Cover"       = "#a8d08d",
  "Shrub Cover"       = "#8c6d31",
  "Flooded Veg"       = "#7b2d8b",
  "Crop Cover"        = "#d9ef8b",
  "Water Cover"       = "#4393c3",
  "Habitat Diversity" = "#d73027"
)

# ── Plot ───────────────────────────────────────────────────────────────────────
p <- ggplot(combined,
            aes(x = land_mean, y = within_mean,
                size = n_sites, color = covariate_label)) +
  
  # Reference lines
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.5) +
  
  # 1:1 concordance line
  geom_abline(slope = 1, intercept = 0, linetype = "dotted",
              color = "gray60", linewidth = 0.5) +
  
  # Bubbles
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
    axis.title       = element_text(size = 11),
    axis.text        = element_text(size = 10),
    legend.position  = "right",
    legend.title     = element_text(size = 9),
    legend.text      = element_text(size = 8),
    panel.background = element_rect(fill = "white"),
    plot.background  = element_rect(fill = "white", color = NA)
  )

ggsave(
  file.path(output_dir, "bubble_plot_scale_comparison.png"),
  p, width = 10, height = 6, dpi = 200, bg = "white"
)

cat("Saved: bubble_plot_scale_comparison.png\n")

# ── Summary by quadrant ────────────────────────────────────────────────────────
cat("\n=== QUADRANT SUMMARY ===\n\n")
quad_summary <- combined %>%
  count(quadrant) %>%
  mutate(pct = round(n / sum(n) * 100, 1)) %>%
  arrange(desc(n))
print(quad_summary)

cat("\n=== BY COVARIATE ===\n\n")
cov_summary <- combined %>%
  group_by(covariate_label, quadrant) %>%
  summarise(n = n(), .groups = "drop") %>%
  arrange(covariate_label, desc(n))
print(cov_summary)