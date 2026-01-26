##############################
#
# Land Cover composition by site
# Ian Becker
# January 2026
#
##############################

# This script makes a figure showing land cover composition by site

library(tidyr)
library(dplyr)
library(ggplot2)

# Load site covariates

landcover_vars <- readRDS("site_covariates.rds")

# Summarize mean cover across all sites
overall_comp <- landcover_vars %>%
  summarise(
    Trees = mean(trees_pct, na.rm = TRUE),
    Grassland = mean(grass_pct, na.rm = TRUE),
    Crops = mean(crops_pct, na.rm = TRUE),
    Shrub = mean(shrub_pct, na.rm = TRUE),
    Wetland = mean(flooded_veg_pct, na.rm = TRUE),
    Urban = mean(urban_pct, na.rm = TRUE),
    Bare = mean(bare_pct, na.rm = TRUE)
  ) %>%
  pivot_longer(everything(), names_to = "cover_type", values_to = "percent")

# Set factor order
overall_comp$cover_type <- factor(overall_comp$cover_type,
                                  levels = c("Urban", "Bare", "Crops", "Grassland", "Shrub", "Trees", "Wetland"))

# Plot
overall_fig <- ggplot(overall_comp, aes(x = "", y = percent, fill = cover_type)) +
  geom_col(width = 0.6) +
  coord_flip() +
  scale_fill_manual(
    values = c(
      "Urban" = "#878787",
      "Bare" = "#d4b996",
      "Crops" = "#f7dc6f",
      "Grassland" = "#abebc6",
      "Shrub" = "#73c6b6",
      "Trees" = "#27ae60",
      "Wetland" = "#5dade2"
    ),
    name = "Land Cover"
  ) +
  labs(title = "Mean Land Cover Composition Across Study Sites",
       x = NULL, y = "Percent Cover") +
  theme_minimal() +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )

print(overall_fig)

ggsave("landcover_overall.png", plot = overall_fig,
       width = 10, height = 4, dpi = 300, bg = "white")