##############################
#
# Figure 4: Bar plot
# Ian Becker
# May 2026
#
##############################

# This script creates both parts of Figure 4 in the manuscript
# This includes:
#  - Bar chart of mean naive occupancy (proportion of area detected) per species
#  - Donut chart showing overall mean % area occupied vs unoccupied across all species

library(dplyr)
library(ggplot2)

setwd("~/Desktop/project_code/UGS_occ/data")

output_dir <- "/Users/ianbecker/Desktop/project_code/UGS_occ/figures_tables"

# LOAD DATA AND PREP------------------------------

# Build summary from individual models

model_dir   <- "within_site_models"
model_files <- list.files(model_dir, pattern = "\\.rds$", full.names = TRUE)
model_files <- model_files[!grepl("all_results", model_files)]

cat("Found", length(model_files), "model files\n")

summary_df <- data.frame()

for (f in model_files) {
  tryCatch({
    result <- readRDS(f)
    row <- data.frame(
      species          = result$species,
      site_id          = result$site_id,
      n_cells          = result$n_cells,
      n_cells_detected = result$n_cells_detected,
      naive_occ        = result$naive_occ,
      stringsAsFactors = FALSE
    )
    summary_df <- rbind(summary_df, row)
  }, error = function(e) {
    cat("Error:", f, "\n")
  })
}


# Calculate species-level summary statistics

species_summary <- summary_df %>%
  group_by(species) %>%
  summarise(
    mean_naive_occ   = mean(naive_occ, na.rm = TRUE),
    sd_naive_occ     = sd(naive_occ, na.rm = TRUE),
    se_naive_occ     = sd(naive_occ, na.rm = TRUE) / sqrt(n()),
    n_sites          = n(),
    mean_cells       = round(mean(n_cells, na.rm = TRUE), 1),
    mean_detected    = round(mean(n_cells_detected, na.rm = TRUE), 1),
    .groups = "drop"
  ) %>%
  arrange(mean_naive_occ) %>%
  mutate(
    species     = factor(species, levels = species),
    pct_label   = paste0(round(mean_naive_occ * 100, 1), "%")
  )


# PLOT BAR CHART ------------------------------

# Create bar plot

p <- ggplot(species_summary,
            aes(x = mean_naive_occ, y = species)) +
  
  # Bars
  
  geom_col(
    fill  = "#1f78b4",
    alpha = 0.85,
    width = 0.7
  ) +
  
  # Error bars (SE across sites)
  
  geom_errorbar(
    aes(xmin = pmax(mean_naive_occ - se_naive_occ, 0),
        xmax = pmin(mean_naive_occ + se_naive_occ, 1)),
    width     = 0.3,
    color     = "gray30",
    linewidth = 0.5
  ) +
  
  # Reference line at 50%
  
  geom_vline(xintercept = 0.5, linetype = "dashed",
             color = "gray50", linewidth = 0.5) +
  
  scale_x_continuous(
    labels = scales::percent_format(accuracy = 1),
    limits = c(0, 1),
    expand = c(0.01, 0)
  ) +
  
  scale_y_discrete(limits = rev(levels(species_summary$species))) +
  
  labs(
    x = "Mean proportion of area detected within site (%)",
    y = NULL
  ) +
  
  theme_classic() +
  theme(
    axis.text.y      = element_text(size = 8, face = "italic"),
    axis.text.x      = element_text(size = 9),
    axis.title.x     = element_text(size = 10),
    panel.grid.major.x = element_line(color = "gray90", linewidth = 0.4)
  )

# Height based on number of species

n_species   <- nrow(species_summary)
plot_height <- max(7, n_species * 0.28)

ggsave(
  file.path(output_dir, "within_site_space_use.png"),
  p, width = 8, height = plot_height, dpi = 200, bg = "white"
)

cat("Saved: within_site_space_use.png\n\n")

# PLOT MEAN PERCENTAGE ------------------------------

# Overall mean percentage of used area

overall_mean <- mean(species_summary$mean_naive_occ)

# Create data frame

percentage_df <- data.frame(
  category = c("Occupied", "Unoccupied"),
  value    = c(overall_mean, 1 - overall_mean),
  label    = c(paste0(round(overall_mean * 100, 1), "%"),
               paste0(round((1 - overall_mean) * 100, 1), "%"))
)

# Calculate positions for labels

percentage_df$ymax <- cumsum(percentage_df$value)
percentage_df$ymin <- c(0, head(percentage_df$ymax, -1))
percentage_df$label_pos <- (percentage_df$ymin + percentage_df$ymax) / 2

# Plot

p_percentage <- ggplot(percentage_df, aes(ymax = ymax, ymin = ymin,
                                xmax = 4, xmin = 2.5,
                                fill = category)) +
  geom_rect() +
  
  scale_fill_manual(
    values = c("Occupied"   = "#1f78b4",
               "Unoccupied" = "gray85"),
    name   = NULL
  ) +
  
  coord_polar(theta = "y") +
  xlim(c(0, 5)) +
  
  labs(title = "Mean Within-Site\nSpace Use") +
  
  theme_void() +
  theme(
    plot.title      = element_text(size = 11, face = "bold",
                                   hjust = 0.5, vjust = 0),
    legend.position = "bottom",
    legend.text     = element_text(size = 10)
  )

ggsave(
  file.path(output_dir, "within_site_space_use_percentage.png"),
  p_percentage, width = 4, height = 4, dpi = 300, bg = "transparent"
)

cat("Saved: within_site_space_use_percentage.png\n")
