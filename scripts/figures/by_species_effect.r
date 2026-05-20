##############################
#
# Per-Species Marginal Effect Figures
# 7 covariates × 15 species
#
##############################

# This script is used to make marginal effect plots by species

library(ggplot2)
library(dplyr)
library(patchwork)

# Set paths

input_dir <- "/Users/ianbecker/Desktop/project_code/UGS_occ/data"
output_dir <- "/Users/ianbecker/Desktop/project_code/UGS_occ/figures_tables"


####################
#   Marginal effect function
####################

plot_marginal_effect <- function(model_result, covariate, site_covs, 
                                 n_points = 100, ci = 0.95) {
  
  model <- model_result$model
  species_name <- model_result$species
  
  # Get original (unscaled) covariate values
  original_values <- site_covs[[covariate]]
  
  # Get scaling parameters
  cov_mean <- mean(original_values, na.rm = TRUE)
  cov_sd <- sd(original_values, na.rm = TRUE)
  
  # Create prediction sequence (original scale)
  x_seq <- seq(min(original_values, na.rm = TRUE), 
               max(original_values, na.rm = TRUE), 
               length.out = n_points)
  
  # Scale for prediction
  x_scaled <- (x_seq - cov_mean) / cov_sd
  
  # Get covariate names from model
  cov_names <- colnames(model$beta.samples)
  
  # Find which column corresponds to this covariate
  cov_index <- which(cov_names == covariate)
  
  if (length(cov_index) == 0) {
    stop(paste("Covariate", covariate, "not found in model. Available:", 
               paste(cov_names, collapse = ", ")))
  }
  
  # Get posterior samples
  beta_samples <- model$beta.samples
  n_samples <- nrow(beta_samples)
  
  # Intercept index
  int_index <- which(cov_names == "(Intercept)")
  
  # Predict occupancy for each posterior sample
  # Hold other covariates at mean (0 on scaled)
  pred_matrix <- matrix(NA, nrow = n_samples, ncol = n_points)
  
  for (i in 1:n_samples) {
    linear_pred <- beta_samples[i, int_index] + beta_samples[i, cov_index] * x_scaled
    pred_matrix[i, ] <- plogis(linear_pred)
  }
  
  # Calculate summaries
  pred_mean <- apply(pred_matrix, 2, mean)
  pred_lower <- apply(pred_matrix, 2, quantile, (1 - ci) / 2)
  pred_upper <- apply(pred_matrix, 2, quantile, 1 - (1 - ci) / 2)
  
  # Create prediction dataframe
  pred_df <- data.frame(
    x = x_seq,
    mean = pred_mean,
    lower = pred_lower,
    upper = pred_upper
  )
  
  # Get naive presence/absence per site
  det_matrix <- model_result$metadata$detection_matrix
  if (is.null(det_matrix)) {
    species_filename <- gsub(" ", "_", species_name)
    det_matrix <- readRDS(file.path(input_dir, "detection_matrices", 
                                    paste0("detection_matrix_", species_filename, ".rds")))
  }
  
  naive_occ <- apply(det_matrix, 1, function(x) as.numeric(any(x == 1, na.rm = TRUE)))
  
  obs_df <- data.frame(
    x = original_values,
    y = naive_occ
  )
  
  # Clean axis label
  x_label <- gsub("_", " ", covariate)
  x_label <- gsub("pct", "(%)", x_label)
  x_label <- tools::toTitleCase(x_label)
  
  # Create plot
  p <- ggplot() +
    geom_ribbon(data = pred_df, aes(x = x, ymin = lower, ymax = upper),
                fill = "#5a9e6f", alpha = 0.3) +
    geom_line(data = pred_df, aes(x = x, y = mean),
              color = "#2d6a4f", linewidth = 1.2) +
    geom_point(data = obs_df, aes(x = x, y = y),
               color = "#5a9e6f", alpha = 0.5, size = 2) +
    labs(
      title = species_name,
      x = x_label,
      y = "Occupancy Probability"
    ) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 17),
      panel.grid.minor = element_blank()
    )
  
  return(p)
}

####################
#   Load data
####################

site_covs <- readRDS(file.path(input_dir, "site_covariates.rds"))

covariates <- c("trees_pct", "grass_pct", "shrub_pct", "flooded_veg_pct",
                "water_pct", "crops_pct", "habitat_diversity", "log_area")

model_files <- list.files(file.path(input_dir, "model_results"), 
                          pattern = "^model_.*\\.rds$", full.names = TRUE)

####################
#   Loop through species
####################

for (file in model_files) {
  
  result <- readRDS(file)
  if (!result$success) next
  
  species_name <- result$species
  species_filename <- gsub(" ", "_", species_name)
  
  cat("\nProcessing:", species_name, "\n")
  
  species_plots <- list()
  
  for (cov in covariates) {
    
    cat("  -", cov, "\n")
    
    p <- plot_marginal_effect(result, cov, site_covs)
    species_plots[[cov]] <- p
    
    # Save individual plot
    ggsave(file.path(output_dir, "by_species", paste0(species_filename, "_", cov, ".png")), 
           p, width = 8, height = 6, dpi = 300)
  }
  
  # Combined 7-panel figure
  combined <- wrap_plots(species_plots, ncol = 3) +
    plot_annotation(
      title = species_name,
      theme = theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5))
    )
  
  ggsave(file.path(output_dir, "by_species", paste0(species_filename, "_ALL_COVARIATES.png")),
         combined, width = 15, height = 14, dpi = 300)
  
  cat("  Saved combined figure\n")
}

cat("\n\nDone! Check", file.path(output_dir, "by_species"), "\n")
