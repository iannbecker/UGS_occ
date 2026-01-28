##############################
#
# Within-Site Marginal Effect Figures
# Per species-site combination
#
##############################

library(ggplot2)
library(dplyr)
library(patchwork)

# Set paths
input_dir <- "/Users/ianbecker/Desktop/project_code/UGS_occ/data"
output_dir <- "/Users/ianbecker/Desktop/project_code/UGS_occ/figures_tables/within_site"

# Create output directory if needed
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

####################
#   Marginal effect function (adapted for within-site)
####################

plot_marginal_effect_within <- function(model, covariate, cov_data, 
                                        n_points = 100, ci = 0.95) {
  
  # Get original (unscaled) covariate values for this site
  original_values <- cov_data[[covariate]]
  
  # Get scaling parameters
  cov_mean <- mean(original_values, na.rm = TRUE)
  cov_sd <- sd(original_values, na.rm = TRUE)
  
  # Handle zero variance
  if (is.na(cov_sd) || cov_sd < 0.001) {
    return(NULL)
  }
  
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
    return(NULL)
  }
  
  # Get posterior samples
  beta_samples <- as.matrix(model$beta.samples)
  n_samples <- nrow(beta_samples)
  
  # Intercept index
  int_index <- which(cov_names == "(Intercept)")
  
  # Predict occupancy for each posterior sample
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
  
  # Clean axis label
  x_label <- gsub("_", " ", covariate)
  x_label <- tools::toTitleCase(x_label)
  
  # Check if effect is significant (CI excludes 0 on logit scale)
  beta_cov <- beta_samples[, cov_index]
  sig_pos <- quantile(beta_cov, 0.025) > 0
  sig_neg <- quantile(beta_cov, 0.975) < 0
  
  # Color based on significance
  if (sig_pos) {
    fill_col <- "#2d6a4f"
    line_col <- "#1b4332"
  } else if (sig_neg) {
    fill_col <- "#ae2012"
    line_col <- "#6a040f"
  } else {
    fill_col <- "#6c757d"
    line_col <- "#495057"
  }
  
  # Create plot
  p <- ggplot() +
    geom_ribbon(data = pred_df, aes(x = x, ymin = lower, ymax = upper),
                fill = fill_col, alpha = 0.3) +
    geom_line(data = pred_df, aes(x = x, y = mean),
              color = line_col, linewidth = 1) +
    labs(
      title = x_label,
      x = NULL,
      y = NULL
    ) +
    scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
      axis.text = element_text(size = 8),
      panel.grid.minor = element_blank()
    )
  
  return(p)
}

####################
#   Load within-site results
####################

# Get all within-site model files
model_files <- list.files(file.path(input_dir, "within_site_models"), 
                          pattern = "\\.rds$", 
                          full.names = TRUE)

# Remove the all_results file
model_files <- model_files[!grepl("all_results", model_files)]

cat("Found", length(model_files), "within-site models\n\n")

####################
#   Loop through each species-site model
####################

for (file in model_files) {
  
  result <- tryCatch(readRDS(file), error = function(e) NULL)
  if (is.null(result)) next
  
  species_name <- result$species
  site_id <- result$site_id
  guild <- result$guild
  n_cells <- result$n_cells
  naive_occ <- round(result$naive_occ, 2)
  ppc_pval <- round(result$ppc_pval, 3)
  
  species_filename <- gsub(" ", "_", species_name)
  
  cat("Processing:", species_name, "- Site", site_id, "\n")
  
  # Get covariate names from model
  cov_names <- result$covariate_names
  if (is.null(cov_names)) {
    cov_names <- colnames(result$model$beta.samples)
    cov_names <- cov_names[cov_names != "(Intercept)"]
  }
  
  # Get covariate data (unscaled)
  cov_data <- result$covariates
  
  # If covariates were stored scaled, we need the unscaled version
  # The covariates in result$covariates should be scaled, so we work with that
  # and just use the scaled range for plotting
  
  # Actually let's use the scaled data directly and label axes as "scaled"
  # Or we can back-transform - but for simplicity, use scaled
  
  species_plots <- list()
  
  for (cov in cov_names) {
    
    # Check if covariate exists in data
    if (!cov %in% names(cov_data)) next
    
    p <- tryCatch({
      plot_marginal_effect_within(
        model = result$model,
        covariate = cov,
        cov_data = cov_data
      )
    }, error = function(e) {
      cat("  Error with", cov, ":", e$message, "\n")
      NULL
    })
    
    if (!is.null(p)) {
      species_plots[[cov]] <- p
    }
  }
  
  if (length(species_plots) == 0) {
    cat("  No valid plots - skipping\n")
    next
  }
  
  # Determine layout
  n_plots <- length(species_plots)
  ncol <- min(4, n_plots)
  nrow <- ceiling(n_plots / ncol)
  
  # Combined figure
  combined <- wrap_plots(species_plots, ncol = ncol) +
    plot_annotation(
      title = paste0(species_name, " - Site ", site_id),
      subtitle = paste0("Guild: ", guild, " | Cells: ", n_cells, 
                        " | Naive occ: ", naive_occ, " | PPC p-value: ", ppc_pval),
      theme = theme(
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 10, hjust = 0.5, color = "gray40")
      )
    )
  
  # Add shared axis labels
  combined <- combined + 
    plot_layout(axes = "collect") &
    theme(axis.title = element_text(size = 9))
  
  # Save
  filename <- paste0(species_filename, "_site", site_id, "_effects.png")
  
  ggsave(file.path(output_dir, filename),
         combined, width = 3 * ncol, height = 3 * nrow, dpi = 300)
  
  cat("  Saved:", filename, "\n")
}

cat("\n\nDone! Check", output_dir, "\n")
