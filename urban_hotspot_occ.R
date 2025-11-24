##############################
#
# Occupancy Modelling 
# 8/25/2025
# Ian Becker 
#
##############################

library(spOccupancy)
library(dplyr)
library(tidyr)
library(forcats)
library(ggplot2)
library(terra)
library(sf)
library(parallel)

setwd("C:/Users/ianbe/OneDrive - The University of Texas-Rio Grande Valley/CampusBirds/ugs")

####################
#   Data Prep
####################

# Loading in species data - now handles any number of species

species_data <- readRDS("urban_hotspot_species_inat.rds")

# Function to convert species data to sf objects

convert_to_sf <- function(species_df, species_name) {
  if (nrow(species_df) == 0) {
    cat("Warning: No data for", species_name, "\n")
    return(NULL)
  }
  
  # Check for coordinate columns (handle different naming conventions)
  coord_cols <- c("longitude", "latitude", "decimalLongitude", "decimalLatitude", 
                  "lon", "lat", "x", "y")
  
  found_coords <- coord_cols[coord_cols %in% names(species_df)]
  if (length(found_coords) < 2) {
    cat("Warning: Cannot find coordinate columns for", species_name, "\n")
    return(NULL)
  }
  
  # Use first two coordinate columns found
  lon_col <- found_coords[1]
  lat_col <- found_coords[2]
  
  # Remove rows with missing coordinates
  species_df <- species_df[!is.na(species_df[[lon_col]]) & !is.na(species_df[[lat_col]]), ]
  
  if (nrow(species_df) == 0) {
    cat("Warning: No valid coordinates for", species_name, "\n")
    return(NULL)
  }
  
  return(st_as_sf(species_df, coords = c(lon_col, lat_col), crs = 4326))
}

# Convert all species to sf objects

species_sf_list <- list()
species_names <- names(species_data)

if (is.null(species_names)) {
  # If no names, create them
  species_names <- paste0("Species_", sprintf("%03d", 1:length(species_data)))
  names(species_data) <- species_names
}

cat("Converting", length(species_data), "species to SF objects...\n")

for (i in seq_along(species_data)) {
  species_name <- species_names[i]
  cat("Processing", species_name, "...\n")
  
  sf_obj <- convert_to_sf(species_data[[i]], species_name)
  if (!is.null(sf_obj)) {
    species_sf_list[[species_name]] <- sf_obj
  }
}

cat("Successfully converted", length(species_sf_list), "out of", length(species_data), "species\n")

# Loading in landcover data and shapefile

landcover <- rast("urban_hotspots_landcover.tif")
urban_hotspot <- st_read("urban_hotspot_shapefile")

# Matching crs

urban_hotspot <- st_transform(urban_hotspot, crs = 32614)
target_crs <- crs(urban_hotspot)  

# Transform all species data to target CRS

for (species_name in names(species_sf_list)) {
  if (st_crs(species_sf_list[[species_name]])$input != target_crs) {
    species_sf_list[[species_name]] <- st_transform(species_sf_list[[species_name]], target_crs)
  }
}

landcover <- project(landcover, target_crs)
cat("Reprojected all data to UTM Zone 14N (meters)\n")

####################
#   Study Grid (unchanged)
####################

cell_size <- 150 # 150 m x 150 m
study_grid <- st_make_grid(urban_hotspot, cellsize = cell_size, square = FALSE)

grid_intersects <- st_intersects(study_grid, urban_hotspot, sparse = FALSE)[,1]
study_grid <- study_grid[grid_intersects]
study_grid <- st_sf(site_id = 1:length(study_grid), geometry = study_grid)

cat("Created", nrow(study_grid), "150m grid cells\n")
saveRDS(study_grid, "study_grid.rds")

####################
#   Landcover Data (unchanged but could be optimized)
####################

# Results dataframe

landcover_vars <- data.frame(
  site_id = integer(),
  forest_pct = numeric(), 
  grass_pct = numeric(), 
  crops_pct = numeric(),
  urban_pct = numeric(),
  # water_pct = numeric(),        # REMOVE THIS LINE
  bare_pct = numeric(),
  shrub_pct = numeric(),
  flooded_veg_pct = numeric(),
  dominant_class = numeric(),
  habitat_diversity = numeric(),
  # edge_density = numeric(),     # REMOVE THIS LINE
  stringsAsFactors = FALSE
)

# Extract data for each cell (consider parallelizing this for very large grids)

for (i in 1:nrow(study_grid)) {
  
  if (i %% 100 == 0) cat("Processing cell", i, "of", nrow(study_grid), "\n")
  
  grid_cell <- study_grid[i,]
  lc_values <- terra::extract(landcover, vect(grid_cell), df = TRUE)
  
  if (nrow(lc_values) == 0) {
    # No data for this cell - REMOVE water_pct and edge_density
    cell_results <- data.frame(
      site_id = grid_cell$site_id,
      forest_pct = 0, grass_pct = 0, crops_pct = 0, urban_pct = 0, 
      # water_pct = 0,              # REMOVE THIS LINE
      bare_pct = 0,
      shrub_pct = 0, flooded_veg_pct = 0,
      dominant_class = NA,
      habitat_diversity = 0
      # edge_density = 0            # REMOVE THIS LINE
    )
  } else {
    
    # Calculate land cover percentages
    lc_column <- names(lc_values)[2]
    total_pixels <- nrow(lc_values)
    
    lc_summary <- lc_values %>%
      count(.data[[lc_column]]) %>%
      mutate(percentage = n / total_pixels * 100)
    
    names(lc_summary)[1] <- "class"
    
    # Initialize variables - REMOVE water_pct calculation
    forest_pct <- ifelse(1 %in% lc_summary$class, 
                         lc_summary$percentage[lc_summary$class == 1], 0)
    grass_pct <- ifelse(2 %in% lc_summary$class, 
                        lc_summary$percentage[lc_summary$class == 2], 0)
    flooded_veg_pct <- ifelse(3 %in% lc_summary$class,
                              lc_summary$percentage[lc_summary$class == 3], 0)
    crops_pct <- ifelse(4 %in% lc_summary$class, 
                        lc_summary$percentage[lc_summary$class == 4], 0)
    shrub_pct <- ifelse(5 %in% lc_summary$class,
                        lc_summary$percentage[lc_summary$class == 5], 0)
    urban_pct <- ifelse(6 %in% lc_summary$class, 
                        lc_summary$percentage[lc_summary$class == 6], 0)
    bare_pct <- ifelse(7 %in% lc_summary$class, 
                       lc_summary$percentage[lc_summary$class == 7], 0)
    # water_pct <- ifelse(0 %in% lc_summary$class,    # REMOVE THIS BLOCK
    #                     lc_summary$percentage[lc_summary$class == 0], 0)
    
    dominant_class <- lc_summary$class[which.max(lc_summary$percentage)]
    
    # Habitat diversity (Shannon index)
    proportions <- lc_summary$percentage / 100
    proportions <- proportions[proportions > 0]
    habitat_diversity <- -sum(proportions * log(proportions), na.rm = TRUE)
    
    # edge_density <- nrow(lc_summary)  # REMOVE THIS LINE
    
    # Complete cell_results - REMOVE water_pct and edge_density
    cell_results <- data.frame(
      site_id = grid_cell$site_id,
      forest_pct = forest_pct,
      grass_pct = grass_pct, 
      crops_pct = crops_pct,
      urban_pct = urban_pct,
      # water_pct = water_pct,          # REMOVE THIS LINE
      bare_pct = bare_pct,
      shrub_pct = shrub_pct,
      flooded_veg_pct = flooded_veg_pct,
      dominant_class = dominant_class,
      habitat_diversity = habitat_diversity
      # edge_density = edge_density     # REMOVE THIS LINE
    )
  }
  
  landcover_vars <- rbind(landcover_vars, cell_results)
}

####################
#   Creating Detection Matrices for All Species
####################

cat("Creating detection histories for", length(species_sf_list), "species...\n")

detection_histories <- list()

for (species_code in names(species_sf_list)) {
  
  cat("Creating detection history for", species_code, "...\n")
  
  current_species_data <- species_sf_list[[species_code]]
  
  # Check if data exists
  if (is.null(current_species_data) || nrow(current_species_data) == 0) {
    cat("No observations found for", species_code, "\n")
    next
  }
  
  # Assign observations to grid cells
  obs_intersect <- st_intersects(current_species_data, study_grid)
  current_species_data$site_id <- sapply(obs_intersect, function(x) ifelse(length(x) > 0, x[1], NA))
  
  # Remove observations not in any grid cell
  current_species_data <- current_species_data[!is.na(current_species_data$site_id),]
  
  if (nrow(current_species_data) == 0) {
    cat("No observations intersect with grid for", species_code, "\n")
    next
  }
  
  # Find date column (flexible naming)
  date_cols <- c("observed_on", "eventDate", "date", "Date", "observation_date")
  date_col <- date_cols[date_cols %in% colnames(current_species_data)][1]
  
  if (is.na(date_col)) {
    cat("No date column found for", species_code, "\n")
    next
  }
  
  # Extract year
  current_species_data$year <- as.numeric(format(as.Date(current_species_data[[date_col]]), "%Y"))
  current_species_data <- current_species_data[!is.na(current_species_data$year), ]
  
  # Define sampling occasions (years)
  sampling_occasions <- sort(unique(current_species_data$year))
  n_occasions <- length(sampling_occasions)
  
  if (n_occasions < 2) {
    cat("Need multiple sampling occasions for", species_code, "(only", n_occasions, "year found)\n")
    next
  }
  
  # Create detection matrix
  n_sites <- nrow(study_grid)
  detection_matrix <- matrix(0, nrow = n_sites, ncol = n_occasions)
  
  # Fill in detections
  for (i in 1:nrow(current_species_data)) {
    site_idx <- current_species_data$site_id[i]
    occasion_idx <- which(sampling_occasions == current_species_data$year[i])
    if (length(occasion_idx) > 0) {
      detection_matrix[site_idx, occasion_idx] <- 1
    }
  }
  
  # Store detection history
  detection_histories[[species_code]] <- list(
    y = detection_matrix,
    species = species_code,
    n_sites = n_sites,
    n_occasions = n_occasions,
    occasions = sampling_occasions,
    obs_count = nrow(current_species_data)
  )
  
  cat("✓", species_code, "- Sites:", n_sites, 
      "| Occasions:", n_occasions, 
      "| Observations:", nrow(current_species_data), "\n")
}

cat("\nSuccessfully created detection histories for", length(detection_histories), "species\n")

####################
#   Batch Model Fitting with Progress Tracking
####################

# Confirming valid species (25 observations)

min_observations <- 25
valid_species <- names(detection_histories)[sapply(detection_histories, function(x) x$obs_count >= min_observations)]

cat("Fitting models for", length(valid_species), "species with >=", min_observations, "observations\n")

# Initialize results storage

occupancy_results <- list()
failed_models <- character()
successful_models <- character()

# Prepare covariates once (same for all species)

occ_covs <- landcover_vars[, c("forest_pct", "grass_pct", "urban_pct", 
                               "crops_pct", "shrub_pct", "flooded_veg_pct", 
                               "habitat_diversity")]
occ_covs_scaled <- scale(occ_covs)

# Fit models with progress tracking

start_time <- Sys.time()

for (i in seq_along(valid_species)) {
  species_code <- valid_species[i]
  det_hist <- detection_histories[[species_code]]
  
  cat(paste0("[", i, "/", length(valid_species), "] Fitting model for ", species_code, "...\n"))
  
  # Detection covariates (intercept only)
  det_covs <- list(
    intercept = matrix(1, nrow = det_hist$n_sites, ncol = det_hist$n_occasions)
  )
  
  # Fit model with error handling
  tryCatch({
    model <- PGOcc(
      occ.formula = ~ forest_pct + grass_pct + urban_pct + crops_pct + 
        shrub_pct + flooded_veg_pct + habitat_diversity,
      det.formula = ~ 1,
      data = list(y = det_hist$y, 
                  occ.covs = data.frame(occ_covs_scaled),
                  det.covs = det_covs),
      n.samples = 8000,
      n.burn = 4000,
      n.thin = 2,
      n.chains = 2,
      n.omp.threads = 1,
      verbose = FALSE,
      n.report = 800
    )
    
    # Store successful result
    occupancy_results[[species_code]] <- list(
      model = model,
      detection_history = det_hist,
      summary = summary(model),
      mean_occupancy = round(mean(model$psi.samples), 3),
      success = TRUE,
      error_message = NULL
    )
    
    successful_models <- c(successful_models, species_code)
    cat("✓ Mean occupancy:", round(mean(model$psi.samples), 3), "\n")
    
  }, error = function(e) {
    
    # Store failed result
    occupancy_results[[species_code]] <- list(
      model = NULL,
      detection_history = det_hist,
      summary = NULL,
      mean_occupancy = NA,
      success = FALSE,
      error_message = e$message
    )
    
    failed_models <- c(failed_models, species_code)
    cat("✗ Failed:", e$message, "\n")
  })
  
  # Progress update
  if (i %% 5 == 0 || i == length(valid_species)) {
    elapsed <- difftime(Sys.time(), start_time, units = "mins")
    cat("Progress:", round(i/length(valid_species)*100, 1), "% | Elapsed:", 
        round(elapsed, 1), "min | Success:", length(successful_models), 
        "| Failed:", length(failed_models), "\n")
  }
}

####################
#   Results Summary
####################

cat("\n=== FINAL RESULTS ===\n")
cat("Total species processed:", length(valid_species), "\n")
cat("Successful models:", length(successful_models), "\n")
cat("Failed models:", length(failed_models), "\n")

if (length(failed_models) > 0) {
  cat("\nFailed species:", paste(failed_models, collapse = ", "), "\n")
}

# Summary statistics for successful models

if (length(successful_models) > 0) {
  mean_occupancies <- sapply(successful_models, function(x) occupancy_results[[x]]$mean_occupancy)
  
  cat("\nOccupancy Summary:\n")
  cat("Mean:", round(mean(mean_occupancies, na.rm = TRUE), 3), "\n")
  cat("Median:", round(median(mean_occupancies, na.rm = TRUE), 3), "\n")
  cat("Range:", round(min(mean_occupancies, na.rm = TRUE), 3), "-", 
      round(max(mean_occupancies, na.rm = TRUE), 3), "\n")
}

# Save results

save(occupancy_results, detection_histories, landcover_vars, 
     file = paste0("multi_species_occupancy_results_", Sys.Date(), ".RData"))

cat("\nResults saved to multi_species_occupancy_results_", Sys.Date(), ".RData\n")

####################
#   Forest Plots
####################

# Function to extract coefficient summaries from spOccupancy model

extract_coef_summary <- function(model_result, species_name) {
  
  if (!model_result$success || is.null(model_result$model)) {
    cat("Skipping", species_name, "- model failed or is NULL\n")
    return(NULL)
  }
  
  model <- model_result$model
  
  # Try to get summary first, if that fails, create it manually
  model_summary <- NULL
  
  # Try to get existing summary
  if (!is.null(model_result$summary)) {
    model_summary <- model_result$summary
  } else {
    # Try to create summary
    tryCatch({
      model_summary <- summary(model)
    }, error = function(e) {
      cat("Could not create summary for", species_name, ":", e$message, "\n")
    })
  }
  
  # If summary approach fails, extract directly from model object
  if (is.null(model_summary) || is.null(model_summary$beta)) {
    
    cat("Extracting coefficients directly from model for", species_name, "...\n")
    
    # Extract beta samples directly
    if (is.null(model$beta.samples)) {
      cat("Skipping", species_name, "- no beta samples found in model\n")
      return(NULL)
    }
    
    beta_samples <- model$beta.samples
    
    # Get coefficient names from model formula or data
    if (!is.null(model$X)) {
      covariate_names <- colnames(model$X)
    } else if (!is.null(colnames(beta_samples))) {
      covariate_names <- colnames(beta_samples)
    } else {
      # Default names if nothing else available
      covariate_names <- paste0("coef_", 1:ncol(beta_samples))
    }
    
    # Calculate summary statistics manually
    beta_mean <- apply(beta_samples, 2, mean, na.rm = TRUE)
    beta_sd <- apply(beta_samples, 2, sd, na.rm = TRUE)
    beta_q025 <- apply(beta_samples, 2, quantile, 0.025, na.rm = TRUE)
    beta_q975 <- apply(beta_samples, 2, quantile, 0.975, na.rm = TRUE)
    
    # Create data frame
    coef_df <- data.frame(
      species = rep(species_name, length(beta_mean)),
      covariate = covariate_names,
      mean = beta_mean,
      sd = beta_sd,
      q025 = beta_q025,
      q975 = beta_q975,
      stringsAsFactors = FALSE
    )
    
  } else {
    # Use summary approach
    beta_summary <- model_summary$beta
    
    # Check if beta summary has rows
    if (is.null(beta_summary) || nrow(beta_summary) == 0) {
      cat("Skipping", species_name, "- empty beta summary\n")
      return(NULL)
    }
    
    # Get covariate names
    covariate_names <- rownames(beta_summary)
    if (is.null(covariate_names) || length(covariate_names) == 0) {
      cat("Skipping", species_name, "- no covariate names in summary\n")
      return(NULL)
    }
    
    # Check required columns
    required_cols <- c("Mean", "SD", "2.5%", "97.5%")
    if (!all(required_cols %in% colnames(beta_summary))) {
      cat("Available columns for", species_name, ":", paste(colnames(beta_summary), collapse = ", "), "\n")
      return(NULL)
    }
    
    # Create data frame
    coef_df <- data.frame(
      species = rep(species_name, nrow(beta_summary)),
      covariate = covariate_names,
      mean = beta_summary[, "Mean"],
      sd = beta_summary[, "SD"],
      q025 = beta_summary[, "2.5%"],
      q975 = beta_summary[, "97.5%"],
      stringsAsFactors = FALSE
    )
  }
  
  # Add significance indicator (95% CI doesn't include 0)
  coef_df$significant <- !(coef_df$q025 <= 0 & coef_df$q975 >= 0)
  
  cat("✓ Extracted", nrow(coef_df), "coefficients for", species_name, "\n")
  return(coef_df)
}

# Extract summaries for all successful models

cat("Extracting coefficient summaries...\n")

all_coef_summaries <- list()

for (species_name in names(occupancy_results)) {
  coef_summary <- extract_coef_summary(occupancy_results[[species_name]], species_name)
  if (!is.null(coef_summary)) {
    all_coef_summaries[[species_name]] <- coef_summary
  }
}

# Combine into single data frame

if (length(all_coef_summaries) > 0) {
  coef_data <- do.call(rbind, all_coef_summaries)
  rownames(coef_data) <- NULL
  
  cat("Extracted coefficients for", length(unique(coef_data$species)), "species\n")
} else {
  stop("No successful models found to extract coefficients from!")
}



###################   Create Forest Plots for Each Covariate

# Get unique covariates (excluding intercept)

covariates <- unique(coef_data$covariate)
covariates <- covariates[covariates != "(Intercept)"]

cat("Creating forest plots for covariates:", paste(covariates, collapse = ", "), "\n")

# Create plots for each covariate using for loop

forest_plots <- list()

for (cov in covariates) {
  
  cat("Creating forest plot for", cov, "...\n")
  
  # Filter for specific covariate
  plot_data <- coef_data %>%
    filter(covariate == cov) %>%
    arrange(mean)
  
  if (nrow(plot_data) == 0) {
    cat("No data found for covariate:", cov, "\n")
    next
  }
  
  # Reorder species by effect size
  plot_data$species <- factor(plot_data$species, levels = plot_data$species)
  
  # Create the plot
  p <- ggplot(plot_data, aes(x = mean, y = species)) +
    
    # Add vertical line at 0 (no effect)
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50", alpha = 0.8) +
    
    # Add confidence intervals
    geom_errorbarh(aes(xmin = q025, xmax = q975, color = significant), 
                   height = 0.3, size = 0.8, alpha = 0.7) +
    
    # Add point estimates
    geom_point(aes(color = significant), size = 2.5, alpha = 0.8) +
    
    # Customize colors
    scale_color_manual(values = c("TRUE" = "darkred", "FALSE" = "gray60"),
                       name = "Significant\n(95% CI)",
                       labels = c("FALSE" = "No", "TRUE" = "Yes")) +
    
    # Labels and theme
    labs(
      title = paste("Effect of", gsub("_", " ", tools::toTitleCase(cov)), "on Occupancy"),
      subtitle = paste("Posterior means and 95% credible intervals for", nrow(plot_data), "species"),
      x = "Coefficient Estimate (log-odds scale)",
      y = "Species"
    ) +
    
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 12, color = "gray60"),
      axis.text.y = element_text(size = 10),
      axis.text.x = element_text(size = 10),
      axis.title = element_text(size = 11, face = "bold"),
      legend.position = "bottom",
      panel.grid.minor = element_blank(),
    )
  
  # Store the plot
  forest_plots[[cov]] <- p
  
  # Display each plot as it's created
  print(p)
}

####################
#   Summary Statistics
####################

# Summary of effects by covariate
cat("\n=== SUMMARY OF COVARIATE EFFECTS ===\n")

for (cov in covariates) {
  cov_data <- coef_data %>% filter(covariate == cov)
  
  n_species <- nrow(cov_data)
  n_positive <- sum(cov_data$mean > 0)
  n_negative <- sum(cov_data$mean < 0)
  n_significant <- sum(cov_data$significant)
  n_pos_sig <- sum(cov_data$mean > 0 & cov_data$significant)
  n_neg_sig <- sum(cov_data$mean < 0 & cov_data$significant)
  
  cat("\n", toupper(gsub("_", " ", cov)), ":\n")
  cat("  Total species:", n_species, "\n")
  cat("  Positive effects:", n_positive, "(", round(n_positive/n_species*100, 1), "%)\n")
  cat("  Negative effects:", n_negative, "(", round(n_negative/n_species*100, 1), "%)\n")
  cat("  Significant effects:", n_significant, "(", round(n_significant/n_species*100, 1), "%)\n")
  cat("  Significant positive:", n_pos_sig, "\n")
  cat("  Significant negative:", n_neg_sig, "\n")
  
  if (n_significant > 0) {
    cat("  Mean effect size (significant only):", round(mean(cov_data$mean[cov_data$significant]), 3), "\n")
  }
}

####################
#   Detailed Species Lists for Each Covariate
####################

cat("\n=== SPECIES WITH STRONGEST ASSOCIATIONS FOR EACH COVARIATE ===\n")

for (cov in covariates) {
  
  cat("\n", rep("=", 60), "\n")
  cat("COVARIATE:", toupper(gsub("_", " ", cov)), "\n")
  cat(rep("=", 60), "\n")
  
  # Filter data for this covariate
  cov_data <- coef_data %>% 
    filter(covariate == cov) %>%
    arrange(desc(mean))
  
  # Top positive associations
  cat("\nTOP 10 POSITIVE ASSOCIATIONS:\n")
  top_positive <- head(cov_data, 10)
  for (i in 1:nrow(top_positive)) {
    sig_indicator <- ifelse(top_positive$significant[i], "***", "")
    cat(sprintf("%2d. %-20s: %6.3f [%6.3f, %6.3f] %s\n", 
                i, top_positive$species[i], top_positive$mean[i], 
                top_positive$q025[i], top_positive$q975[i], sig_indicator))
  }
  
  # Top negative associations
  cat("\nTOP 10 NEGATIVE ASSOCIATIONS:\n")
  top_negative <- tail(cov_data, 10)
  for (i in 1:nrow(top_negative)) {
    sig_indicator <- ifelse(top_negative$significant[i], "***", "")
    cat(sprintf("%2d. %-20s: %6.3f [%6.3f, %6.3f] %s\n", 
                i, top_negative$species[i], top_negative$mean[i], 
                top_negative$q025[i], top_negative$q975[i], sig_indicator))
  }
}

####################
#   Create Multi-Panel Plot
####################

forest_plot_multiple <- if (length(forest_plots) > 1) {
  
  cat("\nCreating multi-panel forest plot...\n")
  
  # Load gridExtra for multiple plots
  if (!require(gridExtra, quietly = TRUE)) {
    install.packages("gridExtra")
    library(gridExtra)
  }
  
  # Create multi-panel plot (adjust based on number of covariates)
  if (length(forest_plots) <= 4) {
    multi_plot <- do.call(grid.arrange, c(forest_plots, ncol = 2))
  } else {
    # Show first 4 plots if more than 4 covariates
    multi_plot <- do.call(grid.arrange, c(forest_plots[1:4], ncol = 2))
  }
}

ggsave(filename = "forest_plot_multiple.png", 
       plot = forest_plot_multiple, 
       width = 15, height = 25,
       dpi = 300)

####################
#   Save Results
####################

# Save coefficient data

write.csv(coef_data, "occupancy_coefficients_summary.csv", row.names = FALSE)

# Save forest plots

for (cov_name in names(forest_plots)) {
  if (!is.null(forest_plots[[cov_name]])) {
    ggsave(filename = paste0("forest_plot_", cov_name, ".png"), 
           plot = forest_plots[[cov_name]], 
           width = 10, height = max(6, nrow(coef_data[coef_data$covariate == cov_name,]) * 0.3),
           dpi = 300)
  }
}

cat("\nResults saved:\n")
cat("- Coefficient summary: occupancy_coefficients_summary.csv\n")
cat("- Forest plots: forest_plot_[covariate].png\n")

# Return the forest plots list for further use
forest_plots$urban_pct

cor_matrix <- cor(landcover_vars[, c("forest_pct", "grass_pct", "urban_pct", 
                                     "crops_pct", "shrub_pct", 
                                     "flooded_veg_pct", "habitat_diversity", 
                                     )], use = "complete.obs")
print(round(cor_matrix, 2))
