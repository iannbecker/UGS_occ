##############################
#
# Simulated NMDS plot for meeting 
# 12/15/2025
# Ian Becker
#
##############################

library(vegan)
library(ggplot2)
library(dplyr)
library(MASS)

set.seed(123)

# Define BCRs in Texas
bcrs <- c("Tamaulipan Brushlands", "Edwards Plateau", "Oaks and Prairies", 
          "Southern Great Plains", "Chihuahuan Desert", "Pineywoods")

# Number of hotspots per BCR (simulate realistic sample sizes)
n_hotspots <- c(15, 12, 20, 8, 6, 10)  # Total = 71 hotspots

##############################
# Scenario 1: STRONG CLUSTERING BY BCR
# (Urban communities retain regional identity)
##############################

simulate_clustered_communities <- function() {
  
  # Create species pool - 200 species total, split into BCR-specific pools
  # Each BCR gets 50 "regional specialists" + 20 shared generalists
  
  communities <- list()
  
  for (i in 1:length(bcrs)) {
    n_sites <- n_hotspots[i]
    
    # Create BCR-specific species matrix
    # Species 1-50 are specific to BCR 1, 51-100 to BCR 2, etc.
    species_matrix <- matrix(0, nrow = n_sites, ncol = 200)
    
    # This BCR's regional specialists (high probability)
    start_col <- (i - 1) * 30 + 1
    end_col <- start_col + 29
    species_matrix[, start_col:end_col] <- matrix(rbinom(n_sites * 30, 1, 0.85), 
                                                  nrow = n_sites)
    
    # Shared urban generalists (columns 181-200, moderate probability)
    species_matrix[, 181:200] <- matrix(rbinom(n_sites * 20, 1, 0.5), 
                                        nrow = n_sites)
    
    # Small chance of other BCRs' species showing up (spillover)
    other_bcr_cols <- setdiff(1:180, start_col:end_col)
    species_matrix[, other_bcr_cols] <- matrix(rbinom(n_sites * length(other_bcr_cols), 
                                                      1, 0.03), nrow = n_sites)
    
    communities[[i]] <- species_matrix
  }
  
  # Combine all BCRs
  full_matrix <- do.call(rbind, communities)
  
  # Remove empty columns
  full_matrix <- full_matrix[, colSums(full_matrix) > 0]
  
  # Create BCR labels
  bcr_labels <- rep(bcrs, times = n_hotspots)
  
  # Run NMDS
  nmds <- metaMDS(full_matrix, distance = "bray", k = 2, trymax = 100)
  
  # Extract scores
  site_scores <- as.data.frame(scores(nmds, display = "sites"))
  site_scores$BCR <- bcr_labels
  
  # Plot
  p1 <- ggplot(site_scores, aes(x = NMDS1, y = NMDS2, color = BCR, shape = BCR)) +
    geom_point(size = 4, alpha = 0.7) +
    stat_ellipse(aes(fill = BCR), alpha = 0.2, geom = "polygon", level = 0.95) +
    scale_color_brewer(palette = "Dark2") +
    scale_fill_brewer(palette = "Dark2") +
    scale_shape_manual(values = c(16, 17, 15, 18, 8, 9)) +
    labs(title = "Scenario 1: STRONG Regional Clustering",
         subtitle = paste0("Stress = ", round(nmds$stress, 3), 
                           " | Urban communities retain regional identity"),
         x = "NMDS1", y = "NMDS2") +
    theme_classic() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 11, color = "gray40"),
      legend.position = "right",
      legend.title = element_text(face = "bold")
    ) +
    coord_fixed()
  
  return(list(plot = p1, stress = nmds$stress, nmds = nmds, data = site_scores))
}

##############################
# Scenario 2: WEAK/NO CLUSTERING (HOMOGENIZATION)
# (Urbanization creates similar communities everywhere)
##############################

simulate_homogenized_communities <- function() {
  
  # Almost all species are shared - minimal regional differences
  communities <- list()
  
  for (i in 1:length(bcrs)) {
    n_sites <- n_hotspots[i]
    
    # Core urban generalists (80 species, very common everywhere)
    # All BCRs have nearly identical species with just random variation
    urban_core <- matrix(rbinom(n_sites * 80, 1, 0.7), nrow = n_sites)
    
    # Tiny regional component (5 species each, low probability)
    regional_tiny <- matrix(rbinom(n_sites * 5, 1, 0.2), nrow = n_sites)
    
    # Some additional shared species (40 species)
    shared_extra <- matrix(rbinom(n_sites * 40, 1, 0.45), nrow = n_sites)
    
    site_matrix <- cbind(urban_core, regional_tiny, shared_extra)
    communities[[i]] <- site_matrix
  }
  
  full_matrix <- do.call(rbind, communities)
  full_matrix <- full_matrix[, colSums(full_matrix) > 0]
  bcr_labels <- rep(bcrs, times = n_hotspots)
  
  nmds <- metaMDS(full_matrix, distance = "bray", k = 2, trymax = 100)
  
  site_scores <- as.data.frame(scores(nmds, display = "sites"))
  site_scores$BCR <- bcr_labels
  
  p2 <- ggplot(site_scores, aes(x = NMDS1, y = NMDS2, color = BCR, shape = BCR)) +
    geom_point(size = 4, alpha = 0.7) +
    stat_ellipse(aes(fill = BCR), alpha = 0.1, geom = "polygon", level = 0.95) +
    scale_color_brewer(palette = "Dark2") +
    scale_fill_brewer(palette = "Dark2") +
    scale_shape_manual(values = c(16, 17, 15, 18, 8, 9)) +
    labs(title = "Scenario 2: WEAK/NO Clustering (Homogenization)",
         subtitle = paste0("Stress = ", round(nmds$stress, 3),
                           " | Urbanization creates similar communities"),
         x = "NMDS1", y = "NMDS2") +
    theme_classic() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 11, color = "gray40"),
      legend.position = "right",
      legend.title = element_text(face = "bold")
    ) +
    coord_fixed()
  
  return(list(plot = p2, stress = nmds$stress, nmds = nmds, data = site_scores))
}

##############################
# Scenario 3: PARTIAL CLUSTERING
# (Some BCRs distinct, others overlap)
##############################

simulate_partial_clustering <- function() {
  
  communities <- list()
  
  for (i in 1:length(bcrs)) {
    n_sites <- n_hotspots[i]
    species_matrix <- matrix(0, nrow = n_sites, ncol = 180)
    
    # STRONG shared urban generalist core (40 species, common everywhere)
    species_matrix[, 141:180] <- matrix(rbinom(n_sites * 40, 1, 0.65), nrow = n_sites)
    
    # Tamaulipan Brushlands (BCR 1) - Moderately distinct
    if (i == 1) {
      # 25 regional specialists (moderate probability)
      species_matrix[, 1:25] <- matrix(rbinom(n_sites * 25, 1, 0.65), nrow = n_sites)
      # Some spillover from adjacent regions (Edwards, Southern GP)
      overlap_cols <- c(26:40, 101:115)
      species_matrix[, overlap_cols] <- matrix(rbinom(n_sites * length(overlap_cols), 1, 0.2), 
                                               nrow = n_sites)
    }
    # Edwards Plateau (BCR 2) - Moderate regional identity
    else if (i == 2) {
      species_matrix[, 26:50] <- matrix(rbinom(n_sites * 25, 1, 0.6), nrow = n_sites)
      # Overlap with neighbors
      overlap_cols <- c(1:15, 51:65, 101:110)
      species_matrix[, overlap_cols] <- matrix(rbinom(n_sites * length(overlap_cols), 1, 0.25), 
                                               nrow = n_sites)
    }
    # Oaks and Prairies (BCR 3) - Overlaps with several regions
    else if (i == 3) {
      species_matrix[, 51:75] <- matrix(rbinom(n_sites * 25, 1, 0.55), nrow = n_sites)
      # High overlap with Edwards and Pineywoods
      overlap_cols <- c(26:45, 76:95)
      species_matrix[, overlap_cols] <- matrix(rbinom(n_sites * length(overlap_cols), 1, 0.3), 
                                               nrow = n_sites)
    }
    # Southern Great Plains (BCR 4) - Similar to Oaks/Prairies
    else if (i == 4) {
      species_matrix[, 101:120] <- matrix(rbinom(n_sites * 20, 1, 0.6), nrow = n_sites)
      # Overlaps with Oaks/Prairies and Edwards
      overlap_cols <- c(26:40, 51:75)
      species_matrix[, overlap_cols] <- matrix(rbinom(n_sites * length(overlap_cols), 1, 0.35), 
                                               nrow = n_sites)
    }
    # Chihuahuan Desert (BCR 5) - Most distinct but not isolated
    else if (i == 5) {
      species_matrix[, 121:140] <- matrix(rbinom(n_sites * 20, 1, 0.7), nrow = n_sites)
      # Some Tamaulipan overlap (adjacent)
      species_matrix[, 1:15] <- matrix(rbinom(n_sites * 15, 1, 0.15), nrow = n_sites)
    }
    # Pineywoods (BCR 6) - Fairly distinct eastern species
    else if (i == 6) {
      species_matrix[, 76:100] <- matrix(rbinom(n_sites * 25, 1, 0.6), nrow = n_sites)
      # Some overlap with Oaks/Prairies
      species_matrix[, 51:70] <- matrix(rbinom(n_sites * 20, 1, 0.25), nrow = n_sites)
    }
    
    communities[[i]] <- species_matrix
  }
  
  full_matrix <- do.call(rbind, communities)
  full_matrix <- full_matrix[, colSums(full_matrix) > 0]
  bcr_labels <- rep(bcrs, times = n_hotspots)
  
  nmds <- metaMDS(full_matrix, distance = "bray", k = 2, trymax = 100)
  
  site_scores <- as.data.frame(scores(nmds, display = "sites"))
  site_scores$BCR <- bcr_labels
  
  p3 <- ggplot(site_scores, aes(x = NMDS1, y = NMDS2, color = BCR, shape = BCR)) +
    geom_point(size = 4, alpha = 0.7) +
    stat_ellipse(aes(fill = BCR), alpha = 0.2, geom = "polygon", level = 0.95) +
    scale_color_brewer(palette = "Dark2") +
    scale_fill_brewer(palette = "Dark2") +
    scale_shape_manual(values = c(16, 17, 15, 18, 8, 9)) +
    labs(title = "Scenario 3: PARTIAL Clustering",
         subtitle = paste0("Stress = ", round(nmds$stress, 3),
                           " | Some BCRs distinct, others homogenized"),
         x = "NMDS1", y = "NMDS2") +
    theme_classic() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 11, color = "gray40"),
      legend.position = "right",
      legend.title = element_text(face = "bold")
    ) +
    coord_fixed()
  
  return(list(plot = p3, stress = nmds$stress, nmds = nmds, data = site_scores))
}

##############################
# Generate all scenarios
##############################

cat("Generating Scenario 1: Strong Clustering...\n")
scenario1 <- simulate_clustered_communities()

cat("Generating Scenario 2: Homogenization...\n")
scenario2 <- simulate_homogenized_communities()

cat("Generating Scenario 3: Partial Clustering...\n")
scenario3 <- simulate_partial_clustering()

##############################
# Display plots
##############################

print(scenario1$plot)
print(scenario2$plot)
print(scenario3$plot)
