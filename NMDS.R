##############################
#
# Urban vs Non-Urban Hotspot NMDS Analysis
# 7/2/2025
# Ian Becker
#
##############################

library(auk)
library(dplyr)
library(terra)
library(sf)
library(vegan)
library(tigris)
library(tidyr)
library(tibble)
library(lubridate)
library(ggplot2)

options(tigris_use_cache = TRUE)

# Load two datasets

urban_hotspots <- readRDS("final_LRGV_hotspots.rds")     
nonurban_hotspots <- readRDS("final_LRGV_hotspots_NOTURBAN.rds") 

####################### Data Preparation for NMDS

# Function to clean and prepare hotspot data

prepare_hotspot_data <- function(data, site_type) {
  data_clean <- data %>%
    filter(
      # Effort-based filtering
      duration_minutes >= 5,           
      duration_minutes <= 300,         
      number_observers <= 10,          
      protocol_type %in% c("Stationary", "Traveling"),
      
      # Temporal filtering - non-breeding season
     month(observation_date) %in% 11:3  # November through March
    ) %>%
    # Handle 'X' counts 
    mutate(observation_count = ifelse(observation_count == "X", 1, 
                                      as.numeric(observation_count))) %>%
    # Add site type identifier
    mutate(site_type = site_type)
  
  return(data_clean)
}

# Clean both datasets

urban_clean <- prepare_hotspot_data(urban_hotspots, "Urban")
nonurban_clean <- prepare_hotspot_data(nonurban_hotspots, "Non-Urban")

# Combine datasets

combined_data <- bind_rows(urban_clean, nonurban_clean)

# Create species matrix with hotspots as sites

species_matrix_freq <- combined_data %>%
  # Use locality_id as the site identifier (hotspot ID)
  group_by(locality_id, scientific_name, site_type) %>%
  summarise(
    n_checklists_with_species = n_distinct(sampling_event_identifier),
    .groups = "drop"
  ) %>%
  # Get total checklists per hotspot
  left_join(
    combined_data %>%
      group_by(locality_id) %>%
      summarise(
        total_checklists = n_distinct(sampling_event_identifier),
        site_type = first(site_type)  # Preserve site type
      ),
    by = "locality_id"
  ) %>%
  # Calculate frequency
  mutate(frequency = n_checklists_with_species / total_checklists) %>%
  select(locality_id, scientific_name, frequency, site_type.x) %>%
  rename(site_type = site_type.x) %>%
  # Convert to wide format
  tidyr::pivot_wider(names_from = scientific_name, 
                     values_from = frequency, 
                     values_fill = 0)

# Prepare data for NMDS
species_matrix <- species_matrix_freq %>%
  filter(!is.na(locality_id)) %>%
  column_to_rownames("locality_id") %>%
  select(-site_type) %>%  # Remove site_type for matrix, but keep for later
  as.matrix()

# Create site type vector for plotting
site_types <- species_matrix_freq %>%
  filter(!is.na(locality_id)) %>%
  select(locality_id, site_type) %>%
  column_to_rownames("locality_id")

# Filter sites with minimum species richness
min_species <- 5  # Minimum number of species per hotspot
species_matrix <- species_matrix[rowSums(species_matrix > 0) >= min_species, ]

# Filter species by frequency across sites
species_freq <- colSums(species_matrix > 0)
min_sites <- 3  # Species must occur in at least 3 hotspots
species_matrix <- species_matrix[, species_freq >= min_sites]

# Convert to presence/absence if desired (optional)
# species_matrix <- (species_matrix > 0) * 1 

# Remove empty rows
species_matrix <- species_matrix[rowSums(species_matrix) > 0, ]

# Match site types to final matrix
site_types_matched <- site_types[rownames(species_matrix), , drop = FALSE]

print(paste("Final matrix dimensions:", nrow(species_matrix), "sites x", ncol(species_matrix), "species"))
print(paste("Urban sites:", sum(site_types_matched$site_type == "Urban")))
print(paste("Non-urban sites:", sum(site_types_matched$site_type == "Non-Urban")))

####################### NMDS Analysis

# Run NMDS
set.seed(123) 
nmds_result <- metaMDS(species_matrix, 
                       distance = "bray",
                       k = 2,  # 2 dimensions
                       try = 100,
                       trymax = 500)

# Check stress
print(paste("NMDS Stress:", round(nmds_result$stress, 3)))

# Stress interpretation:
# < 0.05: excellent
# 0.05-0.1: good
# 0.1-0.2: fair
# > 0.2: poor

####################### Statistical Tests

# PERMANOVA to test for differences between urban and non-urban
permanova_result <- adonis2(species_matrix ~ site_type, 
                            data = site_types_matched, 
                            permutations = 999,
                            method = "bray")

print("PERMANOVA Results:")
print(permanova_result)

# ANOSIM (alternative test)
anosim_result <- anosim(species_matrix, 
                        site_types_matched$site_type, 
                        permutations = 999)

print("ANOSIM Results:")
print(anosim_result)

# Dispersion test (PERMDISP)
disp_result <- betadisper(vegdist(species_matrix, method = "bray"), 
                          site_types_matched$site_type)
permdisp_test <- permutest(disp_result, permutations = 999)

print("PERMDISP Results:")
print(permdisp_test)

####################### Species Indicators

# Find indicator species for each site type
library(indicspecies)

# Convert site types to factor
site_groups <- as.factor(site_types_matched$site_type)

# Run indicator species analysis
indic_result <- multipatt(species_matrix, site_groups, 
                          control = how(nperm = 999))

print("Indicator Species Analysis:")
summary(indic_result, indvalcomp = TRUE)

####################### Plotting

# Extract NMDS scores
site_scores <- as.data.frame(scores(nmds_result, display = "sites"))
species_scores <- as.data.frame(scores(nmds_result, display = "species"))

# Add site type information
site_scores$site_type <- site_types_matched$site_type
site_scores$hotspot_id <- rownames(site_scores)

# Add species names
species_scores$species <- rownames(species_scores)

# Create NMDS plot
p1 <- ggplot() +
  # Add sites
  geom_point(data = site_scores, 
             aes(x = NMDS1, y = NMDS2, color = site_type, shape = site_type),
             size = 3, alpha = 0.7) +
  
  # Add species (optional - only show most influential)
  geom_text(data = species_scores[1:20, ],  # Top 20 species by position
            aes(x = NMDS1, y = NMDS2, label = species),
            size = 2, alpha = 0.6, color = "gray30") +
  
  # Formatting
  scale_color_manual(values = c("Urban" = "red", "Non-Urban" = "blue"),
                     name = "Site Type") +
  scale_shape_manual(values = c("Urban" = 16, "Non-Urban" = 17),
                     name = "Site Type") +
  
  labs(title = "NMDS Ordination: Urban vs Non-Urban Hotspots",
       subtitle = paste("Stress =", round(nmds_result$stress, 3)),
       x = "NMDS1",
       y = "NMDS2") +
  
  theme_classic() +
  theme(
    plot.title = element_text(size = 14, hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(size = 12, hjust = 0.5),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9)
  ) +
  coord_fixed()

print(p1)

# Create dispersion plot
p2 <- ggplot(data.frame(Distance = disp_result$distances,
                        Group = disp_result$group),
             aes(x = Group, y = Distance, fill = Group)) +
  geom_boxplot(alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  scale_fill_manual(values = c("Urban" = "red", "Non-Urban" = "blue")) +
  labs(title = "Beta Diversity Dispersion",
       subtitle = paste("PERMDISP p-value =", round(permdisp_test$tab$`Pr(>F)`[1], 3)),
       x = "Site Type",
       y = "Distance to Centroid") +
  theme_classic() +
  theme(legend.position = "none")

print(p2)

####################### Candidate Species Selection for Patch Occupancy

# 1. Identify species present in both urban and non-urban hotspots
species_occurrence <- combined_data %>%
  filter(locality_id %in% rownames(species_matrix)) %>%
  group_by(scientific_name, site_type) %>%
  summarise(
    n_hotspots = n_distinct(locality_id),
    n_observations = n(),
    .groups = "drop"
  ) %>%
  pivot_wider(names_from = site_type, 
              values_from = c(n_hotspots, n_observations),
              values_fill = 0) %>%
  # Calculate crossover metrics
  mutate(
    total_urban_hotspots = sum(site_types_matched$site_type == "Urban"),
    total_nonurban_hotspots = sum(site_types_matched$site_type == "Non-Urban"),
    urban_occupancy = `n_hotspots_Urban` / total_urban_hotspots,
    nonurban_occupancy = `n_hotspots_Non-Urban` / total_nonurban_hotspots,
    total_occupancy = (`n_hotspots_Urban` + `n_hotspots_Non-Urban`) / nrow(species_matrix),
    # Key metric: minimum occupancy across both site types
    min_occupancy = pmin(urban_occupancy, nonurban_occupancy),
    # Presence in both site types
    in_both_types = (`n_hotspots_Urban` > 0) & (`n_hotspots_Non-Urban` > 0)
  )

# 2. Calculate species positions in NMDS space
species_nmds_positions <- species_scores %>%
  rownames_to_column("scientific_name") %>%
  mutate(
    # Distance from center (0,0) - species in center are more generalist
    distance_from_center = sqrt(NMDS1^2 + NMDS2^2),
    # Absolute position on each axis
    abs_nmds1 = abs(NMDS1),
    abs_nmds2 = abs(NMDS2)
  )

# 3. Combine occurrence and NMDS data
candidate_species <- species_occurrence %>%
  left_join(species_nmds_positions, by = "scientific_name") %>%
  filter(in_both_types == TRUE) %>%  # Only species in both site types
  arrange(distance_from_center)  # Species closer to center first

# 4. Set selection criteria for patch occupancy candidates
selection_criteria <- candidate_species %>%
  filter(
    # Present in both site types with reasonable frequency
    min_occupancy >= 0.15,  # At least 10% occupancy in both site types
    total_occupancy >= 0.25,  # At least 15% overall occupancy
    # Not too rare or too common (good for occupancy modeling)
    total_occupancy <= 1,   # Less than 90% occupancy (avoid ubiquitous species)
    # Sufficient observations for modeling
    (`n_observations_Urban` + `n_observations_Non-Urban`) >= 25,
    # Not extreme specialists (closer to NMDS center)
    distance_from_center <= quantile(candidate_species$distance_from_center, 0.7, na.rm = TRUE)
  ) %>%
  # Rank candidates
  mutate(
    # Occupancy balance score (higher = more balanced between site types)
    occupancy_balance = 1 - abs(urban_occupancy - nonurban_occupancy),
    # Centrality score (higher = more central in ordination space)
    centrality_score = 1 - (distance_from_center / max(distance_from_center, na.rm = TRUE)),
    # Combined score
    candidate_score = (occupancy_balance * 0.4) + (centrality_score * 0.3) + (total_occupancy * 0.3)
  ) %>%
  arrange(desc(candidate_score))

# 5. Select top candidates
top_candidates <- head(selection_criteria, 25)  # Top 15 candidates

print("=== TOP CANDIDATE SPECIES FOR PATCH OCCUPANCY ANALYSIS ===")
print(paste("Selected", nrow(top_candidates), "species out of", 
            sum(species_occurrence$in_both_types), "species present in both site types"))
print("")

# Display candidate species with key metrics
candidate_summary <- top_candidates %>%
  select(scientific_name, urban_occupancy, nonurban_occupancy, total_occupancy,
         distance_from_center, occupancy_balance, candidate_score) %>%
  mutate(across(where(is.numeric), round, 3))

print("Candidate Species Summary:")
print(candidate_summary)

# 6. Create visualization of candidate selection

radius <- quantile(candidate_species$distance_from_center, 0.7, na.rm = TRUE)
circle_data <- data.frame(
  x = radius * cos(seq(0, 2*pi, length.out = 100)),
  y = radius * sin(seq(0, 2*pi, length.out = 100))
)

p3 <- ggplot() +
  # All species (background)
  geom_point(data = species_scores, 
             aes(x = NMDS1, y = NMDS2), 
             color = "lightgray", alpha = 0.5, size = 1) +
  
  # Species present in both site types
  geom_point(data = species_nmds_positions %>% 
               filter(scientific_name %in% candidate_species$scientific_name),
             aes(x = NMDS1, y = NMDS2), 
             color = "blue", alpha = 0.6, size = 2) +
  
  # Top candidates
  geom_point(data = species_nmds_positions %>% 
               filter(scientific_name %in% top_candidates$scientific_name),
             aes(x = NMDS1, y = NMDS2), 
             color = "red", size = 3) +
  
  # Add labels for top candidates
  geom_text(data = species_nmds_positions %>% 
              filter(scientific_name %in% head(top_candidates$scientific_name, 8)),
            aes(x = NMDS1, y = NMDS2, label = scientific_name),
            size = 2.5, hjust = 0, vjust = 0, color = "darkred") +
  
  # Add circle showing selection radius
  geom_path(data = circle_data, aes(x = x, y = y),
            color = "black", linetype = "dashed", alpha = 0.5) +
  
  labs(title = "Candidate Species Selection for Patch Occupancy Analysis",
       subtitle = "Red = Top candidates",
       x = "NMDS1", y = "NMDS2") +
  theme_classic() +
  coord_fixed()

print(p3)

# 7. Export candidate species list
write.csv(top_candidates, "patch_occupancy_candidates.csv", row.names = FALSE)

print("")
print("=== SELECTION CRITERIA SUMMARY ===")
print(paste("Minimum occupancy in both site types: ≥10%"))
print(paste("Total occupancy range: 15-80%"))
print(paste("Minimum total observations: ≥20"))
print(paste("Maximum distance from NMDS center: ≤", 
            round(quantile(candidate_species$distance_from_center, 0.7, na.rm = TRUE), 3)))

print("")
print("=== NEXT STEPS ===")
print("1. Review the candidate species list saved to 'patch_occupancy_candidates.csv'")
print("2. Consider ecological knowledge to refine selection")
print("3. Check detection probability assumptions for selected species")
print("4. Proceed with patch occupancy modeling using selected candidates")

# 8. Additional diagnostic plots
p4 <- ggplot(candidate_species, aes(x = urban_occupancy, y = nonurban_occupancy)) +
  geom_point(alpha = 0.6, size = 2) +
  geom_point(data = top_candidates, aes(x = urban_occupancy, y = nonurban_occupancy),
             color = "red", size = 3) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", alpha = 0.5) +
  labs(title = "Species Occupancy: Urban vs Non-Urban",
       subtitle = "Red points = Selected candidates, Dashed line = Equal occupancy",
       x = "Urban Occupancy", y = "Non-Urban Occupancy") +
  theme_classic()

print(p4)