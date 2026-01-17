##############################
#
# By species site count
# Ian Becker
# January 2026
#
##############################

# This script counts observations per species by site
# to try and use fine-scale single site occupancy models 

library(dplyr)

# Load the already-processed data with site assignments
inat_with_sites <- read.csv("inat_observations_with_sites.csv")

# Define the 15 study species
waterbirds <- c(
  "Black-bellied Whistling-Duck",
  "Green Heron",
  "Neotropic Cormorant",
  "Great Blue Heron",
  "Great Egret"
)

generalists <- c(
  "Great-tailed Grackle",
  "Northern Mockingbird",
  "White-winged Dove",
  "Inca Dove",
  "Mourning Dove"
)

specialists <- c(
  "Great Kiskadee",
  "Tropical Kingbird",
  "Green Jay",
  "Long-billed Thrasher",
  "Golden-fronted Woodpecker"
)

study_species <- c(waterbirds, generalists, specialists)

# Filter to study species and count per site
obs_per_site <- inat_with_sites %>%
  filter(common_name %in% study_species) %>%
  group_by(common_name, site_id) %>%
  summarise(n_obs = n(), .groups = "drop")

# Filter to ≥10 observations
sufficient_data <- obs_per_site %>%
  filter(n_obs >= 10) %>%
  arrange(common_name, desc(n_obs))

# Summary by species
cat("=== SITES WITH ≥10 OBSERVATIONS BY SPECIES ===\n\n")

for (sp in study_species) {
  sp_data <- sufficient_data %>% filter(common_name == sp)
  
  if (nrow(sp_data) > 0) {
    cat(sp, ":", nrow(sp_data), "sites\n")
    cat("  Total obs across these sites:", sum(sp_data$n_obs), "\n")
    cat("  Site IDs:", paste(sp_data$site_id, collapse = ", "), "\n")
    cat("  Obs per site:", paste(sp_data$n_obs, collapse = ", "), "\n\n")
  } else {
    cat(sp, ": 0 sites with ≥10 obs\n\n")
  }
}

# Which sites appear most often?
cat("=== SITES WITH MULTIPLE SPECIES ≥10 OBS ===\n\n")
site_counts <- sufficient_data %>%
  group_by(site_id) %>%
  summarise(
    n_species = n(),
    species = paste(common_name, collapse = ", "),
    total_obs = sum(n_obs)
  ) %>%
  arrange(desc(n_species))

print(site_counts, n = 20)

# Save
write.csv(sufficient_data, "species_sites_10plus_obs.csv", row.names = FALSE)
write.csv(site_counts, "sites_multispecies_10plus.csv", row.names = FALSE)