##############################
#
# Green Space Filtering
# 12/20/2025
# Ian Becker
#
##############################

# I am using the script to manually filter through my OSM green spaces
# as well as manually add in anything missing

#################
#   Filter
################

library(sf)
library(mapview)
library(dplyr)

# Load in green spaces

osm_sites <- st_read("lrgv_osm_green_spaces_filtered")

# Interactive map with labels

m <- mapview(combined_sites, 
             zcol = "type",
             label = paste0("ID:", osm_sites$site_id, " - ", 
                            osm_sites$name, " (", 
                            round(osm_sites$area_ha, 1), " ha)"),
             layer.name = "Site Type")

# View the map

m

# As you click through, write down site IDs to REMOVE in a text file
# Then:

sites_to_remove <- c(132, 134, 154, 187, 197, 157, 155, 156, 202, 198, 201, 200, 203, 199)  # IDs you don't want

osm_final <- osm_sites %>%
  filter(!site_id %in% sites_to_remove)

st_write(osm_final, "lrgv_osm_green_spaces_filtered.shp", delete_dsn = TRUE)


###################
#   Part 2: Add Manually Drawn Sites
####################

cat("\n=== ADDING MANUALLY DIGITIZED SITES ===\n")

# Load manually drawn sites (from geojson.io or QGIS)

manual_sites <- st_read("map.geojson")  # Adjust filename as needed

cat("Loaded", nrow(manual_sites), "manually added sites\n")

# Ensure same CRS as filtered OSM sites
if (st_crs(manual_sites) != st_crs(osm_final)) {
  manual_sites <- st_transform(manual_sites, st_crs(osm_final))
  cat("Transformed manual sites to match OSM CRS\n")
}

# Calculate area and size class for manual sites
manual_sites <- manual_sites %>%
  mutate(
    area_m2 = as.numeric(st_area(geometry)),
    area_ha = area_m2 / 10000,
    size_class = cut(area_ha, 
                     breaks = c(0, 1, 5, 10, 50, 100, Inf),
                     labels = c("<1 ha", "1-5 ha", "5-10 ha", 
                                "10-50 ha", "50-100 ha", ">100 ha"))
  )

cat("\nManual sites summary:\n")
cat("  Total area:", round(sum(manual_sites$area_ha), 1), "hectares\n")
cat("  Size range:", round(min(manual_sites$area_ha), 2), "-", 
    round(max(manual_sites$area_ha), 2), "ha\n\n")

# Standardize columns to match OSM sites
# Check what columns OSM sites have
cat("OSM sites columns:", paste(names(osm_final), collapse = ", "), "\n")
cat("Manual sites columns:", paste(names(manual_sites), collapse = ", "), "\n\n")

# Create standardized manual sites dataframe
manual_sites_std <- manual_sites %>%
  mutate(
    osm_id = paste0("manual_", row_number()),  # Create unique IDs
    site_id = NA,  # Will renumber after combining
    source = "manual"
  ) %>%
  # Make sure it has the same columns as osm_final
  select(any_of(names(osm_final)), geometry)

# Add any missing columns with NA
missing_cols <- setdiff(names(osm_final), names(manual_sites_std))
for (col in missing_cols) {
  manual_sites_std[[col]] <- NA
}

# Reorder columns to match osm_final
manual_sites_std <- manual_sites_std[, names(osm_final)]

####################
#   Combine OSM + Manual
####################

cat("Combining OSM and manual sites...\n")

# Add source column to OSM sites if not present
if (!"source" %in% names(osm_final)) {
  osm_final$source <- "OSM"
}

# Combine
combined_sites <- bind_rows(
  osm_final,
  manual_sites_std
) %>%
  mutate(site_id = row_number())  # Renumber all sites sequentially

cat("\n=== COMBINED SITES SUMMARY ===\n")
cat("OSM sites:", sum(combined_sites$source == "OSM"), "\n")
cat("Manual sites:", sum(combined_sites$source == "manual"), "\n")
cat("Total sites:", nrow(combined_sites), "\n")
cat("Total area:", round(sum(combined_sites$area_ha), 1), "hectares\n\n")

# Summary by source and type
cat("By source and type:\n")
source_summary <- combined_sites %>%
  st_drop_geometry() %>%
  group_by(source, type) %>%
  summarise(
    count = n(),
    total_area_ha = round(sum(area_ha, na.rm = TRUE), 1),
    .groups = "drop"
  ) %>%
  arrange(source, desc(count))

print(source_summary)

# Show manually added sites
cat("\n=== MANUALLY ADDED SITES ===\n")
manual_summary <- combined_sites %>%
  filter(source == "manual") %>%
  st_drop_geometry() %>%
  select(site_id, name, type, area_ha) %>%
  arrange(desc(area_ha))

print(manual_summary)

####################
#   Save Final Results
####################

cat("\nSaving final combined results...\n")

# Save combined sites (FINAL version for occupancy modeling)
st_write(combined_sites, "lrgv_green_spaces_combined.shp", delete_dsn = TRUE)
cat("Saved: lrgv_green_spaces_final.shp\n")

# Save summary table
summary_final <- st_drop_geometry(combined_sites)
write.csv(summary_final, "lrgv_green_spaces_final_summary.csv", row.names = FALSE)
cat("Saved: lrgv_green_spaces_final_summary.csv\n")

cat("\n=== FINAL SITE PREPARATION COMPLETE ===\n")
cat("Total sites ready for occupancy modeling:", nrow(combined_sites), "\n")
cat("Use 'lrgv_green_spaces_final.shp' for all downstream analyses\n")






