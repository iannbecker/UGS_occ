##############################
#
# OSM Green Spaces Data Pull
# 12/20/2025
# Ian Becker
#
##############################

# Pull various types of urban green spaces from OpenStreetMap
# for the LRGV study area

library(osmdata)
library(sf)
library(dplyr)

####################
#   Setup
####################

# Load LRGV urban areas to define bounding box
lrgv_urban <- st_read("lrgv_urban_areas")

# Get bounding box (OSM needs WGS84)
if (st_crs(lrgv_urban)$input != "EPSG:4326") {
  lrgv_urban_wgs84 <- st_transform(lrgv_urban, crs = 4326)
} else {
  lrgv_urban_wgs84 <- lrgv_urban
}

bbox <- st_bbox(lrgv_urban_wgs84)

cat("Study area bounding box (WGS84):\n")
cat("  West:", bbox$xmin, "\n")
cat("  East:", bbox$xmax, "\n")
cat("  South:", bbox$ymin, "\n")
cat("  North:", bbox$ymax, "\n\n")

####################
#   Pull OSM Data
####################

cat("Pulling green spaces from OpenStreetMap...\n")
cat("This may take a few minutes...\n\n")

all_green_spaces <- list()

# 1. Parks and gardens
cat("1. Pulling parks and gardens...\n")
tryCatch({
  parks_osm <- opq(bbox) %>%
    add_osm_feature(key = "leisure", value = c("park", "garden", "nature_reserve")) %>%
    osmdata_sf()
  
  if (!is.null(parks_osm$osm_polygons) && nrow(parks_osm$osm_polygons) > 0) {
    parks <- parks_osm$osm_polygons %>%
      select(osm_id, name, leisure) %>%
      mutate(type = "park/garden")
    all_green_spaces[["parks"]] <- parks
    cat("  Found", nrow(parks), "parks/gardens\n")
  } else {
    cat("  No parks/gardens found\n")
  }
}, error = function(e) {
  cat("  Error:", e$message, "\n")
})

# 2. University/college campuses
cat("2. Pulling university/college campuses...\n")
tryCatch({
  campuses_osm <- opq(bbox) %>%
    add_osm_feature(key = "education", value = c("university", "college")) %>%
    osmdata_sf()
  
  if (!is.null(campuses_osm$osm_polygons) && nrow(campuses_osm$osm_polygons) > 0) {
    campuses <- campuses_osm$osm_polygons %>%
      select(osm_id, name, amenity) %>%
      mutate(type = "campus")
    all_green_spaces[["campuses"]] <- campuses
    cat("  Found", nrow(campuses), "campuses\n")
  } else {
    cat("  No campuses found\n")
  }
}, error = function(e) {
  cat("  Error:", e$message, "\n")
})

# 3. Cemeteries
cat("3. Pulling cemeteries...\n")
tryCatch({
  cemeteries_osm <- opq(bbox) %>%
    add_osm_feature(key = "landuse", value = "cemetery") %>%
    osmdata_sf()
  
  if (!is.null(cemeteries_osm$osm_polygons) && nrow(cemeteries_osm$osm_polygons) > 0) {
    cemeteries <- cemeteries_osm$osm_polygons %>%
      select(osm_id, name) %>%
      mutate(type = "cemetery")
    all_green_spaces[["cemeteries"]] <- cemeteries
    cat("  Found", nrow(cemeteries), "cemeteries\n")
  } else {
    cat("  No cemeteries found\n")
  }
}, error = function(e) {
  cat("  Error:", e$message, "\n")
})

# 4. Golf courses
cat("4. Pulling golf courses...\n")
tryCatch({
  golf_osm <- opq(bbox) %>%
    add_osm_feature(key = "leisure", value = "golf_course") %>%
    osmdata_sf()
  
  if (!is.null(golf_osm$osm_polygons) && nrow(golf_osm$osm_polygons) > 0) {
    golf <- golf_osm$osm_polygons %>%
      select(osm_id, name) %>%
      mutate(type = "golf_course")
    all_green_spaces[["golf"]] <- golf
    cat("  Found", nrow(golf), "golf courses\n")
  } else {
    cat("  No golf courses found\n")
  }
}, error = function(e) {
  cat("  Error:", e$message, "\n")
})

# 5. Recreation grounds, pitches, playgrounds
cat("5. Pulling recreation areas...\n")
tryCatch({
  recreation_osm <- opq(bbox) %>%
    add_osm_feature(key = "leisure", value = c("recreation_ground", "pitch", "playground")) %>%
    osmdata_sf()
  
  if (!is.null(recreation_osm$osm_polygons) && nrow(recreation_osm$osm_polygons) > 0) {
    recreation <- recreation_osm$osm_polygons %>%
      select(osm_id, name, leisure) %>%
      mutate(type = "recreation")
    all_green_spaces[["recreation"]] <- recreation
    cat("  Found", nrow(recreation), "recreation areas\n")
  } else {
    cat("  No recreation areas found\n")
  }
}, error = function(e) {
  cat("  Error:", e$message, "\n")
})

# 6. Sports centres
cat("6. Pulling sports centres...\n")
tryCatch({
  sports_osm <- opq(bbox) %>%
    add_osm_feature(key = "leisure", value = "sports_centre") %>%
    osmdata_sf()
  
  if (!is.null(sports_osm$osm_polygons) && nrow(sports_osm$osm_polygons) > 0) {
    sports <- sports_osm$osm_polygons %>%
      select(osm_id, name) %>%
      mutate(type = "sports_centre")
    all_green_spaces[["sports"]] <- sports
    cat("  Found", nrow(sports), "sports centres\n")
  } else {
    cat("  No sports centres found\n")
  }
}, error = function(e) {
  cat("  Error:", e$message, "\n")
})

# 7. Plazas and town squares
cat("7. Pulling plazas...\n")
tryCatch({
  plazas_osm <- opq(bbox) %>%
    add_osm_feature(key = "place", value = c("square", "plaza")) %>%
    osmdata_sf()
  
  if (!is.null(plazas_osm$osm_polygons) && nrow(plazas_osm$osm_polygons) > 0) {
    plazas <- plazas_osm$osm_polygons %>%
      select(osm_id, name) %>%
      mutate(type = "plaza")
    all_green_spaces[["plazas"]] <- plazas
    cat("  Found", nrow(plazas), "plazas\n")
  } else {
    cat("  No plazas found\n")
  }
}, error = function(e) {
  cat("  Error:", e$message, "\n")
})

####################
#   Combine and Clean
####################

cat("\nCombining all green spaces...\n")

if (length(all_green_spaces) > 0) {
  
  # Standardize columns across all datasets
  standardize_cols <- function(df) {
    if (!"name" %in% names(df)) df$name <- NA_character_
    df %>% select(osm_id, name, type, geometry)
  }
  
  all_green_spaces_clean <- lapply(all_green_spaces, standardize_cols)
  
  # Combine all
  combined_green_spaces <- bind_rows(all_green_spaces_clean)
  
  # Transform to UTM for consistency
  combined_green_spaces_utm <- st_transform(combined_green_spaces, crs = 32614)
  
  # Add area and filter by minimum size
  combined_green_spaces_utm <- combined_green_spaces_utm %>%
    mutate(
      area_m2 = as.numeric(st_area(geometry)),
      area_ha = area_m2 / 10000
    ) %>%
    # Filter to reasonably sized spaces (>0.1 ha = 1000 m2)
    filter(area_ha > 0.1) %>%
    # Add unique site ID
    mutate(site_id = row_number())
  
  cat("\n=== OSM GREEN SPACES SUMMARY ===\n")
  cat("Total features:", nrow(combined_green_spaces_utm), "\n")
  cat("Total area:", round(sum(combined_green_spaces_utm$area_ha), 1), "hectares\n\n")
  
  # Summary by type
  cat("By type:\n")
  type_summary <- combined_green_spaces_utm %>%
    st_drop_geometry() %>%
    group_by(type) %>%
    summarise(
      count = n(),
      total_area_ha = round(sum(area_ha), 1),
      mean_area_ha = round(mean(area_ha), 2),
      .groups = "drop"
    ) %>%
    arrange(desc(count))
  
  print(type_summary)
  
  # Size distribution
  cat("\nSize distribution:\n")
  combined_green_spaces_utm <- combined_green_spaces_utm %>%
    mutate(
      size_class = cut(area_ha, 
                       breaks = c(0, 1, 5, 10, 50, 100, Inf),
                       labels = c("<1 ha", "1-5 ha", "5-10 ha", "10-50 ha", "50-100 ha", ">100 ha"))
    )
  
  size_dist <- table(combined_green_spaces_utm$size_class)
  for (i in 1:length(size_dist)) {
    cat("  ", names(size_dist)[i], ":", size_dist[i], "sites\n")
  }
  
  # List named sites
  named_sites <- combined_green_spaces_utm %>%
    filter(!is.na(name), name != "") %>%
    arrange(desc(area_ha))
  
  if (nrow(named_sites) > 0) {
    cat("\nLargest named sites:\n")
    top_sites <- head(named_sites, 20)
    for (i in 1:nrow(top_sites)) {
      cat(sprintf("  %2d. %-40s (%s, %.1f ha)\n",
                  i,
                  top_sites$name[i],
                  top_sites$type[i],
                  top_sites$area_ha[i]))
    }
  }
  
  ####################
  #   Filter to Manageable Sites
  ####################
  
  cat("\n=== FILTERING TO MANAGEABLE NUMBER OF SITES ===\n\n")
  
  osm_all <- combined_green_spaces_utm  # Keep the full dataset
  
  # 1. Keep ONLY named sites (must have identifiable name)
  cat("1. Filtering to named sites only...\n")
  osm_filtered <- osm_all %>%
    filter(!is.na(name) & name != "")
  
  cat("   Kept", nrow(osm_filtered), "named sites (removed", nrow(osm_all) - nrow(osm_filtered), "unnamed)\n\n")
  
  # 2. Remove plazas and sports centers (not true green spaces)
  cat("2. Removing plazas and sports centers...\n")
  osm_filtered <- osm_filtered %>%
    filter(!type %in% c("plaza", "sports_centre"))
  
  cat("   Remaining:", nrow(osm_filtered), "\n\n")
  
  # 3. Remove ONLY very large wildlife refuges (>200 ha)
  # Keep smaller accessible wildlife refuge parcels
  cat("3. Filtering large wildlife refuges (keeping parcels <200 ha)...\n")
  osm_filtered <- osm_filtered %>%
    filter(!(grepl("Wildlife Refuge", name, ignore.case = TRUE) & area_ha > 200))
  
  cat("   Remaining:", nrow(osm_filtered), "\n\n")
  
  # 4. Size criteria - keep sites between 1-100 ha
  cat("4. Filtering by size (1-100 hectares)...\n")
  osm_filtered <- osm_filtered %>%
    filter(area_ha >= 1, area_ha <= 100)
  
  cat("   Remaining:", nrow(osm_filtered), "\n\n")
  
  # 5. Keep only areas within census urban boundaries
  cat("5. Filtering to areas within urban boundaries...\n")
  
  # Ensure same CRS
  lrgv_urban_check <- lrgv_urban_wgs84
  if (st_crs(osm_filtered) != st_crs(lrgv_urban_check)) {
    lrgv_urban_check <- st_transform(lrgv_urban_check, st_crs(osm_filtered))
  }
  
  # Find intersection with urban areas
  urban_intersects <- st_intersects(osm_filtered, lrgv_urban_check, sparse = FALSE)
  osm_filtered <- osm_filtered[apply(urban_intersects, 1, any), ]
  
  cat("   Remaining:", nrow(osm_filtered), "\n\n")
  
  # Final sites = all that passed filters
  final_sites <- osm_filtered %>%
    mutate(site_id = row_number())
  
  ####################
  #   Final Summary
  ####################
  
  cat("=== FILTERED SITE SUMMARY ===\n")
  cat("Started with:", nrow(osm_all), "sites\n")
  cat("Final count:", nrow(final_sites), "\n")
  cat("Total area:", round(sum(final_sites$area_ha), 1), "hectares\n\n")
  
  # By type
  cat("By type:\n")
  type_summary_final <- final_sites %>%
    st_drop_geometry() %>%
    group_by(type) %>%
    summarise(
      count = n(),
      total_area_ha = round(sum(area_ha), 1),
      mean_area_ha = round(mean(area_ha), 2),
      .groups = "drop"
    ) %>%
    arrange(desc(count))
  
  print(type_summary_final)
  
  # Size distribution
  cat("\nSize distribution:\n")
  size_dist_final <- table(final_sites$size_class)
  for (i in 1:length(size_dist_final)) {
    cat("  ", names(size_dist_final)[i], ":", size_dist_final[i], "sites\n")
  }
  
  # Named vs unnamed
  named_count <- sum(!is.na(final_sites$name) & final_sites$name != "")
  cat("\nNamed sites:", named_count, "\n")
  cat("Unnamed sites:", nrow(final_sites) - named_count, "\n")
  
  # Top sites
  cat("\nTop 20 largest sites:\n")
  top_sites_final <- final_sites %>%
    arrange(desc(area_ha)) %>%
    head(20)
  
  for (i in 1:nrow(top_sites_final)) {
    site_name <- ifelse(is.na(top_sites_final$name[i]) | top_sites_final$name[i] == "", 
                        "Unnamed", top_sites_final$name[i])
    cat(sprintf("  %2d. %-50s (%s, %.1f ha)\n",
                i, site_name, top_sites_final$type[i], top_sites_final$area_ha[i]))
  }
  
  ####################
  #   Save Results
  ####################
  
  cat("\nSaving results...\n")
  
  # Save ALL sites (unfiltered)
  st_write(osm_all, "lrgv_osm_green_spaces_all.shp", delete_dsn = TRUE)
  cat("Saved: lrgv_osm_green_spaces_all.shp (all", nrow(osm_all), "sites)\n")
  
  # Save FILTERED sites (for occupancy modeling)
  st_write(final_sites, "lrgv_osm_green_spaces_filtered.shp", delete_dsn = TRUE)
  cat("Saved: lrgv_osm_green_spaces_filtered.shp (", nrow(final_sites), "sites)\n")
  
  # Save summary tables
  summary_all <- st_drop_geometry(osm_all)
  write.csv(summary_all, "lrgv_osm_green_spaces_all_summary.csv", row.names = FALSE)
  cat("Saved: lrgv_osm_green_spaces_all_summary.csv\n")
  
  summary_filtered <- st_drop_geometry(final_sites)
  write.csv(summary_filtered, "lrgv_osm_green_spaces_filtered_summary.csv", row.names = FALSE)
  cat("Saved: lrgv_osm_green_spaces_filtered_summary.csv\n")
  
  ####################
  #   Recommendation
  ####################
  
  cat("\n=== RECOMMENDATION ===\n")
  if (nrow(final_sites) > 200) {
    cat("⚠ Still have", nrow(final_sites), "sites - consider:\n")
    cat("   - Increasing minimum size to 3 or 5 ha\n")
    cat("   - Focusing on named sites only\n")
    cat("   - Selecting specific park types\n")
  } else if (nrow(final_sites) < 30) {
    cat("⚠ Only", nrow(final_sites), "sites - consider:\n")
    cat("   - Lowering minimum size threshold\n")
    cat("   - Including more unnamed sites\n")
    cat("   - Adding sites from clustering algorithm\n")
  } else {
    cat("✓ Good sample size (", nrow(final_sites), "sites)\n")
    cat("  This should work well for occupancy modeling\n")
  }
  
  cat("\n=== OSM DATA PULL COMPLETE ===\n")
  cat("Use 'lrgv_osm_green_spaces_filtered.shp' for occupancy modeling\n")
  
} else {
  cat("\n⚠ No green spaces found in OSM for this area\n")
  cat("You may need to use the clustering approach instead\n")
}

