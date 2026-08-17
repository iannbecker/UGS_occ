##############################
#
# eBird Data Filter
# Ian Becker
# August 2026
#
##############################

# This script is used to load in raw eBird data,
# filter it to checklists submitted at my study sites,
# and finally filter down using quality control metrics and to
# the 36 species included in the final within-site analysis

# ============================================================================
# 1. LOAD IN PACKAGES
# ============================================================================

library(auk)
library(tidyverse)
library(sf)

# ============================================================================
# 2. DIRECTORIES AND PREP
# ============================================================================

# Setup directories

ebird_dir  <- "/Users/ianbecker/Desktop/project_code/UGS_occ/data/ebird_data"
input_dir  <- "/Users/ianbecker/Desktop/project_code/UGS_occ/data"
output_dir <- "/Users/ianbecker/Desktop/project_code/UGS_occ/data"
tmp_dir    <- "/Users/ianbecker/Desktop/project_code/UGS_occ/data/ebird_tmp"

# This is for storing filtered files from auk

dir.create(tmp_dir, showWarnings = FALSE)

# County codes for LRGV — used to match file pairs

county_codes <- c("US-TX-215", "US-TX-061", "US-TX-427", "US-TX-489")

# Hidalgo = 215, Cameron = 061, Starr = 427, Willacy = 489

# ============================================================================
# 3. LOAD STUDY SITES AND SPECIES LIST
# ============================================================================

# Load study sites shapefile and convert crs

sites <- st_read(file.path(input_dir, "lrgv_green_spaces_gbif"),
                     quiet = TRUE)
sites_wgs <- st_transform(sites, crs = 4326)

# Load in study species

species_csv   <- read.csv(file.path(input_dir, "lrgv_ugs_species_FILTERED.csv"))
study_species <- species_csv$common_name
study_species <- study_species %>%
  recode(
    "Feral Pigeon"         = "Rock Pigeon (Feral Pigeon)",
    "Domestic Muscovy Duck" = "Muscovy Duck (Domestic type)"
  )

# Check numbers here

nrow(sites)
length(study_species)

# ============================================================================
# 4. FIND EBD AND SAMPLING FILE PAIRS
# ============================================================================

# Find all EBD files (non-sampling)

ebd_files <- list.files(ebird_dir,
                        pattern = "ebd_.*\\.txt$",
                        full.names = TRUE,
                        recursive = TRUE)
ebd_files <- ebd_files[!grepl("_sampling", ebd_files)]

# Find all sampling files

smp_files <- list.files(ebird_dir,
                        pattern = ".*_sampling\\.txt$",
                        full.names = TRUE,
                        recursive = TRUE)

# Check numbers here

length(ebd_files)
length(smp_files)

# ============================================================================
# 5. FILTER AND ZEROFILL EACH COUNTY
# ============================================================================

county_data <- list()

# For loop to process each county's EBD and sampling files

for (i in seq_along(ebd_files)) {
  
  ebd_file <- ebd_files[i]
  smp_file <- smp_files[i]
  
  cat("Processing:", basename(ebd_file), "\n")
  
  # Output paths for filtered files
  out_ebd <- file.path(tmp_dir, paste0("filtered_ebd_", i, ".txt"))
  out_smp <- file.path(tmp_dir, paste0("filtered_smp_", i, ".txt"))
  
  # Build auk filter pipeline
  # auk works by writing filtered versions of the raw .txt files to disk
  # before reading into R — this is much faster than reading everything first
  # Each step below adds a filter condition:
  tryCatch({
  data <- auk_ebd(ebd_file, file_sampling = smp_file) %>%
      
      # Complete checklists only — ensures species not reported = not detected
      # rather than not looked for, which is critical for occupancy modeling
      auk_complete() %>%
      
      # Stationary or traveling protocols only — excludes incidental observations
      # and other non-standard effort types
      auk_protocol(c("Stationary", "Traveling")) %>%
      
      # Filter to study species before writing to disk — dramatically reduces
      # file size and zerofill memory footprint vs filtering after
    #  auk_species(study_species) %>%
      
      # Max 6 hours duration — longer checklists cover too large an area
      # and are less comparable across observers
      auk_duration(c(0, 360)) %>%
      
      # Max 10 km distance — caps spatial extent of traveling counts
      auk_distance(c(0, 10)) %>%
      
      # Study period
      auk_date(c("2015-01-01", "2025-12-31")) %>%
      
      # Write filtered EBD and sampling files to tmp directory
      # This is the bottleneck step — reads and filters the raw .txt files
      # line by line, which can take several minutes per county file
      auk_filter(file          = out_ebd,
                 file_sampling = out_smp,
                 overwrite     = TRUE)
    
    cat("  auk_filter complete — reading and zerofilling...\n")
    
    # auk_zerofill joins the filtered EBD and sampling files
    # and adds explicit zero records for species not detected on each checklist
    # collapse = TRUE returns one row per checklist-species combination
      zf = auk::read_ebd(out_ebd)
    
   #   zf <- auk_zerofill(out_ebd, out_smp, collapse = TRUE)
    
    # Deduplicate — remove any duplicate checklist-species records
   zf <- zf %>% distinct(checklist_id, scientific_name, .keep_all = TRUE)
    
   county_data[[i]] <- zf
   # cat("  Checklist-species records after filtering:", nrow(zf), "\n\n")
    
  }, error = function(e) {
    cat("  ERROR:", e$message, "\n\n")
  })
}

# Combine all counties and further filter

ebd_combined <- bind_rows(county_data)

# Cap at 10 observers 

ebd_combined <- ebd_combined %>%
  filter(number_observers <= 10 | is.na(number_observers))

# Get common name from taxonomy

ebd_combined <- ebd_combined %>%
  left_join(ebird_taxonomy %>% 
              select(scientific_name, common_name),
            by = "scientific_name")

# Total observations (not study site filtered)

nrow(ebd_combined)

ebd_summary <- ebd_combined %>%
  group_by(common_name) %>%
  summarise(count = n()) %>%
  filter(count >= 10)

# ============================================================================
# 6. ASSIGN CHECKLISTS TO STUDY SITES
# ============================================================================

# Change to sf 

ebd_sf <- st_as_sf(ebd_combined,
                   coords = c("longitude", "latitude"),
                   crs = 4326)

# Count checklists per site

obs_in_sites <- st_intersects(ebd_sf, sites_wgs)

# Assign site_id to each checklist based on intersection

ebd_sf$site_id <- sapply(obs_in_sites, function(x) {
  if (length(x) > 0) return(sites_wgs$site_id[x[1]])
  else return(NA)
})

# Get rid of checklists not in any sites

ebd_in_sites <- ebd_sf %>% filter(!is.na(site_id))

# Records within sites and distinct sites with data

nrow(ebd_in_sites)
n_distinct(ebd_in_sites$site_id)

# Add date column and drop geometry

ebd_study <- ebd_in_sites %>%
  st_drop_geometry() %>%
  mutate(
    observation_date = as.Date(observation_date),
    year             = as.integer(format(observation_date, "%Y")),
    month            = as.integer(format(observation_date, "%m"))
  )

# ============================================================================
# 8. SUMMARY
# ============================================================================

# Species in study list not found in eBird data

missing_species <- study_species[!study_species %in% ebd_study$common_name]
if (length(missing_species) > 0) {
  cat("Species not found in eBird data:\n")
  for (sp in missing_species) cat(" -", sp, "\n")
  cat("\n")
}

# Detection summary by species

print(ebd_study %>%
        group_by(common_name) %>%
        summarise(
          n_checklists = n(),
          n_detections = sum(species_observed),
          detection_rate = round(mean(species_observed), 3),
          .groups = "drop"
        ) %>%
        arrange(desc(n_detections)))

# ============================================================================
# 9. SAVE
# ============================================================================

# save filtered checklist data

write.csv(ebd_study,
          file.path(output_dir, "ebird_filtered_lrgv_with_sites.csv"),
          row.names = FALSE)
