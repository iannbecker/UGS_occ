##############################
#
# Data prep
# 12/15/2025
# Ian Becker
#
##############################

library(dplyr)
library(ggplot2)
library(readr)
library(sf)
library(terra)
library(tigris)

options(tigris_use_cache = TRUE)

#################
#  BCRs
#################

# Get Texas shapefile from tigris

texas <- states(cb = TRUE) %>% filter(NAME == "Texas")

# Load in BCR shapefile

gdb_path <- "data/NABCI_ecoregion.gdb"
bcr_shapefile <- st_read(gdb_path, layer = "BCR_Terrestrial_Master")

# Transform to match CRS

bcr_shapefile <- st_transform(bcr_shapefile, st_crs(texas))

# Adjusting geometry to clip to Texas

sf_use_s2(FALSE)  # Turn off spherical geometry for intersection
bcr_shapefile <- st_make_valid(bcr_shapefile)
texas <- st_make_valid(texas)

# Clip BCRs to Texas boundary

bcr_texas <- st_intersection(bcr_shapefile, texas)
sf_use_s2(TRUE)  # Turn spherical geometry back on

# Save shapefile

st_write(bcr_texas, "data/bcr_texas.shp", delete_dsn = TRUE)
