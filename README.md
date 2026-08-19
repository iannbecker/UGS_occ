# Data and code for: TBD
**Citation:** TBD

## Data
TBD - These files will likely be stored in zenodo due to managing file size

## Scripts

### Data Prep/Cleaning
These scripts are for data preparation and cleaning that occurred prior to the main analysis in the manuscript.

| Script | Description |
|--------|-------------|
| `UGS_pull_filter.r` | Initial OSM query for urban greenspaces and filtering |
| `covariate_prep.r` | Extracts land cover data for urban green spaces and creates covariate matrix for modeling |
| `ebird_data_filter.r` | Initial filter of raw eBird data for modeling |
| `ebird_CHECKLIST_detection_matrix.r` | Builds detection matrices for 36 species based on eBird data used for landscape-level modeling |
| `within_site_viability.r` | Used to find viable species-site combinations for site-level models |
| `gbif_detection_matrices.r` | Builds detection matrices for 36 species based on iNaturalist data used for site-level modeling |

### Analysis
These scripts are for the main analysis described in the manuscript.

| Script | Description |
|--------|-------------|
| `landscape_ebird_model.r` | Runs landscape-level occupancy models using eBird data for 36 species|
| `site_inat_model.r` | Runs site-level models using iNaturalist data for all viable species-site combinations |
| `scale_comparison_ebird_gbif.r` | Compares eBird landscape-level results to iNat site-level results|
| `within_site_sensitivity_analysis.r` | Sensitivity analysis for site-level models; testing robustness of iNat spatial uncertainty |

### Figures
These scripts are for figures in the main body of the manuscript. 

| Script | Description |
|--------|-------------|
| `Figure1_StudyArea_Scale.r` | Map of study area depicting conceptual scale comparison for our study |
| `Figure2_LandscapeOccupancy.r` | Summary of landscape-level occupancy trends |
| `Figure3_Landscape_vs_Site.r` | Comparison of landscape_level and site_level coefficients by species-covariate-site combination |
| `Figure4_BarPlot.r` | Creates bar plot of mean area used per site by species + pie chart showing overall average area usage |
| `Figure5_SpeciesExample.r` | Makes both Figure 5 and Figure S2 in the manuscript; creates species example showing combined landscape/site-level response + landscape detection map + within-site detection map|


## Abstract
TBD


