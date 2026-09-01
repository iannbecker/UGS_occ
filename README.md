# Data and code for: Integrating community science datasets reveals scale dependent habitat selection in urban bird communities
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
Greenspaces in urban areas contribute to the long-term conservation of biodiversity. Many species rely on these habitat patches for foraging, breeding, and migratory stopover. Despite increasing research quantifying when, and to what extent, birds use urban greenspaces, the factors determining how individuals are distributed within a greenspace have seldom been explored. This information is invaluable for researchers and urban habitat practitioners, who must consider how site level habitat features impact the community of species that use an urban greenspace. We integrate semi-structured, checklist data from eBird with opportunistic, georeferenced observation data from iNaturalist to compare habitat selection patterns between landscape and site-level spatial scales for 36 urban bird species in the Lower Rio Grande Valley in Texas. Our results indicate consistent attraction to tree cover and flooded vegetation as habitat features for a majority of species both among and within urban greenspaces. In contrast, we found evidence of scale reversal in response to water cover, with many species positively associated with water cover at a landscape-scale but negatively associated with water cover once within a greenspace. Species in our study area utilized an average of 30% of the available area of a greenspace, highlighting the potential for these habitat patches to serve as ecological traps. Based on our findings, we recommend habitat management prioritize expanding resources within existing greenspaces and that such management can be evaluated on a species-specific basis. We propose and discuss a novel framework that integrates large-scale community science datasets to develop species-specific conservation plans at both landscape and site levels.


