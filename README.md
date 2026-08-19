# UGS_occ
Leveraging eBird and iNat data to explore patch occupancy dynamics of birds in urban areas

BELOW IS A TEMPLATE FROM ANOTHER REPO - THIS IS A WIP

# Data and code for: Understanding how to implement effective biodiversity monitoring on college campuses through community science 
**Citation:** Becker, I., & Contina, A. (2026). Data and code from: Understanding how to implement effective biodiversity monitoring on college campuses through community science (v1.0.0). Zenodo. https://doi.org/10.5281/zenodo.19930592

## Data
These files can be used to recreate the analyses described below. For privacy purposes we do not include here movement data which included raw location data specific to individual users. Instead here we include an example simulated dataset which can be used to recreate analyses. 

| File Name | Description |
|--------|-------------|
| `campus_data_with_pca.csv` | Master dataset for modelling. Includes all variables by campus included in models including PCA scores |
| `campus_polygons_MASTER` | Shapefile with campus boundaries (n=153) |
| `movement_network_EXAMPLE.csv` | Simulated dataset used as an example for movement analyses to protect eBird user privacy |
| `sensitivity_test.csv` | Dataset used for geographic sensitivity analysis (n=208) |
| `unique_visitor_comparison.csv` | Dataset used for visitor behavior comparison between campuses and hotspots |

## Scripts

### Data Prep/Cleaning
These scripts are for data preparation and cleaning that occurred prior to the main analysis in the manuscript.

| Script | Description |
|--------|-------------|
| `campus_data_pull.R` | Query campus data from the IPEDS database |
| `checklist_extract.R` | Extract checklist data from eBird and count checklists submitted within campus boundaries |
| `environmental_data_prep.R` | Prep environmental covariates for modelling (NLCD Land Cover & %Imperviousness) |
| `location_network.R` | Gather weekly movement window data for observers in our study. All data used has been anonymized to remove PID |
| `osm_polygon_pull.R` | Query OSM database for campus polygons/boundaries |
| `visitor_data.R` | Prep visitor data for behavior comparison between campuses and hotspots |

### Analysis
These scripts are for the main analysis described in the manuscript.

| Script | Description |
|--------|-------------|
| `bird_community_analysis.R` | Calculates species richness, similarity, and rare species occurrence between campuses and hotspots. Also makes the venn diagram for Figure 6 in the manuscript |
| `checklists_visitor_analysis_NOMOVEMENT.R` | Analyzes observer dynamics (not including movement patterns) between campuses and hotspots |
| `movement_stats.R` | Analyzes movement patterns and behavior between campus and hotspot observers |
| `PCA_models.R` | Prep PCA variables, run PCA, and subsequently run models based on competing hypotheses (see manuscript) |

### Figures
These scripts are for figures in the main body of the manuscript as well as the supplemental information indicated as either Figure# (main body) or Figure S# (supplemental).

| Script | Description |
|--------|-------------|
| `Figure1_study_area.r` | Map of study area including campuses sized by checklist count |
| `Figure3_response.r` | Makes 3 figures (3, S6, S7) showing the relationship between checklists and covariates |
| `Figure4_bargraphs.r` | Bar graphs showing submission behavior differences between campus and hotspot observers |
| `Figure5_movement_map.r` | THIS IS STILL A WIP |
| `FigureS2_distribution.r` | WIP - Raw distribution of observers and log transformed distribution of observers |
| `FigureS3_hotspot_map.r` | Map of study area including all included hotspots sized by checklist count |
| `FigureS4_sensitivity_map.r` | Map of campuses used in OR, NJ, and AL for our geographic sensitivity analysis |
| `FigureS5_PCA.r` | PC1 vs. PC2 plotted relationship |

## Abstract
The development of urban areas has rapidly altered native ecosystems, leaving thousands of species with uncertain population trajectories. As a result, understanding the role of green spaces within urban areas has become increasingly more important for maintaining biodiversity. College campuses provide an often overlooked source of urban green space, positioned as an educational body to both contribute to and educate on biodiversity in human-dominated landscapes. Community science databases provide a metric for understanding how humans utilize these unique green spaces from a socio-ecological point of view. Here we focus on eBird, a bird centric community science database, to analyze community science engagement across 153 college campuses through the south central United States. We test 4 competing hypotheses to understand the role of institution-level characteristics and landscape environmental variables on campus engagement. Additionally, we analyze the visit and movement dynamics of campus observers relative to other community science hotspots to explore the significance of site-dependent engagement. Institutional scale emerged as a strong predictor of engagement, while observer patterns differed distinctly between campuses and other nearby green spaces. Campus engagement was largely driven by small communities of highly active observers who tended to move within local patches relative to the more dispersed movement dynamics of other green space observers. However, the campus bird community overlapped heavily with other green spaces, suggesting that ecological aspects of these habitat patches may operate independent of forces influencing people usage. Together, our findings highlight how institutional context drives community science engagement on college campuses. Through targeted academic programming, campuses have the potential to expand community science engagement outside the classroom, involving the community in the active conservation of avian diversity. 
