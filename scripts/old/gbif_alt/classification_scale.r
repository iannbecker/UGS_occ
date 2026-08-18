##############################
#
# Scale Pattern Classification — GBIF Data
# Ian Becker
# May 2026
#
##############################

# Classifies each species x covariate combination as:
#   1. Consistent       — same direction at both scales (at least one sig)
#   2. Scale difference — opposite directions (at least one sig)
#   3. No clear effect  — non-significant at both scales

library(dplyr)

setwd("~/Desktop/project_code/UGS_occ/data")

output_dir <- "/Users/ianbecker/Desktop/project_code/UGS_occ/figures_tables"

####################
#   Load Data
####################

cat("Loading data...\n")

# Landscape level coefficients from GBIF models
all_results <- readRDS(list.files("model_results_gbif",
                                  pattern = "^all_results_gbif_.*\\.rds$",
                                  full.names = TRUE)[1])

# Within-site coefficients from GBIF models
within_coefs <- read.csv("within_site_models_gbif/within_site_results_gbif.csv")

cat("Landscape models loaded:", sum(sapply(all_results, function(x) x$success)), "\n")
cat("Within-site coefficients:", nrow(within_coefs), "\n\n")

####################
#   Extract Landscape Coefficients
####################

cat("Extracting landscape coefficients...\n")

land_coefs <- data.frame()

for (nm in names(all_results)) {
  if (!all_results[[nm]]$success) next
  
  model <- all_results[[nm]]$model
  beta  <- as.matrix(model$beta.samples)
  
  for (p in colnames(beta)) {
    vals <- beta[, p]
    land_coefs <- rbind(land_coefs, data.frame(
      species    = all_results[[nm]]$species,
      parameter  = p,
      land_mean  = mean(vals),
      land_lower = quantile(vals, 0.025),
      land_upper = quantile(vals, 0.975),
      stringsAsFactors = FALSE
    ))
  }
}

land_coefs <- land_coefs %>%
  filter(parameter != "(Intercept)") %>%
  mutate(
    parameter = gsub("[0-9]+$", "", parameter),
    parameter = gsub("_pct$", "", parameter),
    parameter = gsub("flooded_veg", "flooded_veg", parameter),
    land_sig  = land_lower > 0 | land_upper < 0,
    land_dir  = ifelse(land_mean > 0, "positive", "negative")
  )

cat("Landscape coefficients extracted:", nrow(land_coefs), "\n\n")

####################
#   Summarise Within-Site Coefficients
####################

cat("Summarising within-site coefficients...\n")

within_summary <- within_coefs %>%
  filter(parameter != "(Intercept)") %>%
  mutate(
    parameter = gsub("[0-9]+$", "", parameter),
    parameter = gsub("_pct$", "", parameter),
    site_sig  = lower > 0 | upper < 0,
    site_dir  = ifelse(Mean > 0, "positive", "negative")
  ) %>%
  group_by(species, parameter) %>%
  summarise(
    site_mean    = mean(Mean,     na.rm = TRUE),
    n_sites      = n(),
    n_sig_pos    = sum(site_sig & site_dir == "positive"),
    n_sig_neg    = sum(site_sig & site_dir == "negative"),
    n_sig        = sum(site_sig),
    n_pos        = sum(site_dir == "positive"),
    n_neg        = sum(site_dir == "negative"),
    # Majority direction
    site_dir_maj = ifelse(n_pos >= n_neg, "positive", "negative"),
    # Significant if majority of site combinations are significant
    site_sig_any = n_sig > 0,
    .groups = "drop"
  )

cat("Within-site summaries:", nrow(within_summary), "\n\n")

####################
#   Join and Classify
####################

cat("Classifying scale patterns...\n")

combined <- land_coefs %>%
  inner_join(within_summary, by = c("species", "parameter")) %>%
  mutate(
    # At least one scale significant
    any_sig = land_sig | site_sig_any,
    
    pattern = case_when(
      # No clear effect — both non-significant
      !land_sig & !site_sig_any                               ~ "no_clear_effect",
      # Same direction (at least one sig)
      any_sig & land_dir == site_dir_maj                      ~ "consistent",
      # Opposite directions (at least one sig)
      any_sig & land_dir != site_dir_maj                      ~ "scale_difference",
      TRUE                                                     ~ "no_clear_effect"
    ),
    
    pattern_label = case_when(
      pattern == "consistent"      ~ "Consistent",
      pattern == "scale_difference" ~ "Scale Difference",
      pattern == "no_clear_effect" ~ "No Clear Effect"
    )
  )

cat("Classified", nrow(combined), "species-covariate combinations\n\n")

####################
#   Summary
####################

cat("=== OVERALL PATTERN SUMMARY ===\n\n")

pattern_summary <- combined %>%
  count(pattern_label) %>%
  mutate(pct = round(n / sum(n) * 100, 1)) %>%
  arrange(pattern_label)

print(pattern_summary)

cat("\n=== BY COVARIATE ===\n\n")

cov_summary <- combined %>%
  group_by(parameter, pattern_label) %>%
  summarise(n = n(), .groups = "drop") %>%
  arrange(parameter, pattern_label)

print(cov_summary)

cat("\n=== SCALE DIFFERENCES ===\n\n")

scale_diffs <- combined %>%
  filter(pattern == "scale_difference") %>%
  select(species, parameter, land_mean, land_dir, land_sig,
         site_mean, site_dir_maj, site_sig_any, n_sites) %>%
  arrange(parameter, species)

cat(nrow(scale_diffs), "scale differences:\n")
print(scale_diffs)

cat("\n=== CONSISTENT EFFECTS ===\n\n")

consistent <- combined %>%
  filter(pattern == "consistent") %>%
  select(species, parameter, land_mean, land_dir, land_sig,
         site_mean, site_dir_maj, site_sig_any, n_sites) %>%
  arrange(parameter, species)

cat(nrow(consistent), "consistent effects:\n")
print(consistent)

####################
#   Counts for Results Section
####################

cat("\n=== COUNTS FOR RESULTS SECTION ===\n\n")

n_species_land   <- length(unique(land_coefs$species))
n_combos         <- nrow(combined)
covariates       <- unique(combined$parameter)

cat("Landscape species:", n_species_land, "\n")
cat("Species-covariate combinations:", n_combos, "\n\n")

for (cov in sort(covariates)) {
  cov_df <- combined %>% filter(parameter == cov)
  
  n_consistent_pos <- sum(cov_df$pattern == "consistent" &
                            cov_df$land_dir == "positive")
  n_consistent_neg <- sum(cov_df$pattern == "consistent" &
                            cov_df$land_dir == "negative")
  n_scale_diff     <- sum(cov_df$pattern == "scale_difference")
  n_no_effect      <- sum(cov_df$pattern == "no_clear_effect")
  n_land_sig_pos   <- sum(cov_df$land_sig & cov_df$land_dir == "positive")
  n_land_sig_neg   <- sum(cov_df$land_sig & cov_df$land_dir == "negative")
  
  cat(cov, ":\n")
  cat("  Landscape sig positive:", n_land_sig_pos, "/", nrow(cov_df), "\n")
  cat("  Landscape sig negative:", n_land_sig_neg, "/", nrow(cov_df), "\n")
  cat("  Consistent positive:   ", n_consistent_pos, "\n")
  cat("  Consistent negative:   ", n_consistent_neg, "\n")
  cat("  Scale difference:      ", n_scale_diff, "\n")
  cat("  No clear effect:       ", n_no_effect, "\n\n")
}

####################
#   Save
####################

write.csv(combined %>%
            select(species, parameter,
                   land_mean, land_dir, land_sig,
                   site_mean, site_dir_maj, site_sig_any,
                   n_sites, n_sig_pos, n_sig_neg,
                   pattern, pattern_label) %>%
            arrange(pattern, parameter, species),
          file.path(output_dir, "scale_pattern_classification_gbif.csv"),
          row.names = FALSE)
cat("Saved: scale_pattern_classification_gbif.csv\n")

write.csv(pattern_summary,
          file.path(output_dir, "scale_pattern_summary_gbif.csv"),
          row.names = FALSE)
cat("Saved: scale_pattern_summary_gbif.csv\n")

cat("\n=== COMPLETE ===\n")