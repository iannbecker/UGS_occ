##############################
#
# eBird Landscape vs Within-Site GBIF
# Ian Becker
# August 2026
#
##############################

# This script compares landscape-level results from ebird to
# within-site results from iNat

library(dplyr)
library(tidyr)
library(ggplot2)

setwd("~/Desktop/project_code/UGS_occ/data")

output_dir <- "/Users/ianbecker/Desktop/project_code/UGS_occ/figures_tables"

# ============================================================================
# 1. EXTRACT EBIRD LANDSCAPE COEFFICIENTS
# ============================================================================

# Load in eBird landscape-level results

ebird_results <- readRDS(list.files("model_results_ebird_checklist",
                                    pattern = "^all_results_ebird_checklist_.*\\.rds$",
                                    full.names = TRUE)[1])

ebird_coefs <- data.frame()

# Extract landscape-level results by species

for (nm in names(ebird_results)) {
  if (!ebird_results[[nm]]$success) next
  beta <- as.matrix(ebird_results[[nm]]$model$beta.samples)
  
  for (p in colnames(beta)) {
    p_clean <- gsub("[0-9]+$", "", p)
    if (p_clean == "(Intercept)") next
    vals <- beta[, p]
    
    ebird_coefs <- rbind(ebird_coefs, data.frame(
      species   = ebird_results[[nm]]$species,
      parameter = p_clean,
      mean      = mean(vals),
      lower     = quantile(vals, 0.025),
      upper     = quantile(vals, 0.975),
      sig       = quantile(vals, 0.025) > 0 | quantile(vals, 0.975) < 0,
      dir       = ifelse(mean(vals) > 0, "positive", "negative"),
      source    = "eBird Landscape",
      stringsAsFactors = FALSE
    ))
  }
}


# ============================================================================
# 2. EXTRACT WITHIN-SITE AVERAGED POSTERIORS INAT
# ============================================================================

# Get individual inat model files

model_files <- list.files("within_site_models_gbif",
                          pattern = "\\.rds$", full.names = TRUE)
model_files <- model_files[!grepl("all_results", model_files)]

# Extract wihtin site averaged posteriors for each species and covariate

within_avg <- list()

for (f in model_files) {
  tryCatch({
    result <- readRDS(f)
    sp     <- result$species
    beta   <- as.matrix(result$model$beta.samples)
    
    if (is.null(within_avg[[sp]])) within_avg[[sp]] <- list()
    
    for (p in colnames(beta)) {
      p_clean <- gsub("[0-9]+$", "", p)
      if (p_clean == "(Intercept)") next
      
      if (is.null(within_avg[[sp]][[p_clean]])) {
        within_avg[[sp]][[p_clean]] <- list(sum = beta[, p], n = 1)
      } else {
        within_avg[[sp]][[p_clean]]$sum <- within_avg[[sp]][[p_clean]]$sum +
          beta[, p]
        within_avg[[sp]][[p_clean]]$n   <- within_avg[[sp]][[p_clean]]$n + 1
      }
    }
  }, error = function(e) NULL)
}

# Compute mean and credible intervals for each species and covariate

within_coefs <- data.frame()

for (sp in names(within_avg)) {
  for (p in names(within_avg[[sp]])) {
    vals <- within_avg[[sp]][[p]]$sum / within_avg[[sp]][[p]]$n
    
    within_coefs <- rbind(within_coefs, data.frame(
      species   = sp,
      parameter = p,
      mean      = mean(vals),
      lower     = quantile(vals, 0.025),
      upper     = quantile(vals, 0.975),
      sig       = quantile(vals, 0.025) > 0 | quantile(vals, 0.975) < 0,
      dir       = ifelse(mean(vals) > 0, "positive", "negative"),
      source    = "Within-Site",
      stringsAsFactors = FALSE
    ))
  }
}

# ============================================================================
# 3. FILTER TO SHARED SPECIES AND COMBINE RESULTS
# ============================================================================

# Filter to shared species

shared_species <- intersect(unique(ebird_coefs$species),
                            unique(within_coefs$species))

# Species in shared dataset

length(shared_species)

# helper function to clean names 

clean_params <- function(df) {
  df %>% mutate(parameter = gsub("_pct$", "", parameter))
}

# Cleanup parameter names

ebird_coefs  <- clean_params(ebird_coefs)  %>% filter(species %in% shared_species)
within_coefs <- clean_params(within_coefs) %>% filter(species %in% shared_species)


# Combine results

all_coefs <- bind_rows(ebird_coefs, within_coefs) %>%
  mutate(source = factor(source,
                         levels = c("eBird Landscape", "Within-Site")))

# ============================================================================
# 4. SCALE RESULTS CLASSIFICATION
# ============================================================================

# Classify species-covariate combinations based on significance and direction

comparison <- ebird_coefs %>%
  select(species, parameter, ebird_mean = mean, ebird_lower = lower,
         ebird_upper = upper, ebird_sig = sig, ebird_dir = dir) %>%
  left_join(
    within_coefs %>%
      select(species, parameter, within_mean = mean, within_lower = lower,
             within_upper = upper, within_sig = sig, within_dir = dir),
    by = c("species", "parameter")
  ) %>%
  filter(!is.na(within_dir)) %>%
  mutate(
    any_sig = ebird_sig | within_sig,
    direction_agree = ebird_dir == within_dir,
    pattern = case_when(
      !ebird_sig & !within_sig ~ "no_clear_effect",
      any_sig & ebird_dir == within_dir ~ "consistent",
      any_sig & ebird_dir != within_dir ~ "scale_difference",
      TRUE ~ "no_clear_effect"
    ),
    pattern_label = case_when(
      pattern == "consistent" ~ "Consistent",
      pattern == "scale_difference" ~ "Scale Difference",
      pattern == "no_clear_effect" ~ "No Clear Effect"
    )
  )

# ============================================================================
# 5. SUMMARY
# ============================================================================

# Overall summary across covariates

pattern_summary <- comparison %>%
  count(pattern_label) %>%
  mutate(pct = round(n / sum(n) * 100, 1)) %>%
  arrange(pattern_label)

print(pattern_summary)

# Summary by covariate

cov_summary <- comparison %>%
  group_by(parameter, pattern_label) %>%
  summarise(n = n(), .groups = "drop") %>%
  arrange(parameter, pattern_label)

print(cov_summary)

# Direction agreement by covariate

comparison %>%
  group_by(parameter) %>%
  summarise(
    n_species    = n(),
    pct_agree    = round(mean(direction_agree, na.rm = TRUE) * 100, 1),
    n_both_sig   = sum(ebird_sig & within_sig),
    n_ebird_only = sum(ebird_sig & !within_sig),
    n_within_only = sum(!ebird_sig & within_sig),
    .groups = "drop"
  ) %>%
  arrange(desc(pct_agree)) %>%
  print()

# Scale differences 

scale_diffs <- comparison %>%
  filter(pattern == "scale_difference") %>%
  select(species, parameter, ebird_mean, ebird_dir, ebird_sig,
         within_mean, within_dir, within_sig) %>%
  arrange(parameter, species)

print(scale_diffs)

# ============================================================================
# 6. SAVE
# ============================================================================

# Save comparison results

write.csv(comparison,
          file.path(output_dir, "ebird_within_site_comparison_table.csv"),
          row.names = FALSE)
