##############################
#
# Results Summary Numbers + Scale Classification — GBIF Data
# Ian Becker
# May 2026
#
##############################

# For each species x covariate:
#   - Averages posterior samples across viable sites (within-site)
#   - Extracts landscape-level posteriors
#   - Classifies as consistent, scale difference, or no clear effect
#   - Prints all numbers needed for results paragraphs

library(dplyr)

setwd("~/Desktop/project_code/UGS_occ/data")

output_dir <- "/Users/ianbecker/Desktop/project_code/UGS_occ/figures_tables"

####################
#   Load Data
####################

cat("Loading data...\n")

# Landscape models
all_results <- readRDS(list.files("model_results_gbif",
                                  pattern = "^all_results_gbif_.*\\.rds$",
                                  full.names = TRUE)[1])

# Within-site model files
model_files <- list.files("within_site_models_gbif",
                          pattern = "\\.rds$", full.names = TRUE)
model_files <- model_files[!grepl("all_results", model_files)]

# Within-site summary for naive occupancy
within_summary_df <- read.csv("within_site_model_summary_gbif.csv")

n_species_land <- sum(sapply(all_results, function(x) x$success))
cat("Landscape species:", n_species_land, "\n")
cat("Within-site model files:", length(model_files), "\n")
cat("Within-site combinations:", nrow(within_summary_df), "\n\n")

####################
#   Mean Within-Site Occupancy
####################

cat("=== WITHIN-SITE SPACE USE ===\n\n")

cat("Total viable species-site combinations:", nrow(within_summary_df), "\n")

mean_occ <- mean(within_summary_df$naive_occ)
set.seed(42)
boot_means <- replicate(10000, mean(sample(within_summary_df$naive_occ,
                                           replace = TRUE)))
ci_lower <- quantile(boot_means, 0.025)
ci_upper <- quantile(boot_means, 0.975)

cat("Mean within-site occupancy:", round(mean_occ * 100, 1), "%\n")
cat("95% CI:", round(ci_lower * 100, 1), "-",
    round(ci_upper * 100, 1), "%\n\n")

####################
#   Extract Landscape Posteriors
####################

cat("Extracting landscape posteriors...\n")

land_posteriors <- list()  # species -> parameter -> vector of samples

for (nm in names(all_results)) {
  if (!all_results[[nm]]$success) next
  sp   <- all_results[[nm]]$species
  beta <- as.matrix(all_results[[nm]]$model$beta.samples)
  
  land_posteriors[[sp]] <- list()
  
  for (p in colnames(beta)) {
    p_clean <- gsub("[0-9]+$", "", p)
    if (p_clean == "(Intercept)") next
    land_posteriors[[sp]][[p_clean]] <- beta[, p]
  }
}

cat("Landscape posteriors extracted for", length(land_posteriors), "species\n\n")

####################
#   Average Within-Site Posteriors Across Sites Per Species
####################

cat("Averaging within-site posteriors across sites per species...\n")

within_avg_posteriors <- list()  # species -> parameter -> averaged vector

for (f in model_files) {
  tryCatch({
    result <- readRDS(f)
    sp     <- result$species
    model  <- result$model
    beta   <- as.matrix(model$beta.samples)
    
    if (is.null(within_avg_posteriors[[sp]])) {
      within_avg_posteriors[[sp]] <- list()
    }
    
    for (p in colnames(beta)) {
      p_clean <- gsub("[0-9]+$", "", p)
      if (p_clean == "(Intercept)") next
      
      if (is.null(within_avg_posteriors[[sp]][[p_clean]])) {
        # First site — initialize with this site's samples and count
        within_avg_posteriors[[sp]][[p_clean]] <- list(
          sum     = beta[, p],
          n_sites = 1
        )
      } else {
        # Add this site's samples to running sum
        within_avg_posteriors[[sp]][[p_clean]]$sum <- 
          within_avg_posteriors[[sp]][[p_clean]]$sum + beta[, p]
        within_avg_posteriors[[sp]][[p_clean]]$n_sites <- 
          within_avg_posteriors[[sp]][[p_clean]]$n_sites + 1
      }
    }
  }, error = function(e) {
    cat("Error with", f, ":", e$message, "\n")
  })
}

# Divide sum by n_sites to get averaged posterior
for (sp in names(within_avg_posteriors)) {
  for (p in names(within_avg_posteriors[[sp]])) {
    n <- within_avg_posteriors[[sp]][[p]]$n_sites
    s <- within_avg_posteriors[[sp]][[p]]$sum
    within_avg_posteriors[[sp]][[p]] <- s / n
  }
}

cat("Within-site averages computed for", length(within_avg_posteriors),
    "species\n\n")

####################
#   Summarise Into Dataframes
####################

# Landscape summary
land_summary <- data.frame()

for (sp in names(land_posteriors)) {
  for (p in names(land_posteriors[[sp]])) {
    vals <- land_posteriors[[sp]][[p]]
    land_summary <- rbind(land_summary, data.frame(
      species    = sp,
      parameter  = p,
      land_mean  = mean(vals),
      land_lower = quantile(vals, 0.025),
      land_upper = quantile(vals, 0.975),
      land_sig   = quantile(vals, 0.025) > 0 | quantile(vals, 0.975) < 0,
      land_dir   = ifelse(mean(vals) > 0, "positive", "negative"),
      stringsAsFactors = FALSE
    ))
  }
}

# Within-site averaged summary
within_summary <- data.frame()

for (sp in names(within_avg_posteriors)) {
  for (p in names(within_avg_posteriors[[sp]])) {
    vals <- within_avg_posteriors[[sp]][[p]]
    
    within_summary <- rbind(within_summary, data.frame(
      species     = sp,
      parameter   = p,
      site_mean   = mean(vals),
      site_lower  = quantile(vals, 0.025),
      site_upper  = quantile(vals, 0.975),
      site_sig    = quantile(vals, 0.025) > 0 | quantile(vals, 0.975) < 0,
      site_dir    = ifelse(mean(vals) > 0, "positive", "negative"),
      stringsAsFactors = FALSE
    ))
  }
}

# Get n_sites from model files more efficiently
site_counts <- data.frame()
for (f in model_files) {
  tryCatch({
    r <- readRDS(f)
    site_counts <- rbind(site_counts, data.frame(
      species = r$species,
      site_id = r$site_id
    ))
  }, error = function(e) NULL)
}

n_sites_per_species <- site_counts %>%
  group_by(species) %>%
  summarise(n_sites = n(), .groups = "drop")

within_summary <- within_summary %>%
  left_join(n_sites_per_species, by = "species")

cat("Landscape species-covariate rows:", nrow(land_summary), "\n")
cat("Within-site species-covariate rows:", nrow(within_summary), "\n\n")

####################
#   Classify Scale Patterns
####################

cat("Classifying scale patterns...\n")

combined <- land_summary %>%
  inner_join(within_summary, by = c("species", "parameter")) %>%
  mutate(
    any_sig = land_sig | site_sig,
    pattern = case_when(
      !land_sig & !site_sig                          ~ "no_clear_effect",
      any_sig & land_dir == site_dir                 ~ "consistent",
      any_sig & land_dir != site_dir                 ~ "scale_difference",
      TRUE                                           ~ "no_clear_effect"
    ),
    pattern_label = case_when(
      pattern == "consistent"       ~ "Consistent",
      pattern == "scale_difference" ~ "Scale Difference",
      pattern == "no_clear_effect"  ~ "No Clear Effect"
    )
  )

cat("Classified", nrow(combined), "species-covariate combinations\n\n")

####################
#   Scale Classification Summary
####################

cat("=== SCALE CLASSIFICATION SUMMARY ===\n\n")

pattern_summary <- combined %>%
  count(pattern_label) %>%
  mutate(pct = round(n / sum(n) * 100, 1))
print(pattern_summary)

cat("\nBy covariate:\n")
cov_summary <- combined %>%
  group_by(parameter, pattern_label) %>%
  summarise(n = n(), .groups = "drop") %>%
  arrange(parameter, pattern_label)
print(cov_summary)

cat("\nScale differences:\n")
print(combined %>%
        filter(pattern == "scale_difference") %>%
        select(species, parameter, land_mean, land_dir, land_sig,
               site_mean, site_dir, site_sig, n_sites) %>%
        arrange(parameter, species))

####################
#   Landscape Paragraph Numbers
####################

cat("\n=== LANDSCAPE PARAGRAPH NUMBERS ===\n\n")

covariates <- c("grass_pct", "log_area", "trees_pct", "flooded_veg_pct",
                "water_pct", "shrub_pct", "habitat_diversity", "crops_pct")

for (cov in covariates) {
  
  cov_df     <- land_summary %>% filter(parameter == cov)
  n_total    <- nrow(cov_df)
  sig_pos_df <- cov_df %>% filter(land_sig & land_dir == "positive")
  sig_neg_df <- cov_df %>% filter(land_sig & land_dir == "negative")
  
  cat(cov, ":\n")
  cat("  Total species:", n_total, "\n")
  cat("  Sig positive:", nrow(sig_pos_df), "/", n_total, "\n")
  cat("  Sig negative:", nrow(sig_neg_df), "/", n_total, "\n")
  
  if (nrow(sig_pos_df) > 0) {
    strongest <- sig_pos_df %>% arrange(desc(land_mean)) %>% slice(1)
    weakest   <- sig_pos_df %>% arrange(land_mean)       %>% slice(1)
    cat("  β range positive: [", round(min(sig_pos_df$land_mean), 3),
        "to", round(max(sig_pos_df$land_mean), 3), "]\n")
    cat("  Strongest:", strongest$species,
        "(β =", round(strongest$land_mean, 3),
        "; 95% CI =", round(strongest$land_lower, 3), ",",
        round(strongest$land_upper, 3), ")\n")
    cat("  Weakest:", weakest$species,
        "(β =", round(weakest$land_mean, 3),
        "; 95% CI =", round(weakest$land_lower, 3), ",",
        round(weakest$land_upper, 3), ")\n")
  }
  
  if (nrow(sig_neg_df) > 0) {
    strongest <- sig_neg_df %>% arrange(land_mean)       %>% slice(1)
    weakest   <- sig_neg_df %>% arrange(desc(land_mean)) %>% slice(1)
    cat("  β range negative: [", round(min(sig_neg_df$land_mean), 3),
        "to", round(max(sig_neg_df$land_mean), 3), "]\n")
    cat("  Strongest:", strongest$species,
        "(β =", round(strongest$land_mean, 3),
        "; 95% CI =", round(strongest$land_lower, 3), ",",
        round(strongest$land_upper, 3), ")\n")
    cat("  Weakest:", weakest$species,
        "(β =", round(weakest$land_mean, 3),
        "; 95% CI =", round(weakest$land_lower, 3), ",",
        round(weakest$land_upper, 3), ")\n")
  }
  
  cat("\n")
}

####################
#   Within-Site Paragraph Numbers
####################

cat("=== WITHIN-SITE PARAGRAPH NUMBERS ===\n\n")
cat("(Species-level, averaged across sites)\n\n")

within_covariates <- unique(within_summary$parameter)

for (cov in sort(within_covariates)) {
  
  cov_df     <- within_summary %>% filter(parameter == cov)
  n_species  <- nrow(cov_df)
  sig_pos_df <- cov_df %>% filter(site_sig & site_dir == "positive")
  sig_neg_df <- cov_df %>% filter(site_sig & site_dir == "negative")
  
  cat(cov, ":\n")
  cat("  Total species:", n_species, "\n")
  cat("  Sig positive:", nrow(sig_pos_df), "/", n_species, "\n")
  cat("  Sig negative:", nrow(sig_neg_df), "/", n_species, "\n")
  
  if (nrow(sig_pos_df) > 0) {
    cat("  β range positive: [", round(min(sig_pos_df$site_mean), 3),
        "to", round(max(sig_pos_df$site_mean), 3), "]\n")
  }
  
  if (nrow(sig_neg_df) > 0) {
    cat("  β range negative: [", round(min(sig_neg_df$site_mean), 3),
        "to", round(max(sig_neg_df$site_mean), 3), "]\n")
  }
  
  cat("\n")
}

####################
#   Tree Cover Cross-Scale
####################

cat("=== TREE COVER CROSS-SCALE ===\n\n")

tree_compare <- combined %>%
  filter(parameter == "trees_pct", site_sig, site_dir == "positive") %>%
  mutate(site_stronger = abs(site_mean) > abs(land_mean))

cat("Species with sig positive within-site tree effect:",
    nrow(tree_compare), "\n")
cat("Of those, site effect stronger than landscape:",
    sum(tree_compare$site_stronger, na.rm = TRUE), "\n\n")

####################
#   Save
####################

write.csv(combined %>%
            select(species, parameter,
                   land_mean, land_dir, land_sig,
                   site_mean, site_dir, site_sig,
                   n_sites, pattern, pattern_label) %>%
            arrange(pattern, parameter, species),
          file.path(output_dir, "scale_pattern_classification_gbif.csv"),
          row.names = FALSE)
cat("Saved: scale_pattern_classification_gbif.csv\n")

write.csv(land_summary,
          file.path(output_dir, "landscape_coefficients_gbif.csv"),
          row.names = FALSE)
cat("Saved: landscape_coefficients_gbif.csv\n")

write.csv(within_summary,
          file.path(output_dir, "within_site_averaged_coefficients_gbif.csv"),
          row.names = FALSE)
cat("Saved: within_site_averaged_coefficients_gbif.csv\n")

cat("\n=== COMPLETE ===\n")