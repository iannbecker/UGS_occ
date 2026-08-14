##############################
#
# Two-Way Coefficient Comparison
# eBird Landscape vs Within-Site GBIF
# Ian Becker
# August 2026
#
##############################

library(dplyr)
library(tidyr)
library(ggplot2)

setwd("~/Desktop/project_code/UGS_occ/data")

output_dir <- "/Users/ianbecker/Desktop/project_code/UGS_occ/figures_tables"

# ============================================================================
# 1. EXTRACT EBIRD LANDSCAPE COEFFICIENTS
# ============================================================================

cat("Extracting eBird landscape coefficients...\n")

ebird_results <- readRDS(list.files("model_results_ebird_checklist",
                                    pattern = "^all_results_ebird_checklist_.*\\.rds$",
                                    full.names = TRUE)[1])

ebird_coefs <- data.frame()

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

cat("eBird landscape species:", n_distinct(ebird_coefs$species), "\n\n")

# ============================================================================
# 2. EXTRACT WITHIN-SITE AVERAGED POSTERIORS (GBIF)
# ============================================================================

cat("Extracting within-site averaged posteriors...\n")

model_files <- list.files("within_site_models_gbif",
                          pattern = "\\.rds$", full.names = TRUE)
model_files <- model_files[!grepl("all_results", model_files)]

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

cat("Within-site species:", n_distinct(within_coefs$species), "\n\n")

# ============================================================================
# 3. FILTER TO SHARED SPECIES AND COMBINE
# ============================================================================

cat("Filtering to shared species...\n")

shared_species <- intersect(unique(ebird_coefs$species),
                            unique(within_coefs$species))

cat("Species in both datasets:", length(shared_species), "\n\n")

clean_params <- function(df) {
  df %>% mutate(parameter = gsub("_pct$", "", parameter))
}

ebird_coefs  <- clean_params(ebird_coefs)  %>% filter(species %in% shared_species)
within_coefs <- clean_params(within_coefs) %>% filter(species %in% shared_species)

all_coefs <- bind_rows(ebird_coefs, within_coefs) %>%
  mutate(source = factor(source,
                         levels = c("eBird Landscape", "Within-Site")))

# ============================================================================
# 4. CLASSIFICATION
# ============================================================================

cat("Classifying scale patterns...\n")

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
    any_sig         = ebird_sig | within_sig,
    direction_agree = ebird_dir == within_dir,
    pattern = case_when(
      !ebird_sig & !within_sig              ~ "no_clear_effect",
      any_sig & ebird_dir == within_dir     ~ "consistent",
      any_sig & ebird_dir != within_dir     ~ "scale_difference",
      TRUE                                  ~ "no_clear_effect"
    ),
    pattern_label = case_when(
      pattern == "consistent"       ~ "Consistent",
      pattern == "scale_difference" ~ "Scale Difference",
      pattern == "no_clear_effect"  ~ "No Clear Effect"
    )
  )

cat("Classified", nrow(comparison), "species-covariate combinations\n\n")

# ============================================================================
# 5. SUMMARY
# ============================================================================

cat("=== OVERALL PATTERN SUMMARY ===\n\n")

pattern_summary <- comparison %>%
  count(pattern_label) %>%
  mutate(pct = round(n / sum(n) * 100, 1)) %>%
  arrange(pattern_label)

print(pattern_summary)

cat("\n=== BY COVARIATE ===\n\n")

cov_summary <- comparison %>%
  group_by(parameter, pattern_label) %>%
  summarise(n = n(), .groups = "drop") %>%
  arrange(parameter, pattern_label)

print(cov_summary)

cat("\n=== DIRECTION AGREEMENT BY COVARIATE ===\n\n")

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

cat("\n=== SCALE DIFFERENCES ===\n\n")

scale_diffs <- comparison %>%
  filter(pattern == "scale_difference") %>%
  select(species, parameter, ebird_mean, ebird_dir, ebird_sig,
         within_mean, within_dir, within_sig) %>%
  arrange(parameter, species)

cat(nrow(scale_diffs), "scale differences:\n")
print(scale_diffs)

# ============================================================================
# 6. RESULTS NUMBERS
# ============================================================================

cat("\n=== RESULTS NUMBERS BY COVARIATE ===\n\n")

covariates <- unique(comparison$parameter)

for (cov in sort(covariates)) {
  
  cov_df     <- comparison %>% filter(parameter == cov)
  n_species  <- nrow(cov_df)
  
  # Landscape (eBird)
  land_sig_pos <- cov_df %>% filter(ebird_sig & ebird_dir == "positive")
  land_sig_neg <- cov_df %>% filter(ebird_sig & ebird_dir == "negative")
  
  # Within-site
  site_sig_pos <- cov_df %>% filter(within_sig & within_dir == "positive")
  site_sig_neg <- cov_df %>% filter(within_sig & within_dir == "negative")
  
  cat(cov, ":\n")
  cat("  Landscape sig positive:", nrow(land_sig_pos), "/", n_species, "\n")
  cat("  Landscape sig negative:", nrow(land_sig_neg), "/", n_species, "\n")
  if (nrow(land_sig_pos) > 0)
    cat("  Land β range positive: [",
        round(min(land_sig_pos$ebird_mean), 3), "to",
        round(max(land_sig_pos$ebird_mean), 3), "]\n")
  if (nrow(land_sig_neg) > 0)
    cat("  Land β range negative: [",
        round(min(land_sig_neg$ebird_mean), 3), "to",
        round(max(land_sig_neg$ebird_mean), 3), "]\n")
  cat("  Within-site sig positive:", nrow(site_sig_pos), "/", n_species, "\n")
  cat("  Within-site sig negative:", nrow(site_sig_neg), "/", n_species, "\n")
  if (nrow(site_sig_pos) > 0)
    cat("  Site β range positive: [",
        round(min(site_sig_pos$within_mean), 3), "to",
        round(max(site_sig_pos$within_mean), 3), "]\n")
  if (nrow(site_sig_neg) > 0)
    cat("  Site β range negative: [",
        round(min(site_sig_neg$within_mean), 3), "to",
        round(max(site_sig_neg$within_mean), 3), "]\n")
  cat("\n")
}

# ============================================================================
# 7. COEFFICIENT PLOT
# ============================================================================

cat("Generating comparison plot...\n")

cov_labels <- c(
  "grass"             = "Grass",
  "trees"             = "Trees",
  "shrub"             = "Shrub",
  "flooded_veg"       = "Flooded Veg",
  "water"             = "Water",
  "crops"             = "Crops",
  "habitat_diversity" = "Hab. Diversity",
  "log_area"          = "Area (log)"
)

plot_df <- all_coefs %>%
  filter(parameter %in% names(cov_labels)) %>%
  mutate(
    cov_label = recode(parameter, !!!cov_labels),
    cov_label = factor(cov_label, levels = cov_labels)
  )

source_colors <- c(
  "eBird Landscape" = "#9b2335",
  "Within-Site"     = "#1f78b4"
)

p_compare <- ggplot(plot_df,
                    aes(x = mean, y = species,
                        color = source, shape = sig)) +
  geom_vline(xintercept = 0, linetype = "dashed",
             color = "gray50", linewidth = 0.4) +
  geom_errorbarh(aes(xmin = lower, xmax = upper),
                 height = 0, linewidth = 0.35,
                 position = position_dodge(width = 0.7),
                 alpha = 0.5) +
  geom_point(size = 1.8,
             position = position_dodge(width = 0.7)) +
  scale_color_manual(values = source_colors, name = NULL) +
  scale_shape_manual(
    values = c("TRUE" = 16, "FALSE" = 1),
    name   = "Significant"
  ) +
  facet_wrap(~ cov_label, scales = "free_x", ncol = 4,
             strip.position = "bottom") +
  labs(x = "Coefficient (β)", y = NULL) +
  theme_classic() +
  theme(
    axis.text.y        = element_text(size = 5.5, face = "italic"),
    axis.text.x        = element_text(size = 7),
    strip.text         = element_text(size = 9, face = "bold"),
    strip.background   = element_blank(),
    strip.placement    = "outside",
    legend.position    = "bottom",
    panel.grid.major.x = element_line(color = "gray90", linewidth = 0.3),
    panel.spacing      = unit(0.8, "lines")
  )

ggsave(
  file.path(output_dir, "ebird_within_site_comparison.png"),
  p_compare, width = 16, height = 12, dpi = 300, bg = "white"
)
cat("Saved: ebird_within_site_comparison.png\n")

# ============================================================================
# 8. SAVE
# ============================================================================

write.csv(comparison,
          file.path(output_dir, "ebird_within_site_comparison_table.csv"),
          row.names = FALSE)
cat("Saved: ebird_within_site_comparison_table.csv\n")

cat("\n=== COMPLETE ===\n")