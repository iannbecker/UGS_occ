library(dplyr)
library(tidyr)

# Within-site summary by species
within <- read.csv("within_site_results.csv")
within$parameter <- gsub("[0-9]+$", "", within$parameter)

within_summary <- within %>%
  filter(parameter != "(Intercept)") %>%
  mutate(sig_pos = lower > 0,
         sig_neg = upper < 0) %>%
  group_by(species, guild, parameter) %>%
  summarise(
    n_sites = n(),
    within_mean = round(mean(Mean), 2),
    within_se = round(sd(Mean) / sqrt(n()), 2),
    n_sig_pos = sum(sig_pos),
    n_sig_neg = sum(sig_neg),
    .groups = "drop"
  ) %>%
  mutate(
    within_pattern = case_when(
      n_sig_pos > 0 & n_sig_neg == 0 ~ paste0("+sig (", n_sig_pos, "/", n_sites, ")"),
      n_sig_neg > 0 & n_sig_pos == 0 ~ paste0("-sig (", n_sig_neg, "/", n_sites, ")"),
      n_sig_pos > 0 & n_sig_neg > 0 ~ "mixed",
      TRUE ~ "NS"
    )
  )

print(within_summary, n = 101)

############

# Load site-level results
site_results <- readRDS("all_results_2026-01-20.rds")  # your date

# Extract site-level coefficients
site_coefs <- data.frame()

for (sp in names(site_results)) {
  if (!site_results[[sp]]$success) next
  
  model <- site_results[[sp]]$model
  beta <- as.matrix(model$beta.samples)
  
  for (p in colnames(beta)) {
    vals <- beta[, p]
    site_coefs <- rbind(site_coefs, data.frame(
      species = site_results[[sp]]$species,
      parameter = p,
      site_mean = round(mean(vals), 2),
      site_lower = round(quantile(vals, 0.025), 2),
      site_upper = round(quantile(vals, 0.975), 2),
      site_sig = ifelse(quantile(vals, 0.025) > 0, "+sig",
                        ifelse(quantile(vals, 0.975) < 0, "-sig", "NS"))
    ))
  }
}

# Now merge

# Remove _pct suffix from site-level parameters
site_coefs$parameter <- gsub("_pct", "", site_coefs$parameter)

# Now merge again
comparison <- within_summary %>%
  left_join(site_coefs, by = c("species", "parameter")) %>%
  filter(parameter != "(Intercept)") %>%
  mutate(
    sign_flip = sign(site_mean) != sign(within_mean),
    sig_change = (site_sig != "NS" & within_pattern == "NS") | 
      (site_sig == "NS" & within_pattern != "NS" & within_pattern != "mixed"),
    scale_pattern = case_when(
      sign_flip ~ "SIGN FLIP",
      sig_change & site_sig != "NS" ~ "Sig→NS",
      sig_change & site_sig == "NS" ~ "NS→Sig",
      TRUE ~ "Consistent"
    )
  ) %>%
  select(species, guild, parameter, 
         site_mean, site_sig, 
         within_mean, within_pattern, 
         scale_pattern)

# Check it worked
comparison %>% filter(!is.na(site_mean))

# Show the interesting ones
comparison %>% 
  filter(scale_pattern != "Consistent") %>%
  arrange(scale_pattern, species) %>%
  print(n = 101)

