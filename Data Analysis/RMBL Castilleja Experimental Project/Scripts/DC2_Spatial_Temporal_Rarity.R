setwd("/Users/jargrett/Desktop/Castilleja/Data Analysis/RMBL Castilleja Experimental Project")

#----------Packages----------#
library(tidyverse)
library(magrittr)
library(plyr)
library(conflicted)
library(car)
library(lme4)
library(MuMIn)
library(emmeans)
library(indicspecies)
library(pwr)

conflict_prefer("summarize", "dplyr")
conflict_prefer("filter",    "dplyr")
conflict_prefer("mutate",    "dplyr")
conflict_prefer("arrange",   "dplyr")
conflict_prefer("rename",    "dplyr")
conflict_prefer("select",    "dplyr")
conflict_prefer("count",     "dplyr")

#----------Importing Data----------#
rac_cover <- read.csv("Processed Data/Subset.csv")
rac_cover <- as.data.frame(unclass(rac_cover),stringsAsFactors=TRUE)

species_info <- read.csv("Raw Data/EL Species List - EL.csv")
species_info <- as.data.frame(unclass(species_info),stringsAsFactors=TRUE)
cover        <- readRDS("Processed Data/Cleaned Cover.rds")
plot         <- readRDS("Processed Data/Plot Data.rds")

#----------Spatial Rarity (SR)----------#
# Proportional rank of mean cover across all control plots and years
# SR near 1 = rare, SR near 0 = dominant

spatial_rarity <- rac_cover %>%
  filter(removal == "Control") %>%
  group_by(code) %>%
  summarize(mean_cover_control = mean(percent_cover, na.rm = TRUE), .groups = "drop") %>%
  mutate(SR = 1 - percent_rank(mean_cover_control))

#----------Temporal Rarity (TR)----------#
# Proportion of years absent from control plots
# TR near 1 = temporally rare, TR near 0 = temporally persistent

total_years <- rac_cover %>%
  summarize(n = n_distinct(year)) %>%
  pull(n)

temp_rarity <- rac_cover %>%
  filter(removal == "Control") %>%
  group_by(code) %>%
  summarize(present = n_distinct(year), .groups = "drop") %>%
  mutate(TR = 1 - (present / total_years))

#----------Abundance Typology----------#
# Spatial:  SR >= 0.25 = Sparse,       SR < 0.25 = Common
# Temporal: TR >= 0.50 = Intermittent, TR < 0.50 = Persistent

rarity <- spatial_rarity %>%
  left_join(temp_rarity, by = "code") %>%
  mutate(
    spatial  = if_else(SR >= 0.25, "Sparse",        "Common"),
    temporal = if_else(TR >= 0.50, "Intermittent",  "Persistent"),
    class = case_when(
      spatial == "Common" & temporal == "Persistent"   ~ "CP",
      spatial == "Common" & temporal == "Intermittent" ~ "CI",
      spatial == "Sparse" & temporal == "Persistent"   ~ "SP",
      spatial == "Sparse" & temporal == "Intermittent" ~ "SI"
    )
  )

#----------Response Ratio (RR)----------#
# RR,i = (C_removal - C_control) / (C_removal + C_control)
# Pooled across all years and plots — one value per species
# Species absent from control plots excluded

response_ratio <- rac_cover %>%
  group_by(code, removal) %>%
  summarize(mean_cover = mean(percent_cover, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = removal, values_from = mean_cover) %>%
  filter(!is.na(Removal) & !is.na(Control)) %>%
  mutate(RR = (Removal - Control) / (Removal + Control))

#----------Merge RR with Rarity and Species Info----------#
rarity_response <- response_ratio %>%
  rename(mean_cover_control = Control,
         mean_cover_removal = Removal) %>%
  left_join(rarity, by = "code") %>%
  left_join(species_info %>% select(code, functional_group, life_history, growth_form), by = "code") %>%
  filter(!is.na(SR) & !is.na(TR))

#----------Sensitivity Check — Typology Cutoffs ± 0.10----------#
make_rarity <- function(sr_cut, tr_cut) {
  spatial_rarity %>%
    left_join(temp_rarity, by = "code") %>%
    mutate(
      spatial  = if_else(SR >= sr_cut, "Sparse",       "Common"),
      temporal = if_else(TR >= tr_cut, "Intermittent", "Persistent"),
      class = case_when(
        spatial == "Common" & temporal == "Persistent"   ~ "CP",
        spatial == "Common" & temporal == "Intermittent" ~ "CI",
        spatial == "Sparse" & temporal == "Persistent"   ~ "SP",
        spatial == "Sparse" & temporal == "Intermittent" ~ "SI"
      )
    )
}

rarity_low  <- make_rarity(0.15, 0.40)
rarity_high <- make_rarity(0.35, 0.60)

sensitivity_check <- function(rarity_df, label) {
  df   <- response_ratio %>%
    left_join(rarity_df %>% select(code, SR, TR, class), by = "code") %>%
    filter(!is.na(SR) & !is.na(TR))
  m.SR <- lm(RR ~ SR, data = df)
  m.TR <- lm(RR ~ TR, data = df)
  cat("\n---", label, "---\n")
  cat("SR p:", coef(summary(m.SR))["SR", "Pr(>|t|)"], "\n")
  cat("TR p:", coef(summary(m.TR))["TR", "Pr(>|t|)"], "\n")
}

sensitivity_check(rarity_low,  "Cutoffs - 0.10")
sensitivity_check(rarity_high, "Cutoffs + 0.10")

#----------Sensitivity Check — Exclude Transients----------#
# Transient = not present in year 1 of control plots

year_one <- rac_cover %>% filter(removal == "Control", year == min(year)) %>% distinct(code)
no_transients <- rarity_response %>% filter(code %in% year_one$code)

model.SR.nt <- lm(RR ~ SR, data = no_transients)
model.TR.nt <- lm(RR ~ TR, data = no_transients)
model.nt    <- lm(RR ~ SR * TR, data = no_transients)

summary(model.SR.nt); Anova(model.SR.nt, type = "II"); confint(model.SR.nt)
summary(model.TR.nt); Anova(model.TR.nt, type = "II"); confint(model.TR.nt)
summary(model.nt)

pwr.f2.test(u = 1, v = 32, f2 = 0.096 / (1 - 0.096), sig.level = 0.05)# current power 0.45
pwr.f2.test(u = 1, f2 = 0.096 / (1 - 0.096), sig.level = 0.05, power = 0.80) #calulating number of species needed

#----------Delta Occurrence (ΔOcc)----------#
# Proportion of plot-year combinations present per treatment
# ΔOcc = Removed - Present (positive = more common in removal)

delta_occ <- cover %>%
  filter(year != "0") %>%
  complete(plot, year, code, fill = list(cover = 0, removal = NA)) %>%
  group_by(plot) %>%
  mutate(removal = first(na.omit(removal))) %>%
  ungroup() %>%
  filter(!code %in% c("bare", "litter", "rock", "CASE")) %>%
  mutate(present = ifelse(cover > 0, 1, 0)) %>%
  group_by(code, removal) %>%
  summarise(occ_rate = sum(present) / n(), .groups = "drop") %>%
  filter(!is.na(removal)) %>%
  pivot_wider(names_from = removal, values_from = occ_rate, values_fill = 0) %>%
  mutate(delta_occ = Removed - Present)

delta_rarity <- delta_occ %>%
  left_join(rarity_response %>% select(code, SR, TR, RR, class), by = "code") %>%
  filter(!is.na(SR))

#----------Indicator Species Analysis----------#
# Mean cover per plot averaged across years to avoid pseudoreplication
# Permutations constrained within paired plot blocks

ind_pooled_mean <- cover %>%
  filter(year != "0", !code %in% c("bare", "litter", "rock", "CASE")) %>%
  group_by(plot, code) %>%
  summarise(cover = mean(cover, na.rm = TRUE), .groups = "drop")

pooled_matrix <- ind_pooled_mean %>%
  pivot_wider(names_from = code, values_from = cover, values_fill = 0) %>%
  column_to_rownames("plot")

ind_meta <- cover %>%
  filter(year != "0", !code %in% c("bare", "litter", "rock", "CASE")) %>%
  group_by(plot) %>%
  slice(1) %>%
  arrange(plot)

case.removal_pooled <- ind_meta %>% pull(removal)
case.pair_pooled    <- ind_meta %>% pull(pair)

pooled_inv <- multipatt(pooled_matrix, case.removal_pooled, func = "r.g",
                        control = how(blocks = case.pair_pooled, nperm = 9999))

summary(pooled_inv, alpha = 0.1)

pooled_ind_results <- as.data.frame(pooled_inv$sign) %>%
  rownames_to_column("code") %>%
  select(code, stat, p.value) %>%
  rename(indval = stat) %>%
  mutate(significant = p.value < 0.1) %>%
  arrange(p.value)

sig_species <- pooled_ind_results %>% filter(significant) %>% pull(code) %>% as.character()

indval_rarity <- pooled_ind_results %>%
  left_join(delta_rarity %>% select(code, SR, TR, RR, delta_occ, class), by = "code")

#----------Plot-Level Response Classification----------#
# Identifies species that gained or lost plots in one treatment relative to the other
# Combines colonizer and extirpation analyses into a single asymmetry-based classification
# Threshold: response in one treatment must be at least 1 plot AND at least 2x the response in the other

plot_changes <- cover %>%
  filter(year %in% c(1, 3), !code %in% c("bare", "litter", "rock", "CASE")) %>%
  mutate(present = ifelse(cover > 0, 1, 0)) %>%
  group_by(code, removal, year) %>%
  summarise(n_plots = sum(present), .groups = "drop") %>%
  pivot_wider(names_from = c(removal, year), values_from = n_plots, values_fill = 0) %>%
  rename(control_y1 = Present_1, control_y3 = Present_3,
         removal_y1 = Removed_1, removal_y3 = Removed_3) %>%
  mutate(
    change_control = control_y3 - control_y1,
    change_removal = removal_y3 - removal_y1,
    turnover_magnitude = change_control - change_removal,
    turnover_class = case_when(
      change_control >=  1 & change_control >=  2 * abs(change_removal) & change_removal <= 0 ~ "Control colonizer",
      change_removal >=  1 & change_removal >=  2 * abs(change_control) & change_control <= 0 ~ "Removal colonizer",
      change_control <= -1 & abs(change_control) >=  2 * abs(change_removal) & change_removal >= 0 ~ "Control extirpation",
      change_removal <= -1 & abs(change_removal) >=  2 * abs(change_control) & change_control >= 0 ~ "Removal extirpation",
      TRUE                                                                                          ~ "No clear response"
    )
  )


# Tag species in indval_rarity with their response type and plot-level changes
indval_rarity <- indval_rarity %>%
  left_join(
    plot_changes %>% select(code, change_control, change_removal,
                            turnover_magnitude, turnover_class),
    by = "code"
  )
#---------Summary Dataframe---------#
robinhood <- indval_rarity %>%
  left_join(
    species_info %>% select(code, family, genus, species,
                            functional_group, life_history, growth_form),
    by = "code"
  ) %>%
  select(code, family, genus, species, functional_group, life_history, growth_form,
         spatial_rarity = SR, temporal_rarity = TR, response_ratio = RR, dominance_class = class, occupancy_shift = delta_occ,
         indicator = significant,
         turnover_class, turnover_magnitude)
saveRDS(robinhood, "Processed Data/Robinhood Summary.rds")
write.csv(robinhood, "Processed Data/Robinhood Summary.csv", row.names = FALSE)

# 1 — Rarity and Response Ratio
robin_rr <- robinhood %>%
  select(code, spatial_rarity, temporal_rarity, response_ratio, dominance_class) %>%
  filter(!is.na(spatial_rarity) & !is.na(temporal_rarity) & !is.na(response_ratio)) %>%
  arrange(spatial_rarity)

saveRDS(robin_rr, "Processed Data/Species Response Ratio.rds")

# 2 — Delta Occurrence
robin_delta <- delta_occ %>%
  filter(!is.na(delta_occ)) %>%
  select(code, occ_present = Present, occ_removed = Removed, delta_occ) %>%
  arrange(delta_occ)

saveRDS(robin_delta, "Processed Data/Species Change in Occurance.rds")

# 3 — Indicator Species
robin_indicator <- indval_rarity %>%
  select(code, indval, p.value) %>%
  arrange(p.value)

saveRDS(robin_indicator,  "Processed Data/Indicator Species Mean.rds")

# 4 — Turnover (colonizers and extirpated species combined)
robin_turnover <- indval_rarity %>%
  filter(turnover_class != "No clear response") %>%
  select(code, change_control, change_removal, turnover_magnitude, turnover_class) %>%
  arrange(turnover_class, desc(abs(turnover_magnitude)))

saveRDS(robin_turnover, "Processed Data/Species Turnover.rds")
